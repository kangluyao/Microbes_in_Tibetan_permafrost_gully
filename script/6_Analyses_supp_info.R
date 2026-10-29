# ==============================================================================
# Thermokarst gully: supplementary analyses (reorganised script)
# ==============================================================================
# Purpose
#   Compares un-collapsed vs. collapsed (thermokarst gully) plots for
#   (1) dominant MAG genera, (2) functional-gene composition and abundance,
#   (3) soil environmental variables and their heterogeneity, (4) sequencing /
#   MAG-recruitment statistics, and (5) 16S vs. MAG community congruence.
#
# Contents
#   0. Setup (paths, packages, parameters, themes, helper functions)
#   1. MAG taxonomy and dominant genera
#   2. Functional genes
#   3. Environmental factors
#   4. Sequencing statistics, MAG recruitment and MAG quality
#   5. Procrustes analysis of 16S vs. MAG community structure
#   6. Export of combined GLMM diagnostics and session information
#
# Statistical framework (all former LMMs are now GLMMs fitted with glmmTMB)
#   Random effect : (1 | Gully_id) in every model.
#   Fixed effects : Group only (no covariates).
#   Tests         : Type-II Wald chi-square tests (car::Anova) per model term.
#   Effect size   : coefficient of GroupCollapsed on the link scale
#                   (log link -> log ratio Collapsed / Un-collapsed).
#
#   Response type                              Family (link)
#   -----------------------------------------  ---------------------------------
#   Relative abundance / TPM sums (>= 0)       Gamma(log) if all values > 0,
#                                              otherwise Tweedie(log)
#   Bray-Curtis dissimilarity, MAG mapping     Beta(logit)
#   rate, clay + silt (proportions)
#   Euclidean distance (env. heterogeneity)    Gamma(log)
#   Plant richness (counts)                    Poisson(log)
#   Soil pH                                    Gaussian(identity)
#   See section 0.6 and the notes in each section for the rationale.
#
# Inputs (edit the paths in section 0.1)
#   script/0_read_data.R, script/read_data.R, script/Functions.R
#   <wd_fun>/ko_9_levels.txt, <wd_fun>/CN_ko_input.csv
#   <wd_fun>/MAGs/{bin_abundance_coverm.txt, annotation.txt, MAGs_attributes.txt}
#   <project>/revision/data/sequencing_depth.csv
#
# Outputs (written to <project>/revision/result when save_outputs = TRUE)
#   *.pdf / *.png               figures
#   tables/*.csv                GLMM effects, BH-adjusted Group effects, ...
#   diagnostics/*.csv, *.pdf    convergence, DHARMa residual diagnostics
#
# Objects expected from the sourced scripts (not defined here)
#   metadata, mags_tax, mags_abun_tab, ko_tpm_table, ko_count_table, KO_L2_df,
#   otu_16s, mag_comp, PCoA_plot_fun(), similar_deter_fun(), main_theme
# ==============================================================================

rm(list = ls())  # start from a clean workspace

# 0. SETUP ---------------------------------------------------------------------

## 0.1 Paths ---------------------------------------------------------------------
project_dir <- "e:/thermokarst_gully" # EDIT: project root

# The sourced scripts use paths relative to the project root, so the working
# directory is set once here.
setwd(project_dir)

# Input directories (wd_trait and wd_16s are only used by the sourced scripts)
wd_trait <- file.path(project_dir, "data", "trait", "rdp")
wd_16s <- file.path(project_dir, "data", "16S", "rdp")
wd_fun <- file.path(project_dir, "data", "metagenome")
save.dir <- file.path(project_dir, "result") # legacy name, kept for sourced scripts

# Input files read directly in this script
ko_levels_file <- file.path(wd_fun, "ko_9_levels.txt")
sel_ko_file <- file.path(wd_fun, "CN_ko_input.csv")
mag_abundance_file <- file.path(wd_fun, "MAGs", "bin_abundance_coverm.txt")
mag_tax_file <- file.path(wd_fun, "MAGs", "annotation.txt")
mag_attr_file <- file.path(wd_fun, "MAGs", "MAGs_attributes.txt")
seq_depth_file <- file.path(project_dir, "revision", "data", "sequencing_depth.csv")

# Output locations
result_dir <- file.path(project_dir, "revision", "result")
table_dir <- file.path(result_dir, "tables")
diag_dir <- file.path(result_dir, "diagnostics")
save_outputs <- TRUE # FALSE: run the analyses without writing any file

## 0.2 Packages ------------------------------------------------------------------
# The attach order is the same as in the original script because it determines
# function masking (e.g. tidyverse verbs take precedence over tidytable).
# Some packages are not called directly in this script but are kept because the
# sourced scripts may rely on them.
pacman::p_load(
  phyloseq, ape, vegan, Biostrings, tidytable, tidyverse, 
  rstatix, ggrepel, EnhancedVolcano, edgeR, picante, readxl,
  ggpubr, gghalves, ggvenn, cowplot, aplot,
  glmmTMB, DHARMa
)

# Packages called with `::` only (installed, but deliberately not attached)
for (pkg in c("car", "ggpp", "usedist", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

## 0.3 Analysis parameters ---------------------------------------------------------
analysis_seed <- 123 # permutation tests, NMDS, jittered points
n_permutations <- 999
dharma_n_sim <- 1000 # simulations per model for DHARMa residuals

group_levels <- c("Un-collapsed", "Collapsed") # reference level first
group_colours <- c("Un-collapsed" = "#79ceb8", "Collapsed" = "#e95f5c")
effect_colours <- c(Negative = "#79ceb8", Neutral = "grey", Positive = "#e95f5c")

# Every GLMM has the form: response ~ Group + (1 | Gully_id)
glmm_fixed_terms <- "Group"

effect_axis_title <- "Effect size (log ratio)"

## 0.4 Plot themes -----------------------------------------------------------------
# The original script redefined `main_theme` three times with different font
# sizes. They are kept as two explicit themes so that every figure keeps the
# size it had in the original workflow.
main_theme <- theme_bw() +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(size = 0.5),
        strip.text       = element_text(colour = "black", size = 7),
        strip.background = element_rect(colour = "black", fill = "grey"),
        axis.title       = element_text(color = "black", size = 7),
        axis.ticks       = element_line(color = "black", linewidth = 0.5),
        axis.text.y      = element_text(colour = "black", size = 6),
        axis.text.x      = element_text(colour = "black", size = 6),
        legend.position  = "none")

## 0.5 General helper functions ------------------------------------------------------
#' Convert p-values to significance symbols; `ns_label` marks non-significance.
p_to_stars <- function(p, ns_label = "") {
  as.vector(unclass(symnum(
    p,
    corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", ns_label)
  )))
}

#' Text label for a p-value shown on a figure.
format_p_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ "P < 0.001",
    TRUE ~ sprintf("P = %.3f", p)
  )
}

#' Assign Group from the sample name ("_C" = un-collapsed, "_T" = collapsed).
group_from_sample <- function(sample_names) {
  factor(
    dplyr::case_when(
      grepl("_C", sample_names) ~ "Un-collapsed",
      grepl("_T", sample_names) ~ "Collapsed"
    ),
    levels = group_levels
  )
}

#' Convert a percentage (0-100) to a proportion (0-1); values already on the
#' 0-1 scale are returned unchanged.
to_proportion <- function(x) {
  if (max(x, na.rm = TRUE) > 1) x / 100 else x
}

#' Mean, SD, n and SE of one column (works on grouped data frames as well).
summarise_variable <- function(data, column) {
  data %>%
    summarise(
      mean_value = mean({{ column }}, na.rm = TRUE),
      sd_value = sd({{ column }}, na.rm = TRUE),
      n_value = dplyr::n(),
      se_value = sd_value / sqrt(n_value),
      .groups = "drop"
    )
}

#' Save a ggplot to `result_dir` (skipped when save_outputs is FALSE).
save_plot <- function(plot, filename, width, height, units = "in") {
  if (!save_outputs) {
    return(invisible(NULL))
  }
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(result_dir, filename), plot,
         width = width, height = height, units = units
  )
}

#' Write a table to `table_dir` (skipped when save_outputs is FALSE).
save_table <- function(x, filename) {
  if (!save_outputs) {
    return(invisible(NULL))
  }
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file.path(table_dir, filename))
}

## 0.6 GLMM helper functions ------------------------------------------------------------
# Design notes
# * Responses are analysed on their original scale. The original LMMs z-scored
#   the response, which is incompatible with Gamma/Beta/Tweedie likelihoods.
# * Models contain Group as the only fixed effect (no Time, Slope or MAP).
# * Gamma and Tweedie both use a log link, so the Group coefficient is a log
#   ratio under either family and effect sizes stay comparable across responses.

#' Family for non-negative continuous responses (relative abundance, TPM,
#' concentrations): Gamma(log) when all values are positive, otherwise
#' Tweedie(log), which also accommodates exact zeros.
family_nonneg_continuous <- function(y) {
  if (all(y > 0, na.rm = TRUE)) {
    Gamma(link = "log")
  } else {
    glmmTMB::tweedie(link = "log")
  }
}

#' Resolve the family of one response. `family` can be a family object, a named
#' list of family objects, or a function of the response vector.
resolve_family <- function(family, response_name, y) {
  if (inherits(family, "family")) {
    family
  } else if (is.function(family)) {
    family(y)
  } else if (is.list(family) && !is.null(family[[response_name]])) {
    family[[response_name]]
  } else {
    stop("No family defined for response '", response_name, "'.")
  }
}

#' Support of each response (range, zeros, integer-valued); use it to confirm
#' that the chosen family fits the data.
summarise_response_support <- function(data, responses) {
  map_dfr(set_names(responses), \(v) {
    y <- data[[v]]
    tibble(
      n = sum(!is.na(y)),
      n_unique = n_distinct(y),
      min = min(y, na.rm = TRUE),
      max = max(y, na.rm = TRUE),
      n_zero = sum(y == 0, na.rm = TRUE),
      n_negative = sum(y < 0, na.rm = TRUE),
      all_integer = all(y == round(y), na.rm = TRUE)
    )
  }, .id = "id")
}

#' Fit one GLMM per response with glmmTMB.
#' Responses with fewer than `min_unique` distinct values are skipped (NULL).
#' Failed fits return NULL. Warnings raised during fitting are stored in the
#' "fit_warnings" attribute and reported by check_glmm_convergence().
#' The response is copied to `.response` so that unusual column names (spaces,
#' brackets, etc.) never have to be parsed in a formula.
fit_glmm_batch <- function(responses, data,
                           fixed_terms = "Group",
                           random_term = "(1 | Gully_id)",
                           family, ziformula = ~0, min_unique = 1) {
  formula_template <- reformulate(c(fixed_terms, random_term), response = ".response")
  
  map(set_names(responses), function(v) {
    model_data <- data
    model_data$.response <- data[[v]]
    
    if (n_distinct(model_data$.response) < min_unique) {
      message(v, ": skipped (fewer than ", min_unique, " distinct values)")
      return(NULL)
    }
    
    fam <- resolve_family(family, v, model_data$.response)
    fit_warnings <- character()
    fit <- tryCatch(
      withCallingHandlers(
        glmmTMB(formula_template,
                data = model_data, family = fam, ziformula = ziformula
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message(v, ": model failed - ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(fit)) attr(fit, "fit_warnings") <- unique(fit_warnings)
    fit
  })
}

#' Wald chi-square tests (Type II) and Group-coefficient summary of all models,
#' returned as a matrix with one column per response (NA for failed models).
#' Row names follow the convention of the original script:
#'   <term>.mean (estimate), .z (Wald z), .se, .chisq and .P (Type-II test).
extract_glmm_effects <- function(fits) {
  one_fit <- function(fit) {
    tryCatch(
      {
        coefs <- summary(fit)$coefficients$cond
        wald <- car::Anova(fit, type = 2)
        c(
          set_names(coefs[, "Estimate"], paste0(rownames(coefs), ".mean")),
          set_names(coefs[, "z value"], paste0(rownames(coefs), ".z")),
          set_names(coefs[, "Std. Error"], paste0(rownames(coefs), ".se")),
          set_names(wald[, "Chisq"], paste0(rownames(wald), ".chisq")),
          set_names(wald[, "Pr(>Chisq)"], paste0(rownames(wald), ".P"))
        )
      },
      error = \(e) NULL
    )
  }
  
  results <- map(fits, \(fit) if (is.null(fit)) NULL else one_fit(fit))
  if (all(map_lgl(results, is.null))) stop("All models failed; nothing to extract.")
  
  template <- names(compact(results)[[1]])
  results <- map(results, \(r) {
    if (is.null(r)) set_names(rep(NA_real_, length(template)), template) else r
  })
  do.call(cbind, results)
}

#' Group coefficients with Benjamini-Hochberg (FDR) adjustment across all
#' models of a batch. `component` selects the conditional ("cond") or
#' zero-inflation ("zi") part of the model. Supplementary output only; the
#' figures use the unadjusted P-values, as in the original analysis.
extract_group_pvals <- function(fits, term_pattern = "^Group",
                                component = "cond") {
  map_dfr(compact(fits), \(fit) {
    coefs <- summary(fit)$coefficients[[component]]
    tibble(
      term = rownames(coefs),
      estimate = coefs[, "Estimate"],
      std.error = coefs[, "Std. Error"],
      p.value = coefs[, 4]
    )
  }, .id = "id") %>%
    filter(str_detect(term, term_pattern)) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    select(id, term, estimate, std.error, p.value, p.adj)
}

#' Convergence, family and random-effect checks for a batch of fitted models.
check_glmm_convergence <- function(fits) {
  imap_dfr(fits, \(fit, id) {
    if (is.null(fit)) {
      return(tibble(
        id = id, family = NA_character_, converged = FALSE,
        pd_hessian = FALSE, re_variance = NA_real_, singular_re = NA,
        fit_warnings = NA_character_
      ))
    }
    re_variance <- tryCatch(
      as.numeric(VarCorr(fit)$cond[[1]][1, 1]),
      error = \(e) NA_real_
    )
    tibble(
      id = id,
      family = paste0(family(fit)$family, " (", family(fit)$link, ")"),
      converged = fit$fit$convergence == 0,
      pd_hessian = isTRUE(fit$sdr$pdHess),
      re_variance = re_variance,
      singular_re = re_variance < 1e-6,
      fit_warnings = paste(attr(fit, "fit_warnings"), collapse = " | ")
    )
  })
}

#' DHARMa simulated-residual diagnostics for a batch of fitted models.
#' Returns the DHARMa objects (for plotting) and a table with the p-values of
#' the uniformity (KS), dispersion and zero-inflation tests. The zero-inflation
#' test is only informative for responses that can contain zeros.
run_dharma_diagnostics <- function(fits, n_sim = 1000) {
  dharma_res <- map(fits, simulateResiduals, n = n_sim)
  
  diag_table <- imap_dfr(dharma_res, \(res, v) {
    tibble(
      id           = v,
      ks_p         = testUniformity(res, plot = FALSE)$p.value,
      dispersion_p = testDispersion(res, plot = FALSE)$p.value,
      zeroinfl_p   = testZeroInflation(res, plot = FALSE)$p.value
    )
  })
  
  list(residuals = dharma_res, diagnostics = diag_table)
}

#' Save all DHARMa residual plots of a batch to a single PDF.
save_dharma_plots <- function(dharma_res, file) {
  pdf(file, width = 8, height = 6)
  iwalk(dharma_res, \(res, v) {
    tryCatch(plot(res, title = v), error = \(e) plot(res))
  })
  dev.off()
}

#' Fit, check and export one family of GLMMs.
#'
#' @param label Short name used in output file names.
#' @param responses Character vector of response columns in `data`.
#' @param data Data frame holding responses, fixed terms and Gully_id.
#' @param fixed_terms Fixed-effect terms.
#' @param family Family object, named list of families, or function(y).
#' @param drop_invalid Drop responses that are constant or contain NA, as the
#'   original `select(where(~ !any(is.na(.))))` did after scaling.
#' @param min_unique Skip responses with fewer distinct values (NA results).
#' @return List with fits, effects (matrix), group_effects (BH table),
#'   diagnostics (table) and dharma (residual objects).
fit_glmm_batch <- function(responses, data,
                           fixed_terms = "Group",
                           random_term = "(1 | Gully_id)",
                           family, ziformula = ~0) {
  map(set_names(responses), function(v) {
    formula_v <- reformulate(c(fixed_terms, random_term),
                             response = paste0("`", v, "`"))
    tryCatch(
      glmmTMB(formula_v, data = data, family = family, ziformula = ziformula),
      error = function(e) {
        message(v, ": model failed - ", conditionMessage(e))
        NULL
      }
    )
  })
}

## 0.7 Effect-size plot helper (used for pathways and C/N/S/other genes) ------------------
#' Axis limit of an effect-size plot: the original fixed limit (`min_limit`),
#' expanded in steps of 0.5 when log-ratio effects need more room.
effect_axis_limit <- function(effect_table, min_limit) {
  extent <- abs(effect_table$GroupCollapsed.mean) + effect_table$GroupCollapsed.se
  if (all(is.na(extent))) {
    return(min_limit)
  }
  max(min_limit, ceiling(max(extent, na.rm = TRUE) * 1.3 * 2) / 2)
}

#' Effect of collapse (GroupCollapsed coefficient +/- SE) per variable.
#' @param effects Matrix from extract_glmm_effects().
#' @param variable_order Variables to plot, in axis order.
#' @param bands Data frame (xmin, xmax, fill) of background shading.
#' @param style "point" (flipped, labels at the top) or "bar" (horizontal axis).
#' @param label_scale,label_add Position of the significance symbol:
#'   sign(effect) * ((|effect| + SE) * label_scale + label_add).
plot_effect_sizes <- function(effects, variable_order, bands,
                              style = c("point", "bar"),
                              min_limit = 2, band_alpha = 0.1,
                              label_scale = 1.2, label_add = 0,
                              point_size = 2, hline_width = 0.375,
                              errorbar_width = 0.9, ns_label = "",
                              theme_base = main_theme) {
  style <- match.arg(style)
  
  effect_table <- effects %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("variables") %>%
    filter(variables %in% variable_order) %>%
    mutate(
      sig = p_to_stars(Group.P, ns_label),
      variables = factor(variables, levels = variable_order),
      colour = case_when(
        GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
        GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positive",
        Group.P > 0.05 ~ "Neutral"
      ),
      label_y = sign(GroupCollapsed.mean) *
        ((abs(GroupCollapsed.mean) + GroupCollapsed.se) * label_scale + label_add)
    )
  y_limit <- effect_axis_limit(effect_table, min_limit)
  
  if (style == "bar") {
    plot <- ggplot(effect_table, aes(variables, GroupCollapsed.mean, fill = colour)) +
      geom_hline(yintercept = 0, linewidth = hline_width, colour = "gray2") +
      geom_bar(stat = "identity", width = 0.6)
  } else {
    plot <- ggplot(effect_table, aes(variables, GroupCollapsed.mean, colour = colour)) +
      geom_hline(yintercept = 0, linewidth = hline_width, colour = "gray2") +
      geom_point(size = point_size)
  }
  
  plot <- plot +
    geom_errorbar(
      aes(
        ymin = GroupCollapsed.mean - GroupCollapsed.se,
        ymax = GroupCollapsed.mean + GroupCollapsed.se
      ),
      width = 0, position = position_dodge(width = 0.7), linewidth = errorbar_width
    ) +
    geom_text(aes(label = sig, y = label_y),
              position = position_dodge(0.1), vjust = 0.55
    ) +
    labs(x = NULL, y = effect_axis_title) +
    scale_y_continuous(expand = c(0, 0), limits = c(-y_limit, y_limit)) +
    pmap(bands, \(xmin, xmax, fill) {
      annotate("rect",
               xmin = xmin, xmax = xmax, ymin = -y_limit, ymax = y_limit,
               alpha = band_alpha, fill = fill
      )
    }) +
    theme_base
  
  if (style == "bar") {
    plot +
      scale_fill_manual(values = effect_colours) +
      scale_x_discrete(position = "bottom") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
  } else {
    plot +
      scale_colour_manual(values = effect_colours) +
      coord_flip() +
      scale_x_discrete(position = "top") +
      theme(legend.position = "none")
  }
}

#' Background bands from consecutive axis edges and one colour per band.
make_bands <- function(edges, fills) {
  tibble(xmin = head(edges, -1), xmax = tail(edges, -1), fill = fills)
}

# Extract the pairwise distances of the same gully and group from a distance matrix.
similar_deter_fun <- function(dist) {
  vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
            'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
  similar_data <- lapply(vars, function(x) 
    usedist::dist_subset(dist, 
                         grep(x, metadata$Sample_name, value = TRUE))) %>%
    do.call(cbind, .) %>%
    data.frame() %>%
    gather("tem_group", "distance") %>%
    cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
          Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
    select(-tem_group) %>%
    mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
           Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))
}

# 1. MAG TAXONOMY AND DOMINANT GENERA --------------------------------------------------
source("script/0.read_data.R") # provides metadata, mags_tax, mags_abun_tab, ...
if ("package:tidytable" %in% search()) detach("package:tidytable", unload = FALSE)

## 1.1 Number of MAGs per phylum (pie chart) ------------------------------------------
n_top_phyla <- 11 # phyla shown individually; the rest is pooled as "Others"

mag_num_in_phylum <- data.frame(table(mags_tax[, "Phylum"])) %>%
  arrange(desc(Freq))
mag_num_in_phylum <- rbind(
  mag_num_in_phylum[1:n_top_phyla, ],
  data.frame(
    Var1 = "Others",
    Freq = sum(mag_num_in_phylum[-(1:n_top_phyla), 2])
  )
)

mag_count_data <- data.frame(
  Phylum = mag_num_in_phylum$Var1,
  mag_num = mag_num_in_phylum$Freq,
  prop = mag_num_in_phylum$Freq / sum(mag_num_in_phylum$Freq) * 100
) %>%
  mutate(
    lab.ypos = cumsum(prop) - 0.5 * prop, # label position within each slice
    Phylum = factor(Phylum, ordered = TRUE, levels = Phylum)
  )

phylum_colours <- c(
  "#89c5da", "#ffc15c", "#74d944", "#CE50CA", "#5e738f", "#C0717C", "#CBD5ec",
  "#5F7FC7", "#00718b", "#00b0f6", "#a3a500", "#f8766d", "#D7C1B1", "#689030",
  "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248",
  "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"
)

pie_for_mag_num_phylum <- ggplot(
  mag_count_data,
  aes(x = "", y = prop, fill = reorder(Phylum, -lab.ypos))
) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(
    aes(x = 1.35, y = lab.ypos, label = paste0(mag_num, " (", round(prop, 1), "%)")),
    colour = "black", size = 3
  ) +
  scale_fill_manual("Phylum", values = phylum_colours) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_void() +
  theme(legend.position = "left", legend.key.size = unit(0.25, "cm"))
pie_for_mag_num_phylum

## 1.2 Number of MAGs per domain and per order -----------------------------------------
# Original run: Bacteria = 306, Archaea = 24
mags_num_domain <- mags_tax %>%
  count(Domain, name = "number") %>%
  arrange(desc(number)) %>%
  mutate(
    Domain = gsub("d__", "", Domain),
    Domain = factor(Domain, levels = Domain)
  )
mags_num_domain

# Top 15 orders (MAGs without order assignment are excluded)
mags_num_order_plot <- mags_tax %>%
  count(Order, name = "number") %>%
  arrange(desc(number)) %>%
  mutate(
    Order = gsub("o__", "", Order),
    Order = factor(Order, levels = Order)
  ) %>%
  filter(!Order %in% "") %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = Order, y = number)) +
  geom_bar(stat = "identity", fill = "#ffc15c", width = 0.65) +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  labs(x = "Order", y = "Number of MAGs") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 7),
    axis.text.x = element_text(
      colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = 6
    ),
    axis.text.y = element_text(colour = "black", size = 6)
  )
mags_num_order_plot

## 1.3 Ten most abundant genera -------------------------------------------------------
# Relative abundance of each MAG with its genus assignment; sample columns are
# the ones whose names start with "G" (G1_C..., G1_T..., ...).
mags_genus_abun <- mags_abun_tab %>%
  rownames_to_column("Genome_id") %>%
  left_join(mags_tax, by = "Genome_id") %>%
  select("Genus", starts_with("G")) %>%
  mutate(Genus = gsub("g__", "", Genus))

# Ranking criterion (as in the original): mean abundance over all MAGs of the
# genus and all samples; MAGs without genus assignment are excluded.
top10_mags_id <- mags_genus_abun %>%
  group_by(Genus) %>%
  summarise(row_mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  filter(!Genus == "") %>%
  slice_max(row_mean, n = 10) %>%
  pull(Genus)

## 1.4 GLMM: effect of collapse on genus abundance ------------------------------------
# Response : genus relative abundance, summed over the MAGs of the genus (%).
# GLMM     : covariates removed, i.e. abundance ~ Group + (1 | Gully_id).
# GLMM     : glmmTMB, Gamma(log) if all values > 0, otherwise Tweedie(log), so
#            that zeros (genus absent from a sample) are allowed.
top10_mags_meta <- mags_genus_abun %>%
  filter(Genus %in% top10_mags_id) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample_name") %>%
  left_join(metadata, by = "Sample_name")

top10_mags_model_data <- top10_mags_meta %>%
  select(all_of(c("Group", "Gully_id", top10_mags_id))) %>%
  mutate(Group = factor(Group, levels = group_levels))

top10_glmm <- fit_glmm_batch(top10_mags_id, top10_mags_model_data,
                           family = tweedie(link = "log"))

top10_mags_pvals <- extract_group_pvals(top10_glmm)
print(top10_mags_pvals)

top10_mags_diag <- run_dharma_diagnostics(top10_glmm)
print(top10_mags_diag$diagnostics)

## 1.5 Box plot of the top-10 genera ---------------------------------------------------
# Significance symbols from the Group effect (unadjusted P, "n.s" if > 0.05)
top10_sig_labels <- top10_mags_pvals %>%
  as.data.frame() %>%
  mutate(sig = p_to_stars(p.adj, ns_label = "n.s"), y1 = 1.85) %>%
  as_tibble()

# NOTE: the box plot shows individual MAGs (one value per MAG and sample),
# whereas the GLMMs above use genus totals per sample.
top10_mags_abun_df <- mags_genus_abun %>%
  select(-Genome_id) %>%
  filter(Genus %in% top10_mags_id) %>%
  pivot_longer(-Genus, names_to = "Sample_name", values_to = "Rela_abun") %>%
  mutate(
    Group = group_from_sample(Sample_name),
    Genus = factor(Genus, levels = top10_mags_id)
  )

top10_mags_box_plot <- ggplot(top10_mags_abun_df, aes(Genus, Rela_abun)) +
  geom_boxplot(width = 0.5, aes(fill = Group), outlier.colour = NA, lwd = 0.5) +
  scale_y_continuous(limits = c(0, 2)) +
  scale_fill_manual(values = group_colours) +
  geom_text(
    data = top10_sig_labels, aes(x = id, y = y1, label = sig),
    size = 3, inherit.aes = FALSE
  ) +
  labs(x = "Genera", y = "Relative abundance (%)", fill = "Group") +
  theme_bw() +
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(
    axis.title = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    panel.grid = element_blank()
  )
top10_mags_box_plot

mags_num_plot <- cowplot::plot_grid(
  pie_for_mag_num_phylum, top10_mags_box_plot,
  nrow = 2, labels = c("a", "b"), label_size = 7,
  rel_widths = c(0.3, 0.7), align = "hv"
)

# save_plot(pie_for_mag_num_phylum, "pie_for_mag_num_phylum.pdf", 90, 55, "mm")
# save_plot(top10_mags_box_plot, "top10_mags_box_plot.pdf", 115, 55, "mm")


# 2. FUNCTIONAL GENES --------------------------------------------------------------------
## 2.1 Functional composition (Bray-Curtis on KO TPM) ------------------------------------
ko_dist <- vegdist(t(ko_tpm_table), method = "bray")

# PERMANOVA, MRPP and ANOSIM of functional composition between groups
# (permutations are not restricted by Gully_id, as in the original analysis)
set.seed(analysis_seed)
perm_scheme <- how(blocks = metadata$Gully_id, nperm = 999)
perm_test <- adonis2(ko_dist ~ Group,
                     data = metadata,
                     permutations = perm_scheme)
print(perm_test)

# PCoA ordination -------------------------------------------------------
ord.fun <- cmdscale(ko_dist, k = 2, eig = TRUE, add = TRUE)
r2      <- perm_test$R2[1]
p_value <- perm_test$`Pr(>F)`[1]

# Sample scores and group centroids
scores <- data.frame(ord.fun$points) |>
  rownames_to_column("Sample_name") |>
  left_join(metadata[, c("Sample_name", "Group")], by = "Sample_name") |>
  dplyr::rename(Dim1 = X1, Dim2 = X2)

centroids <- aggregate(cbind(Dim1, Dim2) ~ Group, data = scores, FUN = mean)

# Combine sample scores with their group centroids
plot_data <- merge(scores, centroids, by = "Group",
                   suffixes = c("", ".centroid")) |>
  mutate(Group = factor(Group, levels = group_levels))

## Plot ------------------------------------------------------------------
ko_pcoa_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2,
                                           color = Group,
                                           shape = Group)) +
  stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1,
               show.legend = FALSE, level = 0.95) +
  geom_point(size = 1.5) +
  scale_color_manual(name = "Group",
                     values = group_colours,
                     breaks = group_levels) +
  scale_shape_manual(values = c("Un-collapsed" = 16, "Collapsed" = 17)) +
  labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig),
                                   digits = 3), "%)", sep = ""),
       y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig),
                                   digits = 3), "%)", sep = "")) +
  annotate("text", x = Inf, y = Inf,
           label = paste("PERMANOVA R2 =", round(r2, 2),
                         "P =", round(p_value, 3)),
           hjust = 1.05, vjust = 1.1, size = 2) +
  main_theme
ko_pcoa_plot

## 2.2 Within-group dissimilarity of functional composition -------------------------------
# Response : Bray-Curtis dissimilarity between samples of the same gully and
#            group (output of similar_deter_fun(), one row per sample pair).
# Original : lmm_fun("distance", data); the formula inside lmm_fun() is defined
#            in Functions.R, which was not available. The GLMM below is
#            distance ~ Group + (1 | Gully_id).
# GLMM     : glmmTMB, Beta(logit) because Bray-Curtis dissimilarity is bounded
#            in (0, 1); random effect (1 | Gully_id).
distance_ko_data <- similar_deter_fun(ko_dist)

distance_ko_df <- distance_ko_data %>%
  mutate(Group = factor(Group, levels = group_levels))

# Fit the GLMM
dis_ko_glmm <- fit_glmm_batch("distance", distance_ko_df,
                             family = beta_family(link = "logit"))

dis_ko_pvals <- extract_group_pvals(dis_ko_glmm)
print(dis_ko_pvals)

dis_ko_diag <- run_dharma_diagnostics(dis_ko_glmm)
print(dis_ko_diag$diagnostics)

# The P-value shown in the figure is taken from the GLMM
dis_ko_label <- format_p_label(dis_ko_pvals$p.adj)

dis_ko_plot <- distance_ko_df %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(
    position = position_nudge(x = 0.25), side = "r", width = 0.5,
    colour = NA, alpha = 0.65
  ) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.colour = NA, alpha = 0.65) +
  geom_jitter(
    aes(fill = Group, colour = Group),
    shape = 21, size = 0.5, alpha = 0.65,
    position = position_jitter(width = 0.15, seed = analysis_seed)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.5, npcy = 0.95, size = 2, label = dis_ko_label) +
  labs(x = "Group", y = "Dissimilarity in gene composition") +
  scale_fill_manual(values = group_colours) +
  scale_colour_manual(values = group_colours) +
  main_theme +
  theme(legend.position = "none")
dis_ko_plot

## 2.3 Genes detected in only one group (Venn diagram) -----------------------------------
# A KO is "detected" in a group if its summed TPM over the group's samples is > 0.
detected_kos <- function(tpm_table, sample_names) {
  rownames(tpm_table)[which(rowSums(tpm_table[, sample_names, drop = FALSE]) > 0)]
}
ko_sets <- list(
  `Un-collapsed` = detected_kos(ko_tpm_table, grep("_C", metadata$Sample_name, value = TRUE)),
  Collapsed = detected_kos(ko_tpm_table, grep("_T", metadata$Sample_name, value = TRUE))
)

ko.venn <- ggvenn(
  ko_sets,
  fill_color = c("#f8766d", "#a3a500", "#00b0f6"),
  stroke_color = NA,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
)
ko.venn

## 2.4 Relative abundance of the 15 most abundant metabolic pathways (KEGG L2) -----------
ko_9_levels <- data.table::fread(ko_levels_file, header = TRUE, stringsAsFactors = FALSE)

KO_L2_df <- ko_tpm_table %>%
  rownames_to_column("KO") %>%
  left_join(ko_9_levels, by = "KO") %>%
  filter(L1 == "Metabolism") %>%
  select(L2, starts_with("G")) %>%
  distinct()

# Alternative input kept from the original script (currently not used):
KO_L2_df <- read.table(file.path(wd_fun, "KEGG.PathwayL2.raw.txt"), header = TRUE,
                       sep = "\t", comment.char = "", check.names = FALSE) |>
  rename(L2 = PathwayL2)

# NOTE: KO_L2_df is not created in this script (see the commented block above);
# it must come from the sourced scripts. The L3 table built above, KO_L3_df,
# is used for the pathway GLMMs in section 2.6.
if (!exists("KO_L2_df")) {
  stop("KO_L2_df is not defined; load it in read_data.R or see the note above.")
}

aggre_top10_data <- KO_L2_df %>%
  filter(!is.na(L2)) %>%
  group_by(L2) %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(across(where(is.numeric), ~ . / sum(.) * 100)) %>%
  mutate(RowMean = rowMeans(across(where(is.numeric)))) %>%
  slice_max(RowMean, n = 10) %>%
  bind_rows(
    # pool all remaining pathways as "Others" so that each sample sums to 100
    summarise(., across(where(is.numeric) & !RowMean, ~ 100 - sum(.))) %>%
      mutate(RowMean = rowMeans(across(where(is.numeric))), L2 = "Others", .before = 1)
  ) %>%
  select(-RowMean)

pathway_order <- aggre_top10_data %>% pull(L2)

pathway_colours <- c(
  "#FF7F00", "#B3DE69", "#1F78B4", "#FB8072", "#FDB462", "#80B1D3", "#BEBADA",
  "#c988d1", "#FCCDE5", "#A6CEE3", "#8DD3C7", "#B2DF8A", "#33A02C", "#FB9A99",
  "#FFFFB3", "grey"
)

pathway_plot <- aggre_top10_data %>%
  pivot_longer(-L2, names_to = "Sample_id", values_to = "TPM") %>%
  mutate(
    L2 = factor(L2, ordered = TRUE, levels = pathway_order),
    Group = group_from_sample(Sample_id),
    # x position = sample index (column order of the input table)
    sample_id = match(Sample_id, unique(Sample_id))
  ) %>%
  group_by(Group, Sample_id) %>%
  mutate(prop = TPM * 100 / sum(TPM)) %>%
  ungroup() %>%
  ggplot(aes(x = sample_id, y = prop, fill = L2)) +
  geom_area() +
  labs(x = "Sample", y = "Proportion (%)") +
  scale_fill_manual(values = pathway_colours) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  main_theme +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "lines"),
    panel.spacing = unit(0, "lines")
  )
pathway_plot

## 2.5 Combined figure: Venn diagram, PCoA and pathway composition ----------------------
ko_all_plot <- ggdraw() +
  draw_plot(ko.venn, x = 0, y = 0.5, width = 1 / 3, height = 1 / 2) +
  draw_plot(ko_pcoa_plot, x = 0, y = 0, width = 1 / 3, height = 1 / 2) +
  draw_plot(pathway_plot, x = 1 / 3, y = 0, width = 2 / 3, height = 1) +
  draw_plot_label(
    label = c("a", "b", "c"), size = 8,
    x = c(0, 0, 1 / 3), y = c(1, 0.5, 1)
  )
ko_all_plot

## 2.6 GLMM: effect of collapse on metabolic pathways (KEGG L3) ---------------------------
# Response : summed TPM of the KOs of each L3 pathway (Metabolism only).
# Original : lmer(scaled TPM ~ Group + (1 | Gully_id)).
# GLMM     : glmmTMB, Gamma(log) if all values > 0, otherwise Tweedie(log);
#            random effect (1 | Gully_id). Pathways with fewer than 3 distinct
#            values keep an NA result so that axis positions are unchanged.
pathway_data <- KO_L2_df %>%
  filter(!is.na(L2)) %>%
  group_by(L2) %>%
  summarise(across(where(is.numeric), sum)) %>%
  filter(L2 != "") %>%
  column_to_rownames("L2") %>%
  t() %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("Sample_id") %>%
  mutate(
    Gully_id = sapply(stringr::str_split(Sample_id, "_", n = 2), `[`, 1),
    Group = group_from_sample(Sample_id)
  ) %>%
  select(-Sample_id)
pathway_responses <- setdiff(names(pathway_data), c("Gully_id", "Group"))

# Fit the GLMM
pathway_glmm <- fit_glmm_batch(pathway_responses, pathway_data,
                              family = Gamma(link = "log"))

pathway_pvals <- extract_group_pvals(pathway_glmm)
print(pathway_pvals)

pathway_diag <- run_dharma_diagnostics(pathway_glmm)
print(pathway_diag$diagnostics)

pathway_effects <- extract_glmm_effects(pathway_glmm)
pathway_effects[1:min(6, nrow(pathway_effects)), 1:min(6, ncol(pathway_effects))]

## 2.7 Effect-size plot for pathways -------------------------------------------------------
pathway_index <- unique(ko_9_levels[, c("L2")]) %>% pull(L2)
pathway_index <- pathway_index[pathway_index != ""]

# Shaded bands mark the KEGG level-2 categories (positions follow the order of
# pathway_index)
pathway_bands <- make_bands(
  edges = c(0, 12.5, 22.5, 27.5, 30.5, 35.5, 46.5),
  fills = c("#ff6666", "#f19837", "#e56eee", "#5fb236", "#035096", "#8b008b")
)

diff_in_pathway_plot <- plot_effect_sizes(
  pathway_effects,
  variable_order = rev(pathway_index),
  bands = pathway_bands,
  style = "point", min_limit = 0.8, band_alpha = 0.2,
  label_scale = 1, label_add = 0.2,
  point_size = 1.5, hline_width = 0.7, errorbar_width = 0.9
)
diff_in_pathway_plot

# save_plot(diff_in_pathway_plot, "diff_in_pathway_plot.pdf", 4.5, 16)

## 2.8 Genes of the carbon, nitrogen, sulfur and other element cycles ----------------------
sel_ko <- read_csv(sel_ko_file, col_names = TRUE, show_col_types = FALSE)

# Count table of the selected KOs (kept from the original; not used below)
C_N_ko_count_tab <- ko_count_table[rownames(ko_count_table) %in% sel_ko$KO, ]
C_N_ko_count_tab <- cbind(KO = rownames(C_N_ko_count_tab), C_N_ko_count_tab)

C_N_ko_tpm_tab <- ko_tpm_table[rownames(ko_tpm_table) %in% sel_ko$KO, ]
C_N_ko_tpm_tab <- cbind(KO = rownames(C_N_ko_tpm_tab), C_N_ko_tpm_tab)

# Sum KO TPM per enzyme / protein encoded
C_N_path1 <- C_N_ko_tpm_tab %>%
  inner_join(sel_ko, by = "KO") %>%
  select(-1) %>%
  group_by(Category, Subcategory, Enzyme_protein_encoded) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  data.frame() %>%
  mutate(Subcategory = factor(Subcategory, 
                              levels = unique(sel_ko$Subcategory), 
                              ordered = TRUE))

# Enzyme x sample matrix (rows = enzymes, columns = samples)
C_N_path2 <- C_N_path1 %>%
  select(-Category, -Subcategory) %>%
  column_to_rownames("Enzyme_protein_encoded")

# The sample columns of C_N_path2 are matched to metadata by position below
if (!identical(as.character(metadata$Sample_name), colnames(C_N_path2))) {
  stop("Sample order differs between metadata and the functional-gene table.")
}

## 2.9 GLMM: effect of collapse on functional-gene abundance -------------------------------
# Response : summed TPM per enzyme / protein.
# Original : lmer(scaled TPM ~ Group + Time + Slope + MAP + (1 | Gully_id)).
# GLMM     : covariates removed, i.e. TPM ~ Group + (1 | Gully_id).
# GLMM     : glmmTMB, Gamma(log) if all values > 0, otherwise Tweedie(log)
#            (rare genes frequently contain zeros); random effect (1 | Gully_id).
#            Constant responses are dropped and responses with fewer than 3
#            distinct values are skipped, as in the original.
functional_gene_data <- data.frame(
  metadata[, c("Gully_id", "Group")],
  t(C_N_path2),
  check.names = FALSE
) %>%
  mutate(Group = factor(Group, levels = group_levels))

functional_gene_glmm <- fit_glmm_batch(rownames(C_N_path2), functional_gene_data,
                                       family = ziGamma(link = "log"))

functional_gene_pvals <- extract_group_pvals(functional_gene_glmm)
print(functional_gene_pvals)

functional_gene_diag <- run_dharma_diagnostics(functional_gene_glmm)
print(functional_gene_diag$diagnostics)

functional_gene_effects <- extract_glmm_effects(functional_gene_glmm)
functional_gene_effects[1:min(6, nrow(functional_gene_effects)), 1:min(6, ncol(functional_gene_effects))]


functional_gene_glmm <- run_glmm_workflow(
  label = "functional_genes",
  responses = rownames(C_N_path2),
  data = functional_gene_data,
  fixed_terms = glmm_fixed_terms,
  family = family_nonneg_continuous,
  drop_invalid = TRUE,
  min_unique = 3
)
functional_gene_effects <- functional_gene_glmm$effects
functional_gene_effects[1:min(6, nrow(functional_gene_effects)), 1:min(6, ncol(functional_gene_effects))]

## 2.10 Effect-size plots for functional genes ---------------------------------------------
# Gene names in plotting order, per element cycle
C_names <- c(
  "Alpha-amylase", "Glucoamylase", "Pullulanase", "Isopullulanase",
  "Arabinofuranosidase", "Beta_mannanase", "Xylanase", "Xylose isomerase",
  "Beta-glucosidase", "Endoglucanase", "Exoglucanase",
  "Acetylglucosaminidase", "Endochitinase", "Exochitinase", "Pectinase",
  "Aryl-aldehyde oxidase", "Isocitrate lyase", "Limonene-1, 2-epoxide hydrolase",
  "Malate synthase", "Vanillate demethylase",
  "Glyoxal oxidase", "Phenol oxidase (tyrosinase)",
  "Methyl coenzyme M reductase", "Particulate methane monooxygenase",
  "Soluable methane monooxygenase",
  "Pyruvate oxidation", "Pyruvate <=> acetyl-CoA + formate", "Acetogenesis",
  "Lactate utilization", "Acetate to acetyl-CoA", "RuBisCo"
)
N_names <- c(
  "amoA", "amoB", "amoC", "hao", "nxrA", "nxrB", "narG", "narH", "narI", "napA",
  "napB", "nirK", "nirS", "norB", "norC", "nosZ", "HzsC", "HzsB", "HzsA", "hdh",
  "nrfA", "nrfH", "nirB", "nirD",
  "nasA", "nasB", "narB", "NR", "NIT-6", "nirA", "nifD", "nifK", "nifH", "nrtA", "nrtB", "nrtC",
  "nrtD", "nmo", "gdh_K00261", "gdh_K00262", "gdh_K15371", "glsA", "ureA", "ureC", "glnA"
)
S_names <- c(
  "sat", "cysC", "cysD", "cysNC", "csyH", "cysJ", "cysN", "aprA", "aprB", "dsrA", "dsrB",
  "asrA", "asrB", "asrC", "sir", "sor", "sreA", "sreB", "sreC", "hydA", "hydD", "hydB",
  "hydG", "psrA", "psrB", "psrC", "sqr", "fccA", "fccB", "soxD", "soxX", "soxA", "soxB",
  "soxC", "soxY", "soxZ", "ttrA", "ttrB", "ttrC", "phsA", "phsB", "phsC"
)
Other_names <- c(
  "Cyc1", "Cyc2", "MtrA", "MtrB", "MtrC", "arsC (grx)", "arsC (trx)",
  "aioA", "arsM", "ygfM", "xdhD", "YgfK"
)

# Background bands (axis positions) copied from the original script. The bands of
# the nitrogen, sulfur and other-element plots are identical to each other and
# do not span all genes; see the notes delivered with this script.
carbon_bands <- make_bands(
  edges = c(0.5, 4.5, 8.5, 11.5, 14.5, 15.5, 20.5, 22.5, 25.5, 30.5, 31.5),
  fills = c(
    "#E31A1C", "#893f45", "#e56eee", "#5fb236", "#1F78B4",
    "#6A3D9A", "#f19837", "#00bfff", "#ffdb00", "#8b008b"
  )
)
nitrogen_bands <- make_bands(
  edges = c(0.5, 4.5, 8.5, 11.5, 14.5, 15.5, 20.5, 22.5, 25.5, 27.5, 28.5, 29.5, 30.5),
  fills = c(
    "#ffdb00", "#f19837", "#e56eee", "#5fb236", "#9457eb", "#ffdb00",
    "#f19837", "#e56eee", "#5fb236", "#9457eb", "#ffdb00", "#e56eee"
  )
)
other_bands <- make_bands(
  edges = c(0.5, 4.5, 8.5, 11.5, 14.5),
  fills = c("#ffdb00", "#f19837", "#e56eee", "#5fb236")
)

eff.size_carbon_plot <- plot_effect_sizes(
  functional_gene_effects, C_names, carbon_bands,
  style = "bar", errorbar_width = 0.5
)
eff.size_carbon_plot

eff.size_nitrogen_plot <- plot_effect_sizes(
  functional_gene_effects, N_names, nitrogen_bands, style = "point"
)
eff.size_nitrogen_plot

eff.size_sulfur_plot <- plot_effect_sizes(
  functional_gene_effects, S_names, nitrogen_bands, style = "point"
)
eff.size_sulfur_plot

eff.size_other_plot <- plot_effect_sizes(
  functional_gene_effects, Other_names, other_bands, style = "point"
)
eff.size_other_plot

save_plot(eff.size_carbon_plot, "eff.size_carbon_plot.pdf", 5.5, 2.5)


# 3. ENVIRONMENTAL FACTORS -----------------------------------------------------------------
# metadata and the other input objects were loaded in section 2.

## 3.1 Descriptive statistics ------------------------------------------------------------------
# Full variable set of the original analysis, currently reduced to the list below:
# env_vars <- c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt",
#               "WHC", "SOC", "NH4_N", "NO3_N", "AP")
env_vars <- c(
  "Plant_richness", "AGB", "pH", "Soil_moisture", "Clay_Silt", "SOC",
  "NH4_N", "NO3_N", "AP"
)

env_stats <- metadata %>%
  group_by(Group) %>%
  get_summary_stats(env_vars, type = "common") %>% # or type = "mean_sd"
  mutate(Group = factor(Group, levels = group_levels)) %>%
  arrange(variable, Group)
env_stats

## 3.2 GLMMs: effect of collapse on environmental variables --------------------------------------
# GLMM     : variable ~ Group + (1 | Gully_id) (covariates removed).
# GLMM     : glmmTMB with Type-II Wald chi-square tests. Families by variable type:
#   Plant_richness  Poisson(log)       counts; check DHARMa dispersion
#   pH              Gaussian(identity) pH is already a log10 scale; LMM-equivalent
#   Clay_Silt       Beta(logit)        percentage, analysed as a proportion (0-1)
#   all others      Gamma(log), or Tweedie(log) if zeros are present
#                   (AGB, Soil_moisture [can exceed 100% on a dry-mass basis],
#                   SOC, NH4_N, NO3_N, AP; also BGB / WHC if re-enabled above)
env_model_data <- metadata %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Clay_Silt = to_proportion(Clay_Silt)
  )

# env_special_families <- list(
#   Plant_richness = poisson(link = "log"),
#   pH = gaussian(link = "identity"),
#   Clay_Silt = glmmTMB::beta_family(link = "logit")
# )
# env_rule_vars <- setdiff(env_vars, names(env_special_families))
# env_families <- c(
#   env_special_families,
#   map(set_names(env_rule_vars), \(v) family_nonneg_continuous(env_model_data[[v]]))
# )

env_glmm <- run_glmm_workflow(
  label = "environment",
  responses = env_vars,
  data = env_model_data,
  fixed_terms = glmm_fixed_terms,
  family = Gamma(link = "log")
)

env_result_glmm <- env_glmm$group_effects %>%
  transmute(
    Variables = id,
    fixed_factors = "Group",
    estimate, std.error,
    p.value,
    p.adj,
    sig = p_to_stars(p.value, ns_label = "ns")
  )
env_result_glmm
# save_table(env_result_glmm, "glmm_wald_tests_environment.csv")

## 3.3 Bar plot of environmental variables -----------------------------------------------------
facet_label_text <- list(
  "Plant_richness" = "Plant richness",
  "AGB" = expression(paste("AGB ", "(g/", "m"^2, ")")),
  "BGB" = expression(paste("BGB ", "(g/", "m"^2, ")")),
  "pH" = "pH",
  "Soil_moisture" = "Moisture (%)",
  "Clay_Silt" = "Clay+Silt (%)",
  "WHC" = "WHC (%)",
  "SOC" = "SOC (g/kg)",
  "NH4_N" = expression(paste("NH"[4]^"+", "-N", " (mg/kg)")),
  "NO3_N" = expression(paste("NO"[3]^"-", "-N", " (mg/kg)")),
  "AP" = "AP (mg/kg)"
)
facet_labeller <- function(labels) {
  list(lapply(as.character(labels[[1]]), \(v) as.expression(facet_label_text[[v]])))
}

env_plot <- env_stats %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    variable = factor(variable, levels = env_vars)
  ) %>%
  ggplot(aes(Group, mean, fill = Group)) +
  geom_bar(position = position_dodge(0.8), stat = "identity", colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = group_colours) +
  facet_wrap(~variable, scales = "free", ncol = 3, strip.position = "left", labeller = facet_labeller) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_rect(colour = NA),
    strip.text = element_text(colour = "black", size = 6, margin = margin()),
    axis.title = element_text(colour = "black", size = 6),
    axis.text = element_text(colour = "black", size = 6),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(colour = "black", size = 6),
    legend.text = element_text(colour = "black", size = 6),
    legend.key.size = unit(0.5, "cm")
  )

# Significance symbol (Group effect, from the GLMM) and panel letter per facet
env_stars <- tibble(variable = factor(env_vars, levels = env_vars)) %>%
  left_join(
    env_result_glmm %>%
      filter(fixed_factors == "Group") %>%
      transmute(variable = factor(Variables, levels = env_vars), sig.labels = sig),
    by = "variable"
  ) %>%
  mutate(
    sig.labels = replace_na(sig.labels, ""),
    x1 = 0.5, y1 = 0.95, x2 = 0.1, y2 = 1,
    panel.labels = letters[as.numeric(variable)]
  )

# Optional significance brackets between the two groups (y positions were set by
# hand for the original data):
# env_brackets <- tibble(
#   variable = factor(env_vars, levels = env_vars), x = 1, xend = 2,
#   y = c(17, 270, 8, 220, 72, 185, 19, 40, 2.4), yend = y
# )
# env_plot <- env_plot +
#   geom_segment(data = env_brackets, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE)

env_plot <- env_plot +
  ggpp::geom_text_npc(
    data = env_stars, aes(npcx = x1, npcy = y1, label = sig.labels),
    size = 2.5, inherit.aes = FALSE
  ) +
  ggpp::geom_text_npc(
    data = env_stars, aes(npcx = x2, npcy = y2, label = panel.labels),
    size = 2.5, inherit.aes = FALSE
  )
env_plot

# save_plot(env_plot, "env_comparison_plot1.pdf", 145, 90, "mm")

## 3.4 Environmental heterogeneity within groups ------------------------------------------------
# Euclidean distance between standardised environmental profiles of the plots
# of one gully and group (5 plots -> 10 pairwise distances per gully x group).
env_df <- metadata %>% select(all_of(env_vars))
env.dist <- vegdist(scale(env_df), method = "euclidean", na.rm = TRUE)

plot_ids <- c(
  "G1_C", "G1_T", "G2_C", "G2_T", "G3_C", "G3_T",
  "G4_C", "G4_T", "G5_C", "G5_T", "G6_C", "G6_T"
)
gully_names <- c(G1 = "EB", G2 = "ML", G3 = "RS", G4 = "SLH", G5 = "HSX", G6 = "HH")

env_dis_plot_table <- map_dfr(plot_ids, \(id) {
  samples <- grep(id, metadata$Sample_name, value = TRUE)
  tibble(
    Gully_id = unname(gully_names[sub("_.*", "", id)]),
    Group = if (grepl("_C$", id)) "Un-collapsed" else "Collapsed",
    distance = as.vector(usedist::dist_subset(env.dist, samples))
  )
})

env_dis_plot_df <- env_dis_plot_table %>%
  mutate(
    Gully_id = factor(Gully_id, levels = unname(gully_names)),
    Group = factor(Group, levels = group_levels)
  )

## 3.5 GLMM: effect of collapse on environmental heterogeneity -----------------------------------
# GLMM     : distance ~ Group + (1 | Gully_id) (covariates removed).
# Original : lmer(distance ~ Group + Time + Slope + MAP + (1 | Gully_id)).
# GLMM     : glmmTMB, Gamma(log): Euclidean distances are positive and continuous.
#            Random effect (1 | Gully_id); Type-II Wald chi-square tests.
env_dist_glmm <- run_glmm_workflow(
  label = "environmental_heterogeneity",
  responses = "distance",
  data = env_dis_plot_df,
  fixed_terms = glmm_fixed_terms,
  family = Gamma(link = "log")
)
env_dist_glmm_results <- tibble(
  fixed_factors = glmm_fixed_terms,
  Chisq = env_dist_glmm$effects[paste0(glmm_fixed_terms, ".chisq"), "distance"],
  p.value = env_dist_glmm$effects[paste0(glmm_fixed_terms, ".P"), "distance"]
) %>%
  mutate(sig = p_to_stars(p.value, ns_label = "ns"))
env_dist_glmm_results
# save_table(env_dist_glmm_results, "glmm_wald_tests_environmental_heterogeneity.csv")

# The significance symbol in the figure is taken from the GLMM (the original
# script had "**" typed in by hand).
env_dist_label <- env_dist_glmm_results %>%
  filter(fixed_factors == "Group") %>%
  pull(sig)

env_dist_plot <- env_dis_plot_df %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, colour = NA) +
  geom_boxplot(width = 0.4, size = 0.75, outlier.colour = NA) +
  geom_jitter(
    aes(fill = Group),
    shape = 21, size = 1.5,
    position = position_jitter(width = 0.2, seed = analysis_seed)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::geom_text_npc(
    data = tibble(npcx = 0.5, npcy = 0.95, label = env_dist_label),
    aes(npcx = npcx, npcy = npcy, label = label), inherit.aes = FALSE
  ) +
  labs(x = NULL, y = "Environmental heterogeneity") +
  scale_fill_manual(values = group_colours) +
  main_theme +
  theme(legend.position = "none")
env_dist_plot

# save_plot(env_dist_plot, "environmental_heterogeneity_plot.pdf", 45, 45, "mm")


# 4. SEQUENCING STATISTICS, MAG RECRUITMENT AND MAG QUALITY ---------------------------------------

## 4.1 Read MAG tables ----------------------------------------------------------------------------
# First row of the CoverM table is removed, as in the original script.
mag_abundance_tab <- read.delim(
  mag_abundance_file,
  header = TRUE, sep = "\t", row.names = 1, comment.char = "", check.names = FALSE
)[-1, ]
mag_tax_bin <- read.delim( # not used below
  mag_tax_file,
  header = TRUE, sep = "\t", comment.char = "", check.names = FALSE
)
mag_att_raw <- read.delim(
  mag_attr_file,
  header = TRUE, sep = "\t", row.names = 1, comment.char = "", check.names = FALSE
)
mag_abundance_tab[1:5, 1:5]

# Number of samples in which each MAG is detected, joined to the MAG attributes
mags_att <- mag_abundance_tab %>%
  as.data.frame() %>%
  mutate(across(where(is.numeric), ~ as.numeric(. > 0))) %>%
  mutate(rowsum = rowSums(.)) %>%
  select(rowsum) %>%
  rownames_to_column("Genome_id") %>%
  left_join(mag_att_raw, by = "Genome_id")

seq_depth <- read.csv(seq_depth_file, check.names = FALSE) %>%
  dplyr::rename(Sample_name = Sample) %>%
  dplyr::left_join(metadata, by = "Sample_name") %>%
  mutate(Group = factor(Group, levels = group_levels))

## 4.2 GLMM: MAG recruitment rate between groups ----------------------------------------------------
# Response : percentage of reads recruited by the MAGs (analysed as a proportion).
# Original : lmer(MAG_mapping_rate ~ Group + (1 | Gully_id)).
# GLMM     : glmmTMB, Beta(logit), because the rate is bounded in (0, 1);
#            random effect (1 | Gully_id).
mapping_rate_data <- seq_depth %>%
  mutate(MAG_mapping_rate = to_proportion(MAG_mapping_rate))

mapping_rate_glmm <- run_glmm_workflow(
  label = "mag_mapping_rate",
  responses = "MAG_mapping_rate",
  data = mapping_rate_data,
  fixed_terms = "Group",
  family = glmmTMB::beta_family(link = "logit")
)

# P-value of the Group effect from the GLMM
mapping_p <- mapping_rate_glmm$effects["Group.P", "MAG_mapping_rate"]
mapping_sig <- p_to_stars(mapping_p, ns_label = "ns")

# y position of the label and bracket, above the highest observation
y_top <- max(seq_depth$MAG_mapping_rate, na.rm = TRUE)
y_label <- y_top * 1.08

mag_mapping_rate <- ggplot(seq_depth, aes(x = Group,
                                          y = MAG_mapping_rate,
                                          fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    size = 2, alpha = 0.7,
    position = position_jitter(width = 0.08, height = 0, seed = analysis_seed)
  ) +
  # significance symbol from the GLMM
  annotate("text", x = 1.5, y = y_label, label = mapping_sig, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Group", y = "MAG recruitment rate (%)") +
  scale_fill_manual(values = group_colours) +
  theme_classic()
mag_mapping_rate

# save_plot(mag_mapping_rate, "mag_mapping_rate.pdf", 5, 3.5)

## 4.3 Sequencing depth, mapping rate and MAG quality ------------------------------------------------
# NOTE: this block is unchanged from the original. It assumes that `seq_depth`
# still contains a column named `Sample` and that the metrics have the names
# listed in the filter below; see the notes delivered with this script.
sequencing_metrics <- c("Depth", "contig_mapping", "Completeness",
                        "Contamination", "MAG_mapping_rate", 
                        "16S_clean_tags_per_sample")

basic_sequencing_stat <- rbind(
  seq_depth %>%
    select(c("Sample_name", "16S_clean_tags_per_sample", "Depth",
             "contig_mapping", "MAG_mapping_rate")) |>
    pivot_longer(col = -Sample_name, names_to = "variable", values_to = "value") %>%
    select(variable, value),
  data.frame(variable = rep("Completeness", nrow(mags_att)), value = mags_att$Completeness),
  data.frame(variable = rep("Contamination", nrow(mags_att)), value = mags_att$Contamination)
) %>%
  filter(variable %in% sequencing_metrics) %>%
  mutate(variable = factor(variable, levels = sequencing_metrics))

basic_sequencing_stat %>%
  group_by(variable) %>%
  summarise_variable(value)

basic_sequencing_stat_plot <- basic_sequencing_stat %>%
  ggplot(aes(x = 1, y = value)) +
  geom_half_violin(
    fill = "#79ceb8", position = position_nudge(x = 0.01),
    side = "r", width = 0.05, colour = NA
  ) +
  geom_boxplot(
    fill = "#e95f5c", width = 0.01, size = 1.2, lwd = 0.5,
    outlier.colour = NA, position = position_nudge(x = 0.01)
  ) +
  geom_jitter(
    colour = "#e95f5c", shape = 16, size = 1, stroke = 0.5, alpha = 0.5,
    position = position_jitter(width = 0.01, seed = analysis_seed)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  facet_wrap(. ~ variable, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Value") +
  main_theme
basic_sequencing_stat_plot

# save_plot(basic_sequencing_stat_plot, "basic_sequencing_stat_plot.pdf", 120, 80, "mm")


# 5. PROCRUSTES ANALYSIS OF 16S vs. MAG COMMUNITY STRUCTURE --------------------------------------------
# otu_16s (OTU x sample) and mag_comp (sample x MAG) come from the sourced scripts.

## 5.1 Match samples and compute Bray-Curtis dissimilarities -----------------------------------------------
otu_16S <- data.frame(t(otu_16s))
mag_comp <- data.frame(t(mags_abun_tab))
common_samples <- intersect(rownames(otu_16S), rownames(mag_comp))
length(common_samples)

# Keep the same samples in the same order in both tables
otu_16S <- otu_16S[common_samples, , drop = FALSE]
mag_comp <- mag_comp[common_samples, , drop = FALSE]
stopifnot(identical(rownames(otu_16S), rownames(mag_comp)))

# Dissimilarities are based on relative abundances (row totals = 1)
otu_16S_rel <- decostand(otu_16S, method = "total")
mag_rel <- decostand(mag_comp, method = "total")
dist_16s <- vegdist(otu_16S_rel, method = "bray")
dist_MAG <- vegdist(mag_rel, method = "bray")

## 5.2 NMDS ordination ---------------------------------------------------------------------------------------
set.seed(analysis_seed)
nmds_16s <- metaMDS(dist_16s, k = 2, trymax = 999, autotransform = FALSE)
nmds_MAG <- metaMDS(dist_MAG, k = 2, trymax = 999, autotransform = FALSE)
nmds_16s$stress
nmds_MAG$stress

## 5.3 Procrustes rotation and permutation test ---------------------------------------------------------------
procrustes_result <- procrustes(nmds_16s, nmds_MAG, symmetric = TRUE)
summary(procrustes_result)

set.seed(analysis_seed)
protest_result <- protest(nmds_16s, nmds_MAG, permutations = n_permutations)
protest_result

## 5.4 Procrustes plots -------------------------------------------------------------------------------------------
proc_df <- data.frame(
  X1 = procrustes_result$X[, 1], Y1 = procrustes_result$X[, 2],
  X2 = procrustes_result$Yrot[, 1], Y2 = procrustes_result$Yrot[, 2]
)

procrustes_ggplot <- ggplot(proc_df) +
  geom_point(aes(X1, Y1), colour = "blue") +
  geom_point(aes(X2, Y2), colour = "red") +
  geom_segment(
    aes(x = X1, y = Y1, xend = X2, yend = Y2),
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  labs(
    title = "Procrustes comparison of 16S and MAG community structure",
    x = "Dimension 1", y = "Dimension 2"
  ) +
  main_theme
procrustes_ggplot

# save_plot(procrustes_ggplot, "procrustes_plot.pdf", 5, 4)
