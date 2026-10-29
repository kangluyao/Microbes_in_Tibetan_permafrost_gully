# =============================================================================
# MAG functional traits: C-S-R strategy traits of MAGs that increase or
# decrease in collapsed gullies
# =============================================================================
# Purpose : (i) Classify MAGs by the direction of their mean abundance change
#           (Collapsed - Un-collapsed), (ii) compare microtrait-based life-
#           history, resource and stress-tolerance traits between the two
#           MAG groups, (iii) examine correlations among traits across MAGs,
#           and (iv) test community-weighted mean (CWM) traits between
#           Collapsed and Un-collapsed samples with Gaussian GLMMs.
# Design  : Samples are nested in gullies, so Gully_id is a random intercept
#           in the GLMMs.
# Inputs  : - script/0.read_data.R
#           - microtrait results (RDS) and MAG growth-rate table (CSV) under
#             data/metagenome/MAGs
# Outputs : Figures (ggsave calls are commented out) and optional CSV exports
#           of the selected traits and CWM table (commented out).
#
# Sections
#   1. Setup (paths, parameters, packages, plot theme)
#   2. Data input
#   3. Helper functions
#   4. MAG abundance change and effect groups
#   5. MAG traits (comparison between effect groups, trait correlations)
#   6. Community-weighted mean (CWM) traits
#   7. GLMMs of CWM traits
#   8. Combined figure
# =============================================================================

rm(list = ls())  # start from a clean workspace


# 1. Setup ====================================================================

## 1.1 Paths and parameters ---------------------------------------------------
project_dir <- "E:/thermokarst_gully"
setwd(project_dir)  # needed because the source() call below uses a relative path

# Defined before sourcing, as in the original (the sourced script may use them)
wd_16s   <- file.path(getwd(), "data/16S/rdp")
wd_fun   <- file.path(getwd(), "data/metagenome")
save.dir <- file.path(getwd(), "result")

mag_dir    <- file.path(wd_fun, "MAGs")
result_dir <- file.path(project_dir, "revision/result")  # CSV and figures

group_levels <- c("Un-collapsed", "Collapsed")
alpha_level  <- 0.05              # significance threshold used in the effect plot
effect_ylim  <- c(-2.5, 2.5)      # y-range of the effect-size plot

# Colours of the MAG groups (direction of abundance change)
effect_colors <- c("Negative" = "#79ceb8", "Positive" = "#e95f5c")

## 1.2 Packages ---------------------------------------------------------------
library(vegan)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(edgeR)
library(ggplot2)
library(microtrait)   # trait.normalize(), add.metadata()
library(gghalves)     # geom_half_violin()
library(ggpubr)       # stat_compare_means()
library(car)          # logit() transform for bounded [0, 1] data
library(glmmTMB)      # generalized linear mixed models
library(DHARMa)       # simulation-based residual diagnostics for GLMMs
library(purrr)        # functional programming / iteration
library(broom.mixed)  # tidy() extraction of model coefficients
library(dplyr)
library(stringr)
library(tibble)
library(cowplot)

## 1.3 Plot theme -------------------------------------------------------------
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


# 2. Data input ===============================================================
source("script/0.read_data.R")


# 3. Helper functions =========================================================

## 3.1 GLMM fitting and diagnostics -------------------------------------------

#' Fit one glmmTMB model per response variable in `responses`.
#' Each call is wrapped in tryCatch so that one failing model does not stop
#' the loop. Response names are back-tick quoted to safely handle illegal
#' characters (e.g. "Fen-183", which R would otherwise parse as "Fen - 183").
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

#' Extract fixed-effect coefficients for a given term (default: "Group")
#' from a named list of fitted models, with Benjamini-Hochberg (FDR)
#' correction across all models in the list.
#' `component` selects "cond" (conditional model) or "zi" (zero-inflation part).
extract_group_pvals <- function(fits, term_pattern = "^Group",
                                component = "cond") {
  map_dfr(fits, tidy, effects = "fixed", .id = "id") |>
    filter(.data$component == !!component, str_detect(term, term_pattern)) |>
    mutate(p.adj = p.adjust(p.value, method = "BH")) |>
    select(id, term, estimate, std.error, p.value, p.adj)
}

#' Run DHARMa simulated-residual diagnostics on a batch of fitted models.
#' Returns the raw DHARMa objects (for plotting) and a summary table of the
#' uniformity (KS test), dispersion and zero-inflation p-values.
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

#' Save all DHARMa diagnostic plots from a batch to a single PDF.
save_dharma_plots <- function(dharma_res, file) {
  pdf(file, width = 8, height = 6)
  iwalk(dharma_res, \(res, v) plot(res, main = v))
  dev.off()
}

## 3.2 Correlation and plotting helpers ---------------------------------------

#' P-values of all pairwise correlations between the columns of `x` and `y`
#' (rows of the result = columns of `x`). Failed tests return NA.
cross_cor_p <- function(x, y, method = "spearman") {
  x <- as.matrix(x)
  y <- as.matrix(y)
  p <- matrix(NA_real_, ncol(x), ncol(y),
              dimnames = list(colnames(x), colnames(y)))
  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(y))) {
      p[i, j] <- tryCatch(
        cor.test(x[, i], y[, j], method = method, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }
  }
  p
}

#' Significance stars for p-values: *** < 0.001, ** < 0.01, * < 0.05.
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE,
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", "")))
}


# 4. MAG abundance change and effect groups ===================================

## 4.1 MAG abundance table with sample information ----------------------------
# Gully and group labels are derived from the sample names ("G1"-"G6";
# "_C" = Un-collapsed, "_T" = Collapsed).
mags_abun <- mags_abun_tab |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = "Sample_name") |>
  mutate(Gully_id = case_when(grepl("G1", Sample_name) ~ "EB",
                              grepl("G2", Sample_name) ~ "ML",
                              grepl("G3", Sample_name) ~ "RS",
                              grepl("G4", Sample_name) ~ "SLH",
                              grepl("G5", Sample_name) ~ "HSX",
                              grepl("G6", Sample_name) ~ "HH")) |>
  mutate(Group = case_when(grepl("_C", Sample_name) ~ "Un-collapsed",
                           grepl("_T", Sample_name) ~ "Collapsed")) |>
  left_join(metadata |>
              select(Sample_name, Time, Slope, MAP),
            by = "Sample_name") |>
  relocate(c(Group, Time, Slope, MAP, Gully_id), .before = 1) |>
  select(-Sample_name) |>
  data.frame()

## 4.2 MAGs that increase or decrease in collapsed gullies --------------------
# Classification is based only on the sign of the difference between group
# means (Collapsed - Un-collapsed); no statistical test or threshold is used.
mags_result_diff <- mags_abun |>
  select(-Time, -Slope, -MAP, -Gully_id) |>
  group_by(Group) |>
  summarise(across(everything(), mean), .groups = "drop") |>
  column_to_rownames("Group") |>
  t() |>
  as.data.frame() |>
  mutate(diff = Collapsed - `Un-collapsed`,
         eff.group = case_when(diff > 0 ~ "Positive",
                               diff < 0 ~ "Negative",
                               TRUE     ~ "No change")) |>
  rownames_to_column("Genome_id")

nrow(mags_result_diff |> filter(eff.group == "Positive"))
nrow(mags_result_diff |> filter(eff.group == "Negative"))

## 4.3 Genus composition of the positive and negative MAG groups --------------
mags_tax_group_df <- mags_result_diff |>
  left_join(mags_tax, by = "Genome_id") |>
  select(Genus, eff.group) |>
  filter(eff.group %in% c("Positive", "Negative")) |>
  mutate(Number = rep(1, n())) |>
  summarize(across(everything(), sum), .by = c(Genus, eff.group)) |>
  mutate(Genus = gsub("g__", "", Genus)) |>
  mutate(Genus = factor(Genus),
         Genus = fct_reorder(Genus, Number, .desc = TRUE))


# 5. MAG traits ===============================================================

## 5.1 Microtrait results -----------------------------------------------------
microtrait_results <- readRDS(
  file.path(mag_dir, "microtraits/thermokarst_gully.microtraitresults.rds"))

# Normalise the trait matrices by genome length
microtrait_results_metadata_norm <- microtrait_results |>
  trait.normalize(normby = "genome_length")

## 5.2 Additional MAG attributes ----------------------------------------------
# GC content, maximum growth rate (1 / doubling time `d`) and effect group
growth_rate <- read.table(file.path(mag_dir, "growth_rate_results.csv"),
                          header = TRUE, sep = ",", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE,
                          comment.char = "", check.names = FALSE)

mags_all_attrs <- mags.att |>
  left_join(growth_rate, by = "Genome_id") |>
  left_join(mags_result_diff, by = "Genome_id") |>
  mutate(`Max growth rate` = 1 / d) |>
  select(Genome_id, GC, `Max growth rate`, eff.group)

# Attach the MAG attributes to the normalised microtrait results
trait_data <- microtrait::add.metadata(microtrait_results_metadata_norm,
                                       mags_all_attrs,
                                       genome_metadata_idcol = "Genome_id")

## 5.3 Trait selection --------------------------------------------------------
# Traits of interest, in the order used in the figures: life history
# (1-2), resource acquisition (3-10), resource use (11-15) and stress
# tolerance (16-20)
strategy_vars <- c("genome_length", "Max growth rate",
                   "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization",
                   "Resource Acquisition:Substrate degradation:simple compound degradation",
                   "Resource Acquisition:Substrate uptake:carbohydrate transport",
                   "Resource Acquisition:Substrate uptake:N compound transport",
                   "Resource Acquisition:Substrate uptake:organophosphorus transport",
                   "Resource Acquisition:Substrate assimilation:C1 compounds",
                   "Resource Acquisition:Substrate assimilation:N compounds",
                   "Resource Acquisition:Substrate assimilation:P compounds",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration",
                   "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation",
                   "Resource Use:Chemotrophy:chemolithoautotrophy",
                   "Resource Use:Phototrophy:photosystem",
                   "Stress Tolerance:Specific:temperature:low",
                   "Stress Tolerance:Specific:pH stress",
                   "Stress Tolerance:Specific:desiccation/osmotic/salt stress",
                   "Stress Tolerance:Specific:oxidative stress",
                   "Stress Tolerance:Specific:envelope stress")

# Genome length is expressed in Mb and rounded to one decimal
trait_sel_df <- trait_data$trait_matrixatgranularity1 |>
  column_to_rownames("id") |>
  mutate(genome_length = round(genome_length / 1000000, 1)) |>
  select(c("eff.group", all_of(strategy_vars)))

# write.csv(trait_sel_df, file = file.path(result_dir, "trait_sel_df.csv"),
#           row.names = TRUE)

## 5.4 Trait comparison between negative and positive MAGs --------------------
# Wilcoxon rank-sum test per trait (one comparison per panel)
my_comparisons <- list(c("Negative", "Positive"))

microtraits_comparison <- trait_sel_df |>
  filter(eff.group %in% c("Negative", "Positive")) |>
  gather(mag.attrs, value, -c("eff.group")) |>
  mutate(eff.group = factor(eff.group, levels = c("Negative", "Positive"))) |>
  mutate(mag.attrs = factor(mag.attrs, levels = strategy_vars)) |>
  ggplot(aes(eff.group, value, fill = eff.group)) +
  geom_half_violin(position = position_nudge(x = 0.25),
                   side = "r", width = 0.6, color = NA, alpha = 0.65) +
  # geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = eff.group, colour = eff.group),
              shape = 16, size = 0.5, height = 0,
              width = 0.15, alpha = 0.65) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.35, linewidth = 0.25) +
  stat_compare_means(comparisons = my_comparisons, paired = FALSE,
                     p.adjust.method = "BH", label = "p.signif",
                     bracket.size = 0.25, bracket.width = 0.1, size = 2.5,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = effect_colors) +
  scale_color_manual(values = effect_colors) +
  facet_wrap(~mag.attrs, scales = "free_y", ncol = 5) +
  main_theme +
  theme(legend.position = "none",
        strip.text = element_text(colour = "black", size = 5),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
microtraits_comparison

# ggsave(file.path(result_dir, "microtraits_comparison_2026.10.6.pdf"),
#        microtraits_comparison, width = 170, height = 120, units = "mm")

## 5.5 Correlations among traits across MAGs ----------------------------------
# Spearman correlations between two sets of trait columns, with BH-adjusted
# p-values over all pairs.
# WARNING: `trait_sel_df` has `eff.group` in column 1, so columns 3:15 and
# 16:20 are not the trait blocks implied by `strategy_vars` positions; see the
# "potential problems" notes before using these results.
trait_sel_df <- trait_sel_df |> 
  select(-eff.group)  # drop the effect group column for correlation analysi
cross_cor_indi <- cor(trait_sel_df[, 3:15],
                      trait_sel_df[, 16:20],
                      use = "pairwise.complete.obs",
                      method = "spearman")

p_indi <- cross_cor_p(trait_sel_df[, 3:15],
                      trait_sel_df[, 16:20],
                      method = "spearman")

p_indi_adj <- matrix(p.adjust(p_indi, method = "BH"),
                     nrow(p_indi), ncol(p_indi),
                     dimnames = dimnames(p_indi))

# Coefficients displayed in the cells
corr_indiv_plot <- corrplot::corrplot(cross_cor_indi,
                                      method = "color",
                                      tl.col = "black",
                                      tl.srt = 45,
                                      tl.cex = 0.5,
                                      cl.cex = 0.5,
                                      addCoef.col = "black",
                                      number.cex = 0.5)

# Significance levels displayed as stars (BH-adjusted p-values)
corrplot::corrplot(cross_cor_indi,
                   method    = "color",
                   p.mat     = p_indi_adj,
                   sig.level = c(0.001, 0.01, 0.05),
                   insig     = "label_sig",
                   pch.cex   = 1,
                   pch.col   = "black",
                   tl.col    = "black",
                   tl.srt    = 45,
                   tl.cex    = 0.5,
                   cl.cex    = 0.5)


# 6. Community-weighted mean (CWM) traits =====================================
# CWM of trait t in community j = sum_k(p_kj * t_k) / sum_k(p_kj), where p_kj is
# the relative abundance of MAG k in sample j. The trait table is first
# re-ordered to match the rows of the abundance table.
trait_sel_df <- trait_sel_df[rownames(mags_abun_tab), ]
relative_abundance <- mags_abun_tab / 100   # percent -> proportion

cwm_results <- data.frame(matrix(NA, nrow = ncol(relative_abundance),
                                 ncol = ncol(trait_sel_df)))
row.names(cwm_results) <- colnames(relative_abundance)
colnames(cwm_results)  <- colnames(trait_sel_df)

for (i in 1:ncol(trait_sel_df)) {
  for (j in 1:ncol(relative_abundance)) {
    cwm_results[j, i] <- sum(relative_abundance[, j] * trait_sel_df[, i],
                             na.rm = TRUE) / sum(relative_abundance[, j])
  }
}

print(cwm_results)
# write.csv(cwm_results, file = file.path(result_dir, "cwm_results.csv"),
#           row.names = TRUE)

# Keep columns with a positive sum and order the samples as in `metadata`
# (inspected only; not used in the models below)
cwm_results_final <- cwm_results[, colSums(cwm_results[]) > 0]
cwm_results_final <- cwm_results_final[metadata$Sample_name, ]

dim(cwm_results_final)


# 7. GLMMs of CWM traits ======================================================

## 7.1 Model data -------------------------------------------------------------
# CWM traits are z-scored here (scale()), so they can be negative; Gamma or
# Tweedie families are therefore not options and a Gaussian GLMM is used.
# Traits with missing values (e.g. constant traits) are dropped.
trait_df <- cwm_results |>
  cbind(metadata[, c("Gully_id", "Group")]) |>
  select(all_of(c("Group", "Gully_id", strategy_vars))) |>
  mutate(across(where(is.numeric), scale)) |>
  mutate(Group = factor(Group, levels = group_levels)) |>
  select(where(~ !any(is.na(.))))

id.traits <- intersect(strategy_vars, names(trait_df))

# Sanitise illegal column-name characters so that the names can be used in
# model formulas; id.traits is sanitised identically.
colnames(trait_df) <- make.names(colnames(trait_df))
id.traits <- make.names(id.traits)

## 7.2 Gaussian GLMMs ---------------------------------------------------------
fits.traits <- fit_glmm_batch(id.traits, trait_df,
                              family = gaussian(link = "identity"))

trait_pvals <- extract_group_pvals(fits.traits)
print(trait_pvals)

trait_diag <- run_dharma_diagnostics(fits.traits)
print(trait_diag$diagnostics)
# save_dharma_plots(trait_diag$residuals, "DHARMa_diagnostics_traits.pdf")

## 7.3 Effect-size plot -------------------------------------------------------
# Effect (Collapsed vs. Un-collapsed) +/- SE per trait; points are coloured by
# direction when the BH-adjusted p-value is below `alpha_level`.
# Background shading marks the trait categories (see 5.3).
each_trait_comparison <- trait_pvals |>
  as.data.frame() |>
  mutate(sig = as.vector(unlist(lapply(p.adj, p.stars)))) |>
  mutate(id = factor(id, levels = id)) |>
  mutate(colour = case_when(estimate <= 0 & p.adj <= alpha_level ~ "Negative",
                            estimate > 0 & p.adj <= alpha_level ~ "Positive",
                            p.adj > alpha_level ~ "Neutral")) |>
  ggplot(aes(x = id, y = estimate, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375, colour = "gray2") +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = id,
                y = (estimate / abs(estimate)) *
                  (abs(estimate) + std.error) * 1.25),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values = c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limits = effect_ylim) +
  scale_x_discrete(position = "bottom", labels = strategy_vars) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = effect_ylim[1],
           ymax = effect_ylim[2], alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 2.5, xmax = 10.5, ymin = effect_ylim[1],
           ymax = effect_ylim[2], alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 10.5, xmax = 15.5, ymin = effect_ylim[1],
           ymax = effect_ylim[2], alpha = 0.1, fill = "#5cc3e8") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = effect_ylim[1],
           ymax = effect_ylim[2], alpha = 0.1, fill = "#ffdb00") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
each_trait_comparison

# ggsave(file.path(result_dir, "each_trait_comparison_2026.10.6.pdf"),
#        each_trait_comparison, width = 173, height = 110, units = "mm")


# 8. Combined figure ==========================================================
# Top: trait comparison between MAG groups; bottom: GLMM effect sizes
mags_trait_plot <- plot_grid(microtraits_comparison, each_trait_comparison,
                             nrow = 2, rel_heights = c(3, 1))
mags_trait_plot

# ggsave(file.path(result_dir, "mags_trait_plot.pdf"),
#        mags_trait_plot, width = 183, height = 155, units = "mm")