# =============================================================================
# SingleM community analysis: collapsed vs. un-collapsed thermokarst gullies
# =============================================================================
# Purpose : Compare microbial community composition, differential abundance,
#           alpha diversity and beta diversity (SingleM profiles) between
#           Collapsed (sample suffix "_T") and Un-collapsed ("_C") samples.
# Design  : Samples are clustered within gullies (Gully_id), so Gully_id is
#           treated as a random intercept (GLMMs) or as a permutation block
#           (PERMANOVA).
# Inputs  : - script/0.read_data.R, script/0.functions.R (sourced)
#           - SingleM tables: gully-phylum.tsv, gully-genus.tsv
# Outputs : Figures (ggsave calls are commented out) and printed model tables.
#
# Sections
#   1. Setup (paths, parameters, packages)
#   2. Helper functions
#   3. Data input
#   4. Phylum-level community composition
#   5. Phylum-level differential abundance (beta-regression GLMMs)
#   6. Genus-level differential abundance (beta-regression GLMMs)
#   7. Shared genera (Venn diagram)
#   8. Alpha diversity (GLMMs)
#   9. Beta diversity (PERMANOVA and PCoA)
#  10. Combined figure
# =============================================================================

rm(list = ls())  # start from a clean workspace

# 1. Setup ====================================================================

## 1.1 Paths and parameters ---------------------------------------------------
project_dir <- "E:/thermokarst_gully"
setwd(project_dir)  # needed because the source() calls below use relative paths

save.dir    <- file.path(getwd(), "result")  # kept from the original script
singlem_dir <- file.path(project_dir, "data/metagenome/singlem")
fig_dir     <- file.path(project_dir, "revision/result")  # used by ggsave()

group_levels <- c("Un-collapsed", "Collapsed")  # reference level first
n_top_taxa   <- 10                              # taxa shown / modelled
n_perm       <- 999                             # PERMANOVA permutations

## 1.2 Packages ---------------------------------------------------------------
library(car)          # logit() transform for bounded [0, 1] data
library(glmmTMB)      # generalized linear mixed models
library(DHARMa)       # simulation-based residual diagnostics for GLMMs
library(purrr)        # functional programming / iteration
library(broom.mixed)  # tidy() extraction of model coefficients
library(dplyr)
library(stringr)
library(tibble)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(gghalves)
library(ggvenn)
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
group_colors <- c("Un-collapsed" = "#79ceb8", "Collapsed" = "#e95f5c")

# 2. Helper functions =========================================================

## 2.1 Model fitting and diagnostics ------------------------------------------

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

## 2.2 SingleM table handling -------------------------------------------------

#' Read a SingleM taxonomy table, drop unassigned reads and reduce the lineage
#' string to the name at one rank (`rank_prefix`: "p" = phylum, "g" = genus).
#' Returns a data frame with taxa as row names and samples as columns.
read_singlem_table <- function(file, rank_prefix) {
  fread(file) |>
    filter(taxonomy != "unassigned") |>
    mutate(taxonomy = str_extract(taxonomy, paste0(rank_prefix, "__[^;]+"))) |>
    mutate(taxonomy = str_replace_all(taxonomy, ".*__", "")) |>
    column_to_rownames("taxonomy")
}

#' Convert a taxa-by-sample table to long format and assign the land-type
#' group from the sample-name suffix ("_C" = Un-collapsed, "_T" = Collapsed).
to_long_abundance <- function(abundance) {
  abundance |>
    rownames_to_column("taxonomy") |>
    pivot_longer(-taxonomy, names_to = "Sample",
                 values_to = "Relative_abundance") |>
    mutate(Group = case_when(grepl("_C", Sample) ~ "Un-collapsed",
                             grepl("_T", Sample) ~ "Collapsed"),
           Group = factor(Group, levels = group_levels))
}

#' Mean, SD and n of relative abundance for each taxon within each group.
summarise_by_group <- function(long_df) {
  long_df |>
    group_by(Group, taxonomy) |>
    summarise(across(Relative_abundance,
                     list(mean = ~mean(., na.rm = TRUE),
                          sd   = ~sd(., na.rm = TRUE),
                          n    = ~n())),
              .groups = "drop")
}

#' Names of the `n` taxa with the highest summed group-mean abundance.
select_top_taxa <- function(group_mean_df, n = n_top_taxa) {
  group_mean_df |>
    group_by(taxonomy) |>
    summarise(total = sum(Relative_abundance_mean)) |>
    slice_max(total, n = n) |>
    pull(taxonomy)
}


# 3. Data input ===============================================================
source("script/0.read_data.R")

# 4. Phylum-level community composition =======================================
singlem_phyla <- read_singlem_table(file.path(singlem_dir, "gully-phylum.tsv"),
                                    rank_prefix = "p")
singlem_phyla_df <- to_long_abundance(singlem_phyla)
singlem_phyla_group_ave_df <- summarise_by_group(singlem_phyla_df)

# Top phyla ranked by summed group-mean abundance; the remainder is "Others"
top_phyla <- select_top_taxa(singlem_phyla_group_ave_df)

plot_df <- singlem_phyla_group_ave_df |>
  mutate(taxonomy = ifelse(taxonomy %in% top_phyla, taxonomy, "Others")) |>
  mutate(taxonomy = factor(taxonomy, levels = c(top_phyla, "Others")))

# Stacked bar plot of mean relative abundance (stacking sums the "Others" rows)
composition_plot <- ggplot(plot_df, aes(x = Group, y = Relative_abundance_mean,
                                        fill = taxonomy)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_d3("category20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "Relative abundance (%)", x = NULL) +
  main_theme +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_blank())
composition_plot


# 5. Phylum-level differential abundance ======================================
# Beta-regression GLMM (logit link) per phylum, with Gully_id as a random
# intercept. Beta regression requires responses strictly within (0, 1), so
# percentages are converted to proportions (no boundary adjustment is applied
# to phyla, unlike genera in section 7).
singlem_top10_phyla <- data.frame(
  metadata[, c("Gully_id", "Group")],
  t(singlem_phyla[top_phyla, ]),
  check.names = FALSE
) |>
  mutate(Group = factor(Group, levels = group_levels)) |>
  mutate(across(all_of(top_phyla), \(x) x / 100))

fits.phy <- fit_glmm_batch(top_phyla, singlem_top10_phyla,
                           family = beta_family(link = "logit"))

phy_pvals <- extract_group_pvals(fits.phy)
print(phy_pvals)

phy_diag <- run_dharma_diagnostics(fits.phy)
print(phy_diag$diagnostics)
# save_dharma_plots(phy_diag$residuals, "DHARMa_diagnostics_phyla.pdf")


# 6. Genus-level differential abundance =======================================

## 6.1 Genus table and top genera ---------------------------------------------
singlem_genu <- read_singlem_table(file.path(singlem_dir, "gully-genus.tsv"),
                                   rank_prefix = "g")
singlem_genu_df <- to_long_abundance(singlem_genu)
singlem_genu_group_ave_df <- summarise_by_group(singlem_genu_df)

# Top genera within each group (not used downstream; kept for inspection)
singlem_10genu_df_unclo <- singlem_genu_group_ave_df |>
  filter(Group == "Un-collapsed") |>
  slice_max(Relative_abundance_mean, n = n_top_taxa)

singlem_genu_df_clo <- singlem_genu_group_ave_df |>
  filter(Group == "Collapsed") |>
  slice_max(Relative_abundance_mean, n = n_top_taxa)

# Genera modelled below: top genera by summed group-mean abundance
top_genera <- select_top_taxa(singlem_genu_group_ave_df)

## 6.2 Model data -------------------------------------------------------------
singlem_top10_genu <- data.frame(
  metadata[, c("Gully_id", "Group")],
  t(singlem_genu[top_genera, ]),
  check.names = FALSE
) |>
  mutate(across(all_of(top_genera), \(x) x / 100)) |>
  mutate(Group = factor(Group, levels = group_levels))

id.gen <- intersect(top_genera, names(singlem_top10_genu))

# Sanitise illegal column-name characters (e.g. "Fen-183" -> "Fen.183") so the
# names are safe to use in model formulas; id.gen is sanitised identically.
colnames(singlem_top10_genu) <- make.names(colnames(singlem_top10_genu))
id.gen <- make.names(id.gen)

# According to the original analysis, only 2 zero values occur across all
# genera. Instead of a full zero-inflation model (overkill and prone to
# quasi-separation with so few zeros), proportions are shrunk away from the
# boundaries (Smithson & Verkuilen, 2006). The transform is applied to every
# observation of each genus, not only to the zeros.
shrink_beta <- function(x, n = length(x)) (x * (n - 1) + 0.5) / n

singlem_top10_genu <- singlem_top10_genu |>
  mutate(across(all_of(id.gen), shrink_beta))

## 6.3 Beta-regression GLMMs --------------------------------------------------
fits.gen <- fit_glmm_batch(id.gen, singlem_top10_genu,
                           family = beta_family(link = "logit"))

gen_pvals <- extract_group_pvals(fits.gen)
print(gen_pvals)

gen_diag <- run_dharma_diagnostics(fits.gen)
print(gen_diag$diagnostics)
# save_dharma_plots(gen_diag$residuals, "DHARMa_diagnostics_genera.pdf")

## 6.4 Plot: dominant genera --------------------------------------------------
# Genera shown in the figure (selected manually). Original spelling is used
# here to match `singlem_genu_df$taxonomy`, whereas the model IDs were
# sanitised with make.names(); the "Fen.183" ID is therefore restored below.
plot_genera <- c("Defluviicoccus", "Rhodoplanes",
                 "Deferrimicrobium", "Fen-183", "Methyloceanibacter")
gen_pvals$id[gen_pvals$id == "Fen.183"] <- "Fen-183"

# P labels are taken from the model output (BH-adjusted p-values)
sig_labs_gen <- tibble(
  taxonomy = factor(plot_genera, levels = plot_genera),  # facet variable name
  P        = gen_pvals[["p.adj"]][match(plot_genera, gen_pvals[["id"]])],
  label    = ifelse(P < 0.001, "P < 0.001", sprintf("P = %.3f", P)),
  npcx     = 0.5,
  npcy     = 0.95
)

dominant_genera_abun_diff_plot <- singlem_genu_df |>
  filter(taxonomy %in% plot_genera) |>
  mutate(Group    = factor(Group, levels = group_levels),
         taxonomy = factor(taxonomy, levels = plot_genera)) |>
  ggplot(aes(Group, Relative_abundance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25),
                   side = "r", width = 0.65, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.3, size = 0.45, outlier.color = NA, alpha = 0.65) +
  geom_jitter(aes(color = Group), shape = 16, size = 1,
              width = 0.2, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  ggpp::geom_text_npc(data = sig_labs_gen,
                      aes(npcx = npcx, npcy = npcy, label = label)) +
  labs(x = "Group", y = "Relative abundance (%)") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~taxonomy, scales = "free", ncol = 3) +
  main_theme +
  theme(legend.position = "none")
dominant_genera_abun_diff_plot
# ggsave(file.path(fig_dir, "dominant_genera_abun_diff_plot.pdf"),
#        dominant_genera_abun_diff_plot,
#        width = 120, height = 80, units = "mm", dpi = 600)


# 7. Shared genera (Venn diagram) =============================================
# Genera detected (summed abundance > 0) in at least one sample of each group.
genera_detected_in <- function(sample_pattern) {
  samples <- grep(sample_pattern, metadata$Sample_name, value = TRUE)
  singlem_genu |>
    data.frame() |>
    mutate(rowsum = rowSums(pick(all_of(samples)))) |>
    filter(rowsum > 0) |>
    rownames()
}

venn_singlem_df <- list(
  Uncollapsed = genera_detected_in("_C"),
  Collapsed   = genera_detected_in("_T")
)

venn_singlem <- ggvenn(
  venn_singlem_df,
  fill_color = unname(group_colors),
  stroke_color = NA,
  set_name_size = 3,
  text_size = 2.5,
  show_percentage = TRUE,
  auto_scale = FALSE
)
venn_singlem
# ggsave(file.path(fig_dir, "venn_singlem.pdf"),
#        venn_singlem, width = 89, height = 89, units = "mm", dpi = 600)


# 8. Alpha diversity ==========================================================

## 8.1 Diversity indices ------------------------------------------------------
# Matrix form of the genus table (taxa x samples) for vegan functions
singlem_genu_mat <- as.matrix(singlem_genu)

# Overview table (not used downstream; the same indices are recomputed in
# `div_table` below together with the sample metadata)
shannon  <- diversity(t(singlem_genu_mat), index = "shannon")
richness <- specnumber(t(singlem_genu_mat))
alpha_div <- data.frame(
  Richness = richness,
  Shannon  = shannon,
  Evenness = shannon / log(richness),
  stringsAsFactors = FALSE
)

div_index <- c("Richness", "Shannon", "Evenness")
comm <- t(singlem_genu_mat)  # samples x taxa
# Safeguard: sample order of the community matrix must match the metadata
stopifnot(identical(rownames(comm), metadata$Sample_name))

div_table <- metadata |>
  select(Sample_name, Gully_id, Group) |>
  mutate(Group    = factor(Group, levels = group_levels),
         Richness = specnumber(comm),
         Shannon  = diversity(comm, "shannon"),
         Evenness = Shannon / log(Richness))

# Descriptive statistics per gully and group (printed only)
div_table |>
  group_by(Gully_id, Group) |>
  summarise(across(all_of(div_index),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        n    = ~n())),
            .groups = "drop")

## 8.2 GLMMs ------------------------------------------------------------------
# Family choice per index: Gamma (log link) for the continuous, positive
# Shannon index; negative binomial for richness (counts); beta (logit link)
# for evenness (bounded in (0, 1)). All models include Gully_id as a random
# intercept.
fits.div <- list(
  shannon  = glmmTMB(Shannon ~ Group + (1 | Gully_id), data = div_table,
                     family = Gamma(link = "log")),
  richness = glmmTMB(Richness ~ Group + (1 | Gully_id), data = div_table,
                     family = nbinom2),
  evenness = glmmTMB(Evenness ~ Group + (1 | Gully_id), data = div_table,
                     family = beta_family(link = "logit"))
)

div_pvals <- extract_group_pvals(fits.div)
print(div_pvals)

div_diag <- run_dharma_diagnostics(fits.div)
print(div_diag$diagnostics)
# iwalk(div_diag$residuals, \(res, v) plot(res, main = v))

## 8.3 Plot -------------------------------------------------------------------
# P labels are taken from the model output instead of being typed by hand.
# WARNING: `div_pvals` follows the order of `fits.div` (shannon, richness,
# evenness), whereas `div_index` is ordered Richness, Shannon, Evenness. The
# line below pairs them by position, as in the original script; see the
# "potential problems" notes before using these labels.
sig_labs_div <- tibble(
  Index = factor(div_index, levels = div_index),
  P     = div_pvals[["p.adj"]],
  label = ifelse(P < 0.001, "P < 0.001", sprintf("P = %.3f", P)),
  npcx  = 0.5,
  npcy  = 0.95
)

div_plot <- div_table |>
  pivot_longer(all_of(div_index), names_to = "Index", values_to = "value") |>
  mutate(Index = factor(Index, levels = div_index)) |>
  ggplot(aes(Group, value)) +
  geom_half_violin(aes(fill = Group), color = NA, side = "r", width = 0.5,
                   alpha = 0.65, position = position_nudge(x = 0.25)) +
  geom_boxplot(aes(fill = Group), width = 0.35, size = 0.3,
               outlier.color = NA, alpha = 0.65) +
  geom_jitter(aes(colour = Group), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  ggpp::geom_text_npc(data = sig_labs_div,
                      aes(npcx = npcx, npcy = npcy, label = label),
                      size = 2, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~Index, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = NULL) +
  main_theme +
  theme(legend.position = "none")
div_plot


# 9. Beta diversity ==========================================================

## 9.1 Bray-Curtis dissimilarity and PERMANOVA -------------------------------
tax_singlem_dist <- vegdist(t(singlem_genu_mat), "bray")

# Permutations are restricted within gullies (blocks) to respect the clustered
# design.
set.seed(123)
perm_scheme <- how(blocks = metadata$Gully_id, nperm = n_perm)
perm_test <- adonis2(tax_singlem_dist ~ Group,
                     data = metadata,
                     permutations = perm_scheme)
print(perm_test)

## 9.2 PCoA ordination -------------------------------------------------------
ord.fun <- cmdscale(tax_singlem_dist, k = 2, eig = TRUE, add = TRUE)
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

## 9.3 Plot ------------------------------------------------------------------
pcoa_singlem_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2,
                                           color = Group,
                                           shape = Group)) +
  stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1,
               show.legend = FALSE, level = 0.95) +
  geom_point(size = 1.5) +
  scale_color_manual(name = "Group",
                     values = group_colors,
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
pcoa_singlem_plot


# 10. Combined figure =========================================================
# a: alpha diversity, b: shared genera, c: PCoA, d: phylum composition
combined_singlem_plot <- ggdraw() +
  draw_plot(div_plot,           x = 0,     y = 0.5, width = 3/5,   height = 0.5) +
  draw_plot(venn_singlem,       x = 0,     y = 0,   width = 1.5/5, height = 0.5) +
  draw_plot(composition_plot,   x = 3/5,   y = 0,   width = 2/5,   height = 1) +
  draw_plot(pcoa_singlem_plot,  x = 1.5/5, y = 0,   width = 1.5/5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 9,
                  x = c(0, 0, 1.5/5, 3/5), y = c(1, 0.5, 0.5, 1))
combined_singlem_plot

# ggsave(file.path(fig_dir, "combined_singlem_plot.pdf"),
#        combined_singlem_plot, width = 183,
#        height = 90, units = "mm", dpi = 600)