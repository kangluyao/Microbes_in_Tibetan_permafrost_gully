# =============================================================================
# SingleM beta diversity: dispersion and partitioning of beta-diversity change
# =============================================================================
# Purpose : (i) Test whether multivariate dispersion (beta dispersion) differs
#           between Collapsed and Un-collapsed samples, and (ii) partition the
#           change in beta diversity (Collapsed vs. Un-collapsed) within each
#           gully into additive/subtractive and differentiation/homogenization
#           components, including the contribution of individual genera.
# Design  : Samples are nested in gullies (Gully_id); permutations are
#           restricted within gullies, and the partitioning is done per gully.
# Inputs  : - script/0.read_data.R, script/0.Functions.R 
#           - SingleM genus table: gully-genus.tsv
#           - Cached partitioning results (CSV) in `result_dir`, read in 5.5
# Outputs : Figures (ggsave calls are commented out); the partitioning results
#           can be exported to CSV (commented out in 5.5).
#
# Sections
#   1. Setup (paths, parameters, packages)
#   2. Helper functions
#   3. Data input
#   4. Beta dispersion (betadisper + block-restricted permutation test)
#   5. Partitioning of beta-diversity change (ecopart, parallel)
#   6. Summary and plots of the components
#   7. Combined figure
#   8. Relative abundance of the top-contributing genera
# =============================================================================

rm(list = ls())  # start from a clean workspace


# 1. Setup ====================================================================

## 1.1 Paths and parameters ---------------------------------------------------
project_dir <- "E:/thermokarst_gully"
setwd(project_dir)  # needed because the source() calls below use relative paths

save.dir    <- file.path(getwd(), "result")  # kept from the original script
singlem_dir <- file.path(project_dir, "data/metagenome/singlem")
result_dir  <- file.path(project_dir, "revision/result")  # CSV and figures

group_levels      <- c("Un-collapsed", "Collapsed")
n_perm            <- 999  # permutations for the dispersion test
n_sites_per_cell  <- 5    # samples per gully and treatment (hard-coded design)
n_reserved_cores  <- 4    # cores left free when running in parallel

# Labels (Unicode escapes keep the script independent of file encoding)
total_change_label <- "Total change (\u0394\u03b2)"          # "Total change (Δβ)"
component_levels   <- c(total_change_label,
                        "Additive component",
                        "Subtractive component",
                        "Additive differentiation",
                        "Additive homogenization",
                        "Subtractive differentiation",
                        "Subtractive homogenization")
# Order of the four components returned per species by ecopart.multi(
# components = "sp"); assumed to match the output of the package.
four_component_labels <- c("Subtractive homogenization",
                           "Subtractive differentiation",
                           "Additive homogenization",
                           "Additive differentiation")

## 1.2 Packages ---------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsci)
library(data.table)
library(gghalves)
library(lme4)
library(betapart)    # beta.pair.abund()
library(parallel)
library(foreach)
library(doParallel)
library(ecopart)     # ecopart.multi()
library(gtools)      # permutations()
library(cowplot)

ecopart_pkgs <- c("ecopart", "gtools", "purrr")  # loaded on each worker

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
group_colors     <- c("Un-collapsed" = "#79ceb8", "Collapsed" = "#e95f5c")
component_colors <- c("#9966cc", "#79ceb8", "#e95f5c",   # total, additive, subtractive
                      "#007fff", "#5cc3e8", "#ffdb00",   # four-component palette:
                      "#21abcd")                         # elements 4-7 (alphabetical)

# 2. Helper functions =========================================================

#' Standard error of the mean.
se <- function(x) {
  sd(x) / sqrt(length(x))
}

#' Partition the beta-diversity change for one gully.
#' The Un-collapsed samples are re-ordered in every row of `site_comb` and
#' ecopart.multi() is run against the Collapsed samples each time; the results
#' are averaged over all orderings so that they do not depend on an arbitrary
#' pairing of samples. Extra arguments (`...`, e.g. `components = "two"`) are
#' passed to ecopart.multi(); `n_out` is the length of its output.
#' Returns a list holding the mean result named by gully, or NULL if the gully
#' lacks either treatment or the two treatments have different sample sizes.
partition_gully_beta <- function(g, comm, meta, site_comb, n_out, ...) {
  idx1 <- which(meta$block == g & meta$treatment == "Uncollapsed")
  idx2 <- which(meta$block == g & meta$treatment == "Collapsed")
  if (length(idx1) < 1 || length(idx2) < 1) return(NULL)
  
  d1_blockg <- comm[meta$sample[idx1], , drop = FALSE]
  d2_blockg <- comm[meta$sample[idx2], , drop = FALSE]
  # Both treatments must have the same number of samples (here 5)
  if (nrow(d1_blockg) != nrow(d2_blockg)) return(NULL)
  
  set.seed(123)
  res_matrix <- vapply(
    seq_len(nrow(site_comb)),
    function(x, ...) {
      ecopart.multi(d1_blockg[site_comb[x, ], ], d2_blockg,
                    index = "baselga", ...)
    },
    FUN.VALUE = numeric(n_out), ...
  )
  
  result <- list(rowMeans(res_matrix))
  names(result) <- as.character(g)
  result
}


# 3. Data input ===============================================================
# Provide `metadata` (Sample_name, Gully_id, Group) and `main_theme`.
source("script/0.read_data.R")

# SingleM genus table: drop unassigned reads and keep the genus name only
singlem_genu <- fread(file.path(singlem_dir, "gully-genus.tsv")) |>
  filter(taxonomy != "unassigned") |>
  mutate(taxonomy = str_extract(taxonomy, "g__[^;]+")) |>
  mutate(taxonomy = str_replace_all(taxonomy, ".*__", "")) |>
  column_to_rownames("taxonomy") |>
  as.data.frame()

singlem_comm <- t(singlem_genu)   # samples x genera
treatment    <- metadata$Group
block        <- metadata$Gully_id


# 4. Beta dispersion ==========================================================

## 4.1 Bray-Curtis dissimilarity and dispersion -------------------------------
# beta.pair.abund() decomposes abundance-based Bray-Curtis dissimilarity; the
# total dissimilarity (`beta.bray`) is used for the dispersion analysis.
bpair_singlem_abund <- beta.pair.abund(singlem_comm, index.family = "bray")
dist_singlem_total  <- as.dist(bpair_singlem_abund$beta.bray)

# Distance of each sample to its group centre (betadisper default: spatial
# median)
bd_singlem_total <- betadisper(dist_singlem_total, treatment)

## 4.2 Permutation test restricted within gullies -----------------------------
# Treatment labels are shuffled freely, but only within gullies (blocks)
ctrl_how <- how(nperm = n_perm, blocks = block, within = Within(type = "free"))

set.seed(123)
pt_singlem_total <- permutest(bd_singlem_total, permutations = ctrl_how)
pt_singlem_total

## 4.3 Plot -------------------------------------------------------------------
df_singlem_disp <- data.frame(
  disp_total = bd_singlem_total$distances,
  treatment  = treatment
) |>
  pivot_longer(-treatment, names_to = "disp", values_to = "value") |>
  mutate(treatment = factor(treatment, levels = group_levels))

# One-row data frame for the P label (avoids drawing the label once per row)
disp_label <- tibble(
  npcx  = 0.5,
  npcy  = 0.95,
  label = paste0("P = ", signif(pt_singlem_total$tab[1, "Pr(>F)"], 3))
)

disp_singlem_plot <- df_singlem_disp |>
  ggplot(aes(treatment, value)) +
  geom_half_violin(aes(fill = treatment), color = NA,
                   position = position_nudge(x = 0.25),
                   side = "r", width = 0.5, alpha = 0.65) +
  geom_boxplot(aes(fill = treatment), width = 0.35,
               size = 0.3, outlier.color = NA, alpha = 0.65) +
  geom_jitter(aes(colour = treatment), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  ggpp::geom_text_npc(data = disp_label,
                      aes(npcx = npcx, npcy = npcy, label = label),
                      size = 2, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Group", y = "\u03b2 dispersion") +   # "β dispersion"
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  main_theme +
  theme(legend.position = "none")
disp_singlem_plot


# 5. Partitioning of beta-diversity change ====================================

## 5.1 Inputs -----------------------------------------------------------------
# `comm`: samples (rows) x genera (columns); `meta`: sample, treatment, block.
# Hyphens are removed from the treatment labels ("Un-collapsed" ->
# "Uncollapsed") to match the labels used in partition_gully_beta().
comm <- singlem_comm
comm[1:5, 1:5]  # quick inspection

meta <- data.frame(sample    = metadata$Sample_name,
                   treatment = metadata$Group,
                   block     = metadata$Gully_id) |>
  mutate(treatment = gsub("-", "", treatment))
head(meta)

## 5.2 Permutation scheme and parallel back end -------------------------------
# All possible orderings of the samples within a gully and treatment
# (5! = 120 rows)
SiteComb <- permutations(n = n_sites_per_cell, r = n_sites_per_cell,
                         v = seq_len(n_sites_per_cell))

n_cores <- detectCores() - n_reserved_cores
cl <- makeCluster(n_cores)
# cl <- makeCluster(n_cores, type = "PSOCK", port = 20000)
registerDoParallel(cl)

blocks <- unique(metadata$Gully_id)  # gullies

## 5.3 Components of the change in beta diversity -----------------------------
# Two components: additive vs. subtractive
DBeta_blocks_two.comp <- (
  foreach(g = blocks, .packages = ecopart_pkgs, .combine = "c") %dopar% {
    partition_gully_beta(g, comm, meta, SiteComb,
                         n_out = 2, components = "two")
  }
) |>
  imap_dfr(~ data.frame(gully      = .y,
                        components = names(.x),
                        value      = as.numeric(.x)))

# Four components (ecopart.multi() default): additive/subtractive x
# differentiation/homogenization
DBeta_blocks_four.comp <- (
  foreach(g = blocks, .packages = ecopart_pkgs, .combine = "c") %dopar% {
    partition_gully_beta(g, comm, meta, SiteComb, n_out = 4)
  }
) |>
  imap_dfr(~ data.frame(gully      = .y,
                        components = names(.x),
                        value      = as.numeric(.x)))

## 5.4 Contribution of individual genera --------------------------------------
# Output length = 4 components x number of genera
DBeta_blocks_sp <- (
  foreach(g = blocks, .packages = ecopart_pkgs, .combine = "c") %dopar% {
    partition_gully_beta(g, comm, meta, SiteComb,
                         n_out = 4 * ncol(comm), components = "sp")
  }
) |>
  imap_dfr(~ data.frame(gully      = .y,
                        taxonomy   = rep(colnames(comm), each = 4),
                        components = rep(four_component_labels,
                                         times = length(colnames(comm))),
                        value      = .x))

stopCluster(cl)

## 5.5 Export and re-import of the results ------------------------------------
# NOTE: the objects above are overwritten by the stored CSV files below, so
# the downstream analysis uses the cached results.
# write.csv(DBeta_blocks_two.comp,
#           file = file.path(result_dir, "singlem_DBeta_blocks_two.comp.csv"))
# write.csv(DBeta_blocks_four.comp,
#           file = file.path(result_dir, "singlem_DBeta_blocks_four.comp.csv"))
# write.csv(DBeta_blocks_sp,
#           file = file.path(result_dir, "singlem_DBeta_blocks_sp.csv"))
DBeta_blocks_two.comp <- read.csv(
  file.path(result_dir, "singlem_DBeta_blocks_two.comp.csv"), row.names = 1)
DBeta_blocks_four.comp <- read.csv(
  file.path(result_dir, "singlem_DBeta_blocks_four.comp.csv"), row.names = 1)
DBeta_blocks_sp <- read.csv(
  file.path(result_dir, "singlem_DBeta_blocks_sp.csv"), row.names = 1)

DBeta_blocks_four.comp[1:5, ]  # quick inspection
DBeta_blocks_sp[1:5, ]


# 6. Summary and plots of the components ======================================

## 6.1 Mean change per component (across gullies) -----------------------------
DBeta_blocks_four.comp |>
  group_by(components) |>
  summarise(mean_value = mean(value, na.rm = TRUE))

# Mean and SE per component, including the total change per gully (sum of the
# four components). The SE of negative means is sign-flipped so that the error
# bar extends away from zero (ymin = mean, ymax = mean + se in the plot).
plot_component_df <- rbind(
  DBeta_blocks_four.comp |>
    select(components, value),
  DBeta_blocks_four.comp |>
    group_by(gully) |>
    summarise(value = sum(value, na.rm = TRUE)) |>
    mutate(components = total_change_label) |>
    select(components, value),
  DBeta_blocks_two.comp |>
    select(components, value)
) |>
  group_by(components) |>
  summarise(mean = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") |>
  mutate(se = case_when(mean < 0 ~ se * (-1),
                        TRUE ~ se),
         components = factor(components, levels = component_levels))

bar_ylim <- c(-0.15, 0.18)

beta_singlem_barplot <- plot_component_df |>
  ggplot(aes(components, mean, fill = components)) +
  geom_bar(position = position_dodge(0.8), stat = "identity",
           colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = bar_ylim) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Components", y = "Changes in \u03b2 diversity") +  # "β diversity"
  scale_fill_manual(values = component_colors) +
  # Shaded backgrounds: total | two-component | four-component groups
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = bar_ylim[1],
           ymax = bar_ylim[2], alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = bar_ylim[1],
           ymax = bar_ylim[2], alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = bar_ylim[1],
           ymax = bar_ylim[2], alpha = 0.1, fill = "#5cc3e8") +
  main_theme +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4))
beta_singlem_barplot

## 6.2 Contribution of individual genera --------------------------------------
# Mean and SE of each genus's contribution across gullies (zero values are
# excluded before averaging, as in the original analysis)
summary_DBeta_sp_data <- DBeta_blocks_sp |>
  filter(value != 0) |>
  group_by(taxonomy, components) |>
  summarise(mean_contribution = mean(value, na.rm = TRUE),
            se_contribution = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# The 20 genera with the largest absolute mean contribution per component
top20_data <- summary_DBeta_sp_data |>
  group_by(components) |>
  slice_max(order_by = abs(mean_contribution), n = 20) |>
  mutate(taxonomy = fct_reorder(taxonomy, mean_contribution)) |>
  ungroup()

# The interaction() with the contribution rank gives each genus a unique x
# position within every free-scale facet; the rank suffix is stripped from the
# axis labels.
top20_4process_singlem_plot <- ggplot(top20_data) +
  geom_segment(aes(x = interaction(taxonomy, rank(mean_contribution)),
                   xend = interaction(taxonomy, rank(mean_contribution)),
                   y = 0, yend = mean_contribution, colour = components),
               linewidth = 0.25) +
  geom_point(aes(x = interaction(taxonomy, rank(mean_contribution)),
                 y = mean_contribution, color = components),
             size = 1.25, shape = 16) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = component_colors[4:7]) +
  facet_wrap(~ components, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
  coord_flip() +
  labs(x = "Genera",
       y = "Contribution (mean \u00b1 SE)") +
  main_theme +
  theme(strip.text = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank())
top20_4process_singlem_plot


# 7. Combined figure ==========================================================
# a: beta dispersion, b: components of change, c: top genera per component
singlem_disp_parti_4.comp_plots <- ggdraw() +
  draw_plot(disp_singlem_plot,
            x = 0, y = 2/3, width = 1/3, height = 1/3) +
  draw_plot(beta_singlem_barplot,
            x = 1/3, y = 2/3, width = 2/3, height = 1/3) +
  draw_plot(top20_4process_singlem_plot,
            x = 0, y = 0, width = 1, height = 2/3) +
  draw_plot_label(label = c("a", "b", "c"), size = 8,
                  x = c(0, 1/3, 0), y = c(1, 1, 2/3))
singlem_disp_parti_4.comp_plots

# ggsave(file.path(result_dir, "singlem_disp_parti_4_11.comp_plots.pdf"),
#        singlem_disp_parti_4.comp_plots,
#        width = 150, height = 155, units = "mm")


# 8. Relative abundance of the top-contributing genera ========================
# Mean, SD and SE of the relative abundance of the top-20 genera of each
# component, in Un-collapsed samples ("_C") only.
top20_tax_control_df <- top20_data |>
  select(taxonomy, components) |>
  left_join(singlem_genu |> rownames_to_column("taxonomy"),
            by = "taxonomy") |>
  select(taxonomy, components,
         all_of(grep("_C", metadata$Sample_name, value = TRUE))) |>
  pivot_longer(-c(taxonomy, components), names_to = "Sample_name",
               values_to = "abundance") |>
  group_by(taxonomy, components) |>
  summarise(across(abundance, list(mean = ~ mean(.x, na.rm = TRUE),
                                   sd   = ~ sd(.x, na.rm = TRUE),
                                   se   = ~ se(.x)))) |>
  left_join(top20_data, by = c("taxonomy", "components"))

top20_tax_abun_plot <- ggplot(top20_tax_control_df,
                              aes(x = interaction(taxonomy,
                                                  rank(mean_contribution)),
                                  y = abundance_mean,
                                  fill = components)) +
  geom_bar(position = position_dodge(0.8), stat = "identity",
           colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = abundance_mean,
                    ymax = abundance_mean + abundance_se), width = 0.2) +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$"),
                   limits = rev) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Components", y = "Relative abundance (%)") +
  scale_fill_manual(values = component_colors[4:7]) +
  facet_wrap(~ components, scales = "free_x") +
  main_theme +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4))
top20_tax_abun_plot

# ggsave(file.path(result_dir, "top20_tax_abun_plot.pdf"),
#        top20_tax_abun_plot, width = 120, height = 98, units = "mm")