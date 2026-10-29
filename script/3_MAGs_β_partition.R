# =============================================================================
# MAG-based beta diversity: dispersion, PERMANOVA/PCoA and partitioning of
# beta-diversity change (Collapsed vs. Un-collapsed gullies)
# =============================================================================
# Purpose : (i) Test whether multivariate dispersion (beta dispersion) differs
#           between Collapsed and Un-collapsed samples, (ii) test and
#           visualise differences in community composition (PERMANOVA, PCoA),
#           and (iii) partition the change in beta diversity within each gully
#           into additive/subtractive and differentiation/homogenization
#           components, including the contribution of individual MAGs.
# Design  : Samples are nested in gullies (Gully_id); permutations are
#           restricted within gullies, and the partitioning is done per gully.
# Inputs  : - script/0.read_data.R, script/0.Functions.R)
#           - Cached partitioning results (CSV) in `result_dir`, read in 6.5
# Outputs : Figures (ggsave calls are commented out); the partitioning results
#           can be exported to CSV (commented out in 6.5).
#
# Sections
#   1. Setup (paths, parameters, packages)
#   2. Helper functions
#   3. Data input
#   4. Beta dispersion (betadisper + block-restricted permutation test)
#   5. PERMANOVA and PCoA (Hellinger-transformed data)
#   6. Partitioning of beta-diversity change (ecopart, parallel)
#   7. Summary and plots of the components
# =============================================================================

rm(list = ls())  # start from a clean workspace


# 1. Setup ====================================================================

## 1.1 Paths and parameters ---------------------------------------------------
project_dir <- "E:/thermokarst_gully"
setwd(project_dir)  # needed because the source() calls below use relative paths
result_dir  <- file.path(project_dir, "revision/result")  # CSV and figures

group_levels     <- c("Un-collapsed", "Collapsed")  # reference level first
n_sites_per_cell <- 5   # samples per gully and treatment (hard-coded design)
n_reserved_cores <- 4   # cores left free when running in parallel
cluster_port     <- 20000

# Colours
group_colors     <- c("Un-collapsed" = "#79ceb8", "Collapsed" = "#e95f5c")
component_colors <- c("#9966cc", "#79ceb8", "#e95f5c",   # total, additive, subtractive
                      "#007fff", "#5cc3e8", "#ffdb00",   # four-component palette:
                      "#21abcd")                         # elements 4-7 (alphabetical)

# Labels (Unicode escapes keep the script independent of file encoding)
total_change_label <- "Total change (\u0394\u03b2)"          # "Total change (Δβ)"
component_levels   <- c(total_change_label,
                        "Additive component",
                        "Subtractive component",
                        "Additive differentiation",
                        "Additive homogenization",
                        "Subtractive differentiation",
                        "Subtractive homogenization")
# Order of the four components returned per MAG by ecopart.multi(
# components = "sp"); assumed to match the output of the package.
four_component_labels <- c("Subtractive homogenization",
                           "Subtractive differentiation",
                           "Additive homogenization",
                           "Additive differentiation")

## 1.2 Packages ---------------------------------------------------------------
library(betapart)    # beta.pair.abund()
library(vegan)
library(tidyverse)
library(gghalves)    # geom_half_violin(); not loaded in the original script
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

# 2. Helper functions =========================================================

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
# Provide `metadata` (Sample_name, Gully_id, Group), `mags_abun_tab`,
# `mags_tax` (columns ID, Genus) and `main_theme`.
source("script/0.read_data.R")

mag_comm  <- t(mags_abun_tab)   # samples x MAGs
treatment <- metadata$Group
block     <- metadata$Gully_id


# 4. Beta dispersion ==========================================================

## 4.1 Bray-Curtis dissimilarity and dispersion -------------------------------
# beta.pair.abund() decomposes abundance-based Bray-Curtis dissimilarity; the
# total dissimilarity (`beta.bray`) is used for the dispersion analysis.
bpair_mag_abund <- beta.pair.abund(mag_comm, index.family = "bray")
dist_mag_total  <- as.dist(bpair_mag_abund$beta.bray)

# Distance of each sample to its group centre (betadisper default: spatial
# median)
bd_mag_total <- betadisper(dist_mag_total, treatment)

## 4.2 Permutation test restricted within gullies -----------------------------
# Treatment labels are shuffled freely, but only within gullies (blocks)
ctrl_how <- how(nperm = 999, blocks = block, within = Within(type = "free"))

pt_mag_total <- permutest(bd_mag_total, permutations = ctrl_how)
pt_mag_total

## 4.3 Plot -------------------------------------------------------------------
df_mag_disp <- data.frame(
  disp_total = bd_mag_total$distances,
  treatment  = treatment
) |>
  pivot_longer(-treatment, names_to = "disp", values_to = "value") |>
  mutate(treatment = factor(treatment, levels = group_levels))

# One-row data frame for the P label (avoids drawing the label once per row)
disp_label <- tibble(
  npcx  = 0.5,
  npcy  = 0.95,
  label = paste0("P = ", signif(pt_mag_total$tab[1, "Pr(>F)"], 3))
)

disp_mag_plot <- df_mag_disp |>
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
disp_mag_plot


# 5. PERMANOVA and PCoA =======================================================

## 5.1 Dissimilarity and PERMANOVA --------------------------------------------
# Hellinger transformation reduces the influence of highly abundant MAGs
# before computing Bray-Curtis dissimilarity.
mag_results_trans <- decostand(mag_comm, method = "hellinger")
mag_dist <- vegdist(mag_results_trans, "bray")

# PERMANOVA on compositional differences; permutations restricted within
# gullies (blocks)
perm_scheme <- how(blocks = metadata$Gully_id, nperm = 999)
set.seed(123)
perm_test <- adonis2(mag_dist ~ Group,
                     data = metadata,
                     permutations = perm_scheme)
print(perm_test)

## 5.2 PCoA ordination --------------------------------------------------------
ord.fun <- cmdscale(mag_dist, k = 2, eig = TRUE, add = TRUE)
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

## 5.3 Plot -------------------------------------------------------------------
pcoa_mag_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2,
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
                         "p =", round(p_value, 3)),
           hjust = 1.05, vjust = 1.1, size = 2) +
  main_theme
pcoa_mag_plot


# 6. Partitioning of beta-diversity change ====================================

## 6.1 Inputs -----------------------------------------------------------------
# `comm`: samples (rows) x MAGs (columns); `meta`: sample, treatment, block.
# Hyphens are removed from the treatment labels ("Un-collapsed" ->
# "Uncollapsed") to match the labels used in partition_gully_beta().
comm <- mag_comm
comm[1:5, 1:5]  # quick inspection

meta <- data.frame(sample    = metadata$Sample_name,
                   treatment = metadata$Group,
                   block     = metadata$Gully_id) |>
  mutate(treatment = gsub("-", "", treatment))
head(meta)

## 6.2 Permutation scheme and parallel back end -------------------------------
# All possible orderings of the samples within a gully and treatment
# (5! = 120 rows)
SiteComb <- permutations(n = n_sites_per_cell, r = n_sites_per_cell,
                         v = seq_len(n_sites_per_cell))

n_cores <- detectCores() - n_reserved_cores
# cl <- makeCluster(n_cores)  # redundant: replaced by the call below, and the
#                             # first cluster was never stopped in the original
cl <- makeCluster(n_cores, type = "PSOCK", port = cluster_port)
registerDoParallel(cl)

blocks <- unique(metadata$Gully_id)  # gullies

## 6.3 Components of the change in beta diversity -----------------------------
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

## 6.4 Contribution of individual MAGs ----------------------------------------
# Output length = 4 components x number of MAGs
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

## 6.5 Export and re-import of the results ------------------------------------
# NOTE: the objects above are overwritten by the stored CSV files below, so
# the downstream analysis uses the cached results.
# write.csv(DBeta_blocks_two.comp,
#           file = file.path(result_dir, "MAGs_DBeta_blocks_two.comp.csv"))
# write.csv(DBeta_blocks_four.comp,
#           file = file.path(result_dir, "MAGs_DBeta_blocks_four.comp.csv"))
# write.csv(DBeta_blocks_sp,
#           file = file.path(result_dir, "MAGs_DBeta_blocks_sp.csv"))
DBeta_blocks_two.comp <- read.csv(
  file.path(result_dir, "MAGs_DBeta_blocks_two.comp.csv"), row.names = 1)
DBeta_blocks_four.comp <- read.csv(
  file.path(result_dir, "MAGs_DBeta_blocks_four.comp.csv"), row.names = 1)
DBeta_blocks_sp <- read.csv(
  file.path(result_dir, "MAGs_DBeta_blocks_sp.csv"), row.names = 1)

DBeta_blocks_four.comp[1:5, ]  # quick inspection
DBeta_blocks_sp[1:5, ]


# 7. Summary and plots of the components ======================================

## 7.1 Mean change per component (across gullies) -----------------------------
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

bar_ylim <- c(-0.23, 0.3)

beta_mag_barplot <- plot_component_df |>
  ggplot(aes(components, mean, fill = components)) +
  geom_bar(position = position_dodge(0.5), stat = "identity",
           colour = "black", width = 0.45) +
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
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(color = "black", linewidth = 0.4))
beta_mag_barplot

## 7.2 Combined figure --------------------------------------------------------
# a: PCoA, b: beta dispersion, c: components of change
all_beta_mag_plot <- ggdraw() +
  draw_plot(pcoa_mag_plot,
            x = 0, y = 3/5, width = 1/2, height = 2/5) +
  draw_plot(disp_mag_plot,
            x = 1/2, y = 3/5, width = 1/2, height = 2/5) +
  draw_plot(beta_mag_barplot,
            x = 0, y = 0, width = 1, height = 3/5) +
  draw_plot_label(label = c("a", "b", "c"), size = 9,
                  x = c(0, 1/2, 0), y = c(1, 1, 3/5))
all_beta_mag_plot

# ggsave(file.path(result_dir, "all_\u03b2_mag_plot.pdf"),  # "all_β_mag_plot.pdf"
#        all_beta_mag_plot, width = 135, height = 120, units = "mm", dpi = 600)

## 7.3 Contribution of individual MAGs ----------------------------------------
# Mean and SE of each MAG's contribution across gullies (zero values are
# excluded before averaging, as in the original analysis)
summary_DBeta_sp_mags_data <- DBeta_blocks_sp |>
  filter(value != 0) |>
  group_by(taxonomy, components) |>
  summarise(mean_contribution = mean(value, na.rm = TRUE),
            se_contribution = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# The 20 MAGs with the largest absolute mean contribution per component
top20_mags_data <- summary_DBeta_sp_mags_data |>
  group_by(components) |>
  slice_max(order_by = abs(mean_contribution), n = 20) |>
  mutate(taxonomy = fct_reorder(taxonomy, mean_contribution)) |>
  ungroup()

# Add the genus to each MAG label, e.g. "<Genome_id> (<Genus>)"
top20_mags_plot_data <- top20_mags_data |>
  dplyr::rename(Genome_id = taxonomy) |>
  left_join(mags_tax |> dplyr::select(Genome_id, Genus), by = "Genome_id") |>
  mutate(Genus = gsub("g__", "", Genus),
         taxonomy = paste0(Genome_id, " (", Genus, ")"))

# The interaction() with the contribution rank gives each MAG a unique x
# position within every free-scale facet; the rank suffix is stripped from the
# axis labels.
top20_4process_amgs_plot <- ggplot(top20_mags_plot_data) +
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
top20_4process_amgs_plot

# ggsave(file.path(result_dir, "top20_4process_amgs_plot.pdf"),
#        top20_4process_amgs_plot, width = 120, height = 155, units = "mm")