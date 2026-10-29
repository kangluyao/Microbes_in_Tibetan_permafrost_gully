# =============================================================================
# 16S rRNA gene community analysis: composition, diversity and beta diversity
# of Collapsed vs. Un-collapsed gullies
# =============================================================================
# Purpose : Describe the phylum-level composition, shared genera, alpha
#           diversity, community composition (PERMANOVA, PCoA) and beta
#           dispersion of the 16S communities in Collapsed and Un-collapsed
#           samples.
# Design  : Samples are nested in gullies (Gully_id), so Gully_id is a random
#           intercept in the GLMMs and a permutation block in PERMANOVA and
#           the dispersion test.
# Inputs  : script/0.read_data.R (sourced)
# Outputs : Figures (ggsave calls are commented out) and printed model tables.
#
# Sections
#   1. Setup (paths, parameters, packages, plot theme)
#   2. Data input
#   3. Helper functions
#   4. Phylum-level community composition
#   5. Alpha diversity (GLMMs)
#   6. Community composition (PERMANOVA and PCoA)
#   7. Beta dispersion
#   8. Combined figure
# =============================================================================

rm(list = ls())  # start from a clean workspace

# 1. Setup ====================================================================

## 1.1 Paths and parameters ---------------------------------------------------
project_dir <- "E:/thermokarst_gully"
setwd(project_dir)  # needed because the source() call below uses a relative path

save.dir   <- file.path(getwd(), "result")  # kept from the original script
result_dir <- file.path(project_dir, "revision/result")  # figures

group_levels <- c("Un-collapsed", "Collapsed")  # reference level first
n_top_taxa   <- 10                              # phyla shown in the bar plot
n_perm       <- 999                             # permutations (PERMANOVA, dispersion)

## 1.2 Packages ---------------------------------------------------------------
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsci)
library(tidyr)        # pivot_longer(), gather(); not loaded in the original
library(tibble)       # rownames_to_column()
library(betapart)     # beta.pair.abund(); not loaded in the original
library(ggvenn)
library(car)          # logit() transform for bounded [0, 1] data
library(glmmTMB)      # generalized linear mixed models
library(DHARMa)       # simulation-based residual diagnostics for GLMMs
library(purrr)        # functional programming / iteration
library(broom.mixed)  # tidy() extraction of model coefficients
library(gghalves)     # geom_half_violin()
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

# Colours of the two groups
group_colors <- c("Un-collapsed" = "#79ceb8", "Collapsed" = "#e95f5c")

# 2. Data input ===============================================================
# Provide `metadata` (Sample_name, Gully_id, Group), `phylo_16s` (phyloseq
# object), `genu` (genus table) and `otu_16s_rare` (rarefied OTU table).
source("script/0.read_data.R")


# 3. Helper functions =========================================================

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


# 4. Phylum-level community composition =======================================

## 4.1 Relative abundance per phylum ------------------------------------------
# Aggregate to phylum level (taxa without a phylum are kept) and convert the
# counts to percent per sample.
subphylo_16s <- tax_glom(phylo_16s, "Phylum", NArm = FALSE)
subphylo_16s.rel <- transform_sample_counts(subphylo_16s,
                                            function(x) 100 * x / sum(x))
ntaxa(subphylo_16s.rel)

ra.tab_16s <- otu_table(subphylo_16s.rel)
sum(ra.tab_16s[, 1])   # check: a sample should sum to 100
subtaxa_tab_16s <- tax_table(subphylo_16s.rel)[, 2]   # phylum names
rownames(ra.tab_16s) <- subtaxa_tab_16s

phyla_df <- ra.tab_16s |>
  as.data.frame() |>
  rownames_to_column("taxonomy") |>
  pivot_longer(-taxonomy, names_to = "Sample",
               values_to = "Relative_abundance") |>
  mutate(Group = case_when(grepl("_C", Sample) ~ "Un-collapsed",
                           grepl("_T", Sample) ~ "Collapsed"),
         Group = factor(Group, levels = group_levels))

# Mean, SD and n of relative abundance per phylum and group
phyla_group_ave_df <- phyla_df |>
  group_by(Group, taxonomy) |>
  summarise(across(Relative_abundance,
                   list(mean = ~mean(., na.rm = TRUE),
                        sd   = ~sd(., na.rm = TRUE),
                        n    = ~n())),
            .groups = "drop")

## 4.2 Top phyla --------------------------------------------------------------
# Top phyla by summed group-mean abundance; "Unassigned" is excluded from the
# ranking (it is then pooled into "Others" in the plot)
top_taxa <- phyla_group_ave_df |>
  group_by(taxonomy) |>
  summarise(total = sum(Relative_abundance_mean)) |>
  filter(taxonomy != "Unassigned") |>
  slice_max(total, n = n_top_taxa) |>
  pull(taxonomy)

plot_df <- phyla_group_ave_df |>
  mutate(taxonomy = ifelse(taxonomy %in% top_taxa, taxonomy, "Others")) |>
  mutate(taxonomy = factor(taxonomy, levels = c(top_taxa, "Others")))

## 4.3 Plot -------------------------------------------------------------------
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
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))
composition_plot


# 5. Alpha diversity ==========================================================

## 5.1 Diversity indices ------------------------------------------------------
# Rarefied OTU table, transposed to samples x OTUs
otu_16s_rare_mat <- as.matrix(t(otu_16s_rare))

# Overview table (not used downstream; the same indices are recomputed in
# `div_table` below together with the sample metadata)
shannon  <- diversity(otu_16s_rare_mat, index = "shannon")
richness <- specnumber(otu_16s_rare_mat)
alpha_div <- data.frame(
  Richness = richness,
  Shannon  = shannon,
  Evenness = shannon / log(richness),
  stringsAsFactors = FALSE
)

div_index <- c("Richness", "Shannon", "Evenness")
comm <- as.data.frame(otu_16s_rare_mat)  # samples x OTUs
# Safeguard: sample order of the community matrix must match the metadata
stopifnot(identical(rownames(comm), metadata$Sample_name))

div_table <- metadata |>
  select(Sample_name, Gully_id, Group) |>
  mutate(Group    = factor(Group, levels = group_levels),
         Richness = specnumber(comm),
         Shannon  = diversity(comm, "shannon"),
         Evenness = Shannon / log(Richness))

## 5.2 GLMMs ------------------------------------------------------------------
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

## 5.3 Plot -------------------------------------------------------------------
# WARNING: the P labels are typed by hand (as in the original script) and are
# not linked to `div_pvals`; check them against the printed table. Order:
# Richness, Shannon, Evenness.
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


# 6. Community composition (PERMANOVA and PCoA) ===============================

## 6.1 Bray-Curtis dissimilarity and PERMANOVA --------------------------------
otu_16s_comm <- as.data.frame(t(otu_16s))
otu_dist <- vegdist(otu_16s_comm, "bray")

# Permutations are restricted within gullies (blocks) to respect the clustered
# design.
set.seed(123)
perm_scheme <- how(blocks = metadata$Gully_id, nperm = n_perm)
perm_test <- adonis2(otu_dist ~ Group,
                     data = metadata,
                     permutations = perm_scheme)
print(perm_test)

## 6.2 PCoA ordination --------------------------------------------------------
ord.fun <- cmdscale(otu_dist, k = 2, eig = TRUE, add = TRUE)
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

## 6.3 Plot -------------------------------------------------------------------
pcoa_16s_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2,
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
pcoa_16s_plot


# 7. Beta dispersion ==========================================================

## 7.1 Bray-Curtis dissimilarity and dispersion -------------------------------
treatment    <- metadata$Group
block        <- metadata$Gully_id

# beta.pair.abund() decomposes abundance-based Bray-Curtis dissimilarity; the
# total dissimilarity (`beta.bray`) is used for the dispersion analysis.
bpair_16s_abund <- beta.pair.abund(otu_16s_comm, index.family = "bray")
dist_16s_total  <- as.dist(bpair_16s_abund$beta.bray)

# Distance of each sample to its group centre (betadisper default: spatial
# median)
bd_16s_total <- betadisper(dist_16s_total, treatment)

## 7.2 Permutation test restricted within gullies -----------------------------
# Treatment labels are shuffled freely, but only within gullies (blocks)
ctrl_how <- how(nperm = n_perm, blocks = block, within = Within(type = "free"))

set.seed(123)
pt_16s_total <- permutest(bd_16s_total, permutations = ctrl_how)
pt_16s_total

## 7.3 Plot -------------------------------------------------------------------
df_16s_disp <- data.frame(
  disp_total = bd_16s_total$distances,
  treatment  = treatment
) |>
  pivot_longer(-treatment, names_to = "disp", values_to = "value") |>
  mutate(treatment = factor(treatment, levels = group_levels))

# One-row data frame for the P label (avoids drawing the label once per row)
disp_label <- tibble(
  npcx  = 0.5,
  npcy  = 0.95,
  label = paste0("P = ", signif(pt_16s_total$tab[1, "Pr(>F)"], 3))
)

disp_16s_plot <- df_16s_disp |>
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
disp_16s_plot


# 8. Combined figure ==========================================================
# a: phylum composition, b: alpha diversity, c: PCoA, d: beta dispersion
combined_16S_plot <- ggdraw() +
  draw_plot(div_plot,           x = 0,     y = 0.5, width = 3/5,   height = 0.5) +
  draw_plot(pcoa_16s_plot,       x = 0,     y = 0,   width = 1.5/5, height = 0.5) +
  draw_plot(composition_plot,   x = 3/5,   y = 0,   width = 2/5,   height = 1) +
  draw_plot(disp_16s_plot,  x = 1.5/5, y = 0,   width = 1.5/5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 9,
                  x = c(0, 0, 1.5/5, 3/5), y = c(1, 0.5, 0.5, 1))
combined_16S_plot

# ggsave(file.path(result_dir, "combined_16S_plot.pdf"),
#        combined_16S_plot, width = 183, height = 90, units = "mm", dpi = 600)
