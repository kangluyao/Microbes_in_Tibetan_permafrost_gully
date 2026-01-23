# ================
# 1) Inputs
# ================
# comm : community table, sites (rows) x species (columns)
# meta: data.frame with columns: sample, block, treatment
library(tidyverse)
comm <- t(singlem_genu)
comm[1:5, 1:5]
meta <- data.frame(sample = metadata$Sample_name,
                   treatment = metadata$Group,
                   block = metadata$Gully_id) %>%
  mutate(treatment = gsub("-", "",  treatment))
head(meta)
# ================
# 2) Parallelly compute DBeta per gully
# ================
library(parallel)
library(foreach)
library(doParallel)

# Set number of parallel cores
n_cores <- detectCores() - 2 # here use 8 cores
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Pre-calculate permutations
library(ecopart)
library(gtools)
SiteComb <- permutations(n=5, r=5, v=1:5) # 5 quarts in each gully and treatment

# Parallel computation to partition β diversity changes into 4 components
blocks <- unique(metadata$Gully_id) # set the gully id as blocks
DBeta_blocks_four.comp <- foreach(g = blocks, 
                                  .packages = c("ecopart", "gtools", "purrr"),
                                  .combine = 'c') %dopar% {
                                    # indices of samples for this block and treatment
                                    idx1 <- which(meta$block == g & meta$treatment == "Uncollapsed")
                                    idx2 <- which(meta$block == g & meta$treatment == "Collapsed")
                                    # check counts
                                    if(length(idx1) < 1 || length(idx2) < 1) return(NULL)
                                    # extract matrices
                                    d1_blockg <- comm[meta$sample[idx1], , drop = FALSE]
                                    d2_blockg <- comm[meta$sample[idx2], , drop = FALSE]
                                    # require same N (here both should be 5)
                                    if(nrow(d1_blockg) != nrow(d2_blockg)) return(NULL)
                                    # Tatsumi's ecopart.multi function with baselga index that returns DBeta
                                    set.seed(123)
                                    ResMatrix <- vapply(1:nrow(SiteComb), function(x) {
                                      ecopart.multi(d1_blockg[SiteComb[x,], ], d2_blockg, index = "baselga")
                                    }, FUN.VALUE = numeric(4)) # 4 components
                                    # output
                                    result <- list(rowMeans(ResMatrix))
                                    names(result) <- as.character(g)
                                    result
                                  } %>% imap_dfr(~ data.frame(
                                    gully = .y,
                                    components = names(.x),
                                    value = as.numeric(.x)
                                  ))


# Determine the individual species contribution to changes in β diversity
DBeta_blocks_sp <- foreach(g = blocks, 
                           .packages = c("ecopart", "gtools", "purrr"),
                           .combine = 'c') %dopar% {
                             # indices of samples for this block and treatment
                             idx1 <- which(meta$block == g & meta$treatment == "Uncollapsed")
                             idx2 <- which(meta$block == g & meta$treatment == "Collapsed")
                             # check counts
                             if(length(idx1) < 1 || length(idx2) < 1) return(NULL)
                             # extract matrices
                             d1_blockg <- comm[meta$sample[idx1], , drop = FALSE]
                             d2_blockg <- comm[meta$sample[idx2], , drop = FALSE]
                             # require same N (here both should be 5)
                             if(nrow(d1_blockg) != nrow(d2_blockg)) return(NULL)
                             # Tatsumi's ecopart.multi function with baselga index that returns DBeta
                             set.seed(123)
                             ResMatrix <- vapply(1:nrow(SiteComb), function(x) {
                               ecopart.multi(d1_blockg[SiteComb[x,], ], d2_blockg, 
                                             index = "baselga", components = "sp")
                             }, FUN.VALUE = numeric(4*ncol(comm))) # this value should be set as 4 (components) * S (species number)
                             # output
                             result <- list(rowMeans(ResMatrix))
                             names(result) <- as.character(g)
                             result
                           } %>% imap_dfr(~ data.frame(gully = .y, # gully
                                                       taxonomy = rep(colnames(comm), each = 4),
                                                       components = rep(c("Subtractive homogenization",
                                                                          "Subtractive differentiation",
                                                                          "Additive homogenization",
                                                                          "Additive differentiation"),
                                                                        times = length(colnames(comm))),
                                                       value = .x # value
                           ))

stopCluster(cl)
# Saving the results
# write.csv(DBeta_blocks_four.comp, file ="E:/thermokarst_gully/revision/result/DBeta_blocks_four.comp.csv")
# write.csv(DBeta_blocks_sp, file ="E:/thermokarst_gully/revision/result/DBeta_blocks_sp.csv")
DBeta_blocks_four.comp[1:5, ]
DBeta_blocks_sp[1:5, ]

# ================
# 3) Aggregate across gullies
# ================
# Calculate mean contribution and standard error for each genera in each component
summary_DBeta_sp_data <- DBeta_blocks_sp %>%
  filter(!value == 0) %>%
  group_by(taxonomy, components) %>%
  summarise(
    mean_contribution = mean(value, na.rm = TRUE),
    se_contribution = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Select top 20 genera with highest absolute mean contribution in each component
top20_data <- summary_DBeta_sp_data %>%
  group_by(components) %>%
  slice_max(order_by = abs(mean_contribution), n = 20) %>%
  mutate(taxonomy = fct_reorder(taxonomy, mean_contribution)) %>%
  ungroup()

# Step 3: Create plot
top20_4process_plot <- ggplot(top20_data, aes(x = interaction(taxonomy, rank(mean_contribution)), y = mean_contribution)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_contribution - se_contribution, 
                    ymax = mean_contribution + se_contribution),
                width = 0.3, linewidth = 0.25) +
  geom_point(size = 1, color = "#2C7BB6", shape = 16) +
  facet_wrap(~ components, scales = "free") +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
  coord_flip() +
  labs(x = "Genera", 
       y = "Contribution (mean ± SE)") +
  main_theme +
  theme(strip.text = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
  )

top20_4process_plot

library(ggpmisc)
summary_DBeta_sp_data  %>%
  # top20_data %>%
  left_join(singlem_genu_group_ave_df %>% 
              filter(Group == "Un-collapsed") %>%
              select(taxonomy,  Relative_abundance_mean), by = "taxonomy") %>%
  ggplot(aes(Relative_abundance_mean, mean_contribution)) +
  geom_point(alpha = 0.6, size = 1.25, color = "#79ceb8") +
  geom_smooth(method = "lm", 
              se = TRUE, color = "black", linewidth = 0.9) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), 
                                 after_stat(p.value.label), sep = "*\", \"*")),
               size = 2.4) +
  facet_grid( ~ components, scales = "free") +
  labs(x = "Mean relative abundance",
       y = "Contribution to beta diversity") +
  theme_bw() +
  main_theme +
  theme(legend.position = "none")

# ggsave(filename = "E:/thermokarst_gully/revision/result/top20_4process_plot.pdf",
#        top20_4process_plot, width = 89, height = 89, units = "mm", dpi = 600)


# ================
# 4) Genomic traits
# ================

# Compare the difference in genera traits
bac_gtdb <- fread("E:/thermokarst_gully/data/metagenome/singlem/bac120_metadata_r226.tsv") %>%
  select(gtdb_taxonomy, ncbi_species_taxid, gc_percentage, 
         genome_size, ncbi_rrna_count, checkm2_completeness)

arc_gtdb <- fread("E:/thermokarst_gully/data/metagenome/singlem/ar53_metadata_r226.tsv") %>%
  select(gtdb_taxonomy, ncbi_species_taxid, gc_percentage, 
         genome_size, ncbi_rrna_count, checkm2_completeness)

gtdb_trait <- rbind(bac_gtdb, arc_gtdb) %>%
  mutate(
    gtdb_taxonomy = str_extract(gtdb_taxonomy, "g__[^;]+"),
  ) %>%
  mutate(across(c(gtdb_taxonomy),  ~ str_replace_all(., ".*__", ""))) %>%
  mutate(
    # Replace any non-numeric entries (e.g., text like 'NA' or 'none') with R's NA
    ncbi_rrna_count = na_if(ncbi_rrna_count, "none"), # Replace specific text
    ncbi_rrna_count = as.numeric(ncbi_rrna_count)
  ) %>% # correct the genome size and rrn copies using completeness
  mutate(genome_size = genome_size/(checkm2_completeness/100),
         ncbi_rrna_count = ncbi_rrna_count/(checkm2_completeness/100))

MAGs_genera <- unique(tax_bin$Genus) %>% gsub("g__", "", .)
singlem_genera <- unique(top20_data$taxonomy)
intersect(MAGs_genera, singlem_genera)
setdiff(MAGs_genera, singlem_genera)
setdiff(singlem_genera, MAGs_genera)



gtdb_trait_top20_genera_ave <- gtdb_trait %>%
  rename("taxonomy" = "gtdb_taxonomy") %>%
  group_by(taxonomy) %>%
  summarise(across(c(gc_percentage, genome_size, 
                     ncbi_rrna_count), mean, na.rm = TRUE)) %>%
  inner_join(top20_data, by = "taxonomy") %>%
  mutate(genome_size = genome_size/1000000)


head(gtdb_trait_top20_genera_ave)


traits_id <- c("gc_percentage", "genome_size", "ncbi_rrna_count")
compo_comparisons <- list(c("Subtractive differentiation", "Subtractive homogenization"),
                          c("Subtractive homogenization", "Additive differentiation"),
                          c("Additive differentiation", "Additive homogenization"),
                          c("Subtractive differentiation", "Additive differentiation"),
                          c("Subtractive homogenization", "Additive homogenization"),
                          c("Subtractive differentiation", "Additive homogenization"))
singlem_traits_compo_comparison <- gtdb_trait_top20_genera_ave %>% 
  select(c("components", all_of(traits_id))) %>%
  pivot_longer(-components, names_to = "traits", values_to = "value") %>%
  mutate(components = factor(components, 
                             levels = c("Subtractive differentiation",
                                        "Subtractive homogenization",
                                        "Additive differentiation", 
                                        "Additive homogenization"))) %>%
  mutate(traits = factor(traits, levels = traits_id)) %>%
  ggplot(aes(components, value)) +
  geom_half_violin(aes(fill = components), color = NA, 
                   position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, alpha = 0.65) +
  # geom_boxplot(aes(fill = component), width = 0.35, 
  #              size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(colour = components), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.35, linewidth = 0.15) +
  stat_compare_means(comparisons = compo_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", 
                     bracket.size = 0.25, bracket.width = 0.1,  size = 2,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~traits, scales = "free_y", ncol = 6) + #
  main_theme +
  theme(legend.position = "none",
        axis.ticks.length=unit(.1, "cm"),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
singlem_traits_compo_comparison

library(patchwork)
combined_beita_par_plot <- disp_plot + 
  top20_4process_plot +
  plot_layout(ncol = 1, heights = c(1, 2.5))
# ggsave(filename = "E:/thermokarst_gully/revision/result/combined_beita_par_plot.pdf",
#        combined_beita_par_plot, width = 165, height = 150, units = "mm", dpi = 600)


top20_4process_traits_plot <- plot_grid(top20_4process_plot, 
                                        singlem_traits_compo_comparison,
                                        ncol = 1, rel_heights = c(2,1))
# ggsave(filename = "E:/thermokarst_gully/revision/result/top20_4process_traits_plot.pdf",
#        top20_4process_traits_plot, width = 150, height = 150, units = "mm", dpi = 600)


genera_traits_bar_plot <- gtdb_trait_top20_genera_ave %>% 
  select(c("taxonomy", all_of(traits_id))) %>%
  pivot_longer(-taxonomy, names_to = "traits", values_to = "value") %>%
  group_by(taxonomy, traits) %>%
  summarise(value = mean(value, na.rm =T), .groups = "drop") %>%
  mutate(traits = factor(traits, levels = traits_id)) %>%
  ggplot(aes(taxonomy, value)) +
  geom_bar(stat = "identity", aes(fill = "#5cc3e8"), color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~traits, scales = "free_y", ncol = 6) + #
  main_theme +
  theme(legend.position = "none",
        axis.ticks.length=unit(.1, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
genera_traits_bar_plot

# ggsave(filename = "E:/thermokarst_gully/revision/result/genera_traits_bar_plot.pdf",
#        genera_traits_bar_plot, width = 89, height = 89, units = "mm", dpi = 600)





df <- gtdb_trait %>%
  rename("taxonomy" = "gtdb_taxonomy") %>%
  group_by(taxonomy) %>%
  summarise(across(c(gc_percentage, genome_size, 
                     ncbi_rrna_count), mean, na.rm = TRUE)) %>%
  inner_join(summary_DBeta_sp_data , by = "taxonomy") %>%
  mutate(genome_size = genome_size/1000000)

# 前20用gtdb_trait_top20_genera_ave
plot_correlation(
  data = df %>%
    filter(components == "Subtractive differentiation"),
  x_var = "genome_size",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = df %>%
    filter(components == "Subtractive differentiation"),
  x_var = "gc_percentage",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = df %>%
    filter(components == "Subtractive differentiation"),
  x_var = "ncbi_rrna_count",
  y_var = "mean_contribution",
  method = "spearman")


plot_correlation(
  data = df %>%
    filter(components == "Additive differentiation"),
  x_var = "genome_size",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = gtdb_trait_top20_genera_ave %>%
    filter(components == "Additive differentiation"),
  x_var = "gc_percentage",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = gtdb_trait_top20_genera_ave %>%
    filter(components == "Additive differentiation"),
  x_var = "ncbi_rrna_count",
  y_var = "mean_contribution",
  method = "spearman")
