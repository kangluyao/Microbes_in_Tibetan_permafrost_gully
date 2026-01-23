source("script/read_data.R")
source("script/Functions.R")
# --- β dispersion analysis ---
mag_comm <- t(mags_abun_tab)
treatment <- metadata$Group
block <- metadata$Gully_id

# --- 1) Beta-diversity decomposition ---
# For abundance data
library(betapart)
library(vegan)
bpair_mag_abund <- beta.pair.abund(mag_comm, index.family = "bray")

# Convert to dist objects for betadisper()
dist_mag_total <- as.dist(bpair_mag_abund$beta.bray)
dist_mag_turn  <- as.dist(bpair_mag_abund$beta.bray.bal)
dist_mag_nest  <- as.dist(bpair_mag_abund$beta.bray.gra)

# --- 2) Calculate dispersion (distance to group centroid) for each treatment ---
bd_mag_total <- betadisper(dist_mag_total, treatment)
bd_mag_turn  <- betadisper(dist_mag_turn, treatment)
bd_mag_nest  <- betadisper(dist_mag_nest, treatment)

# --- 3) Permutation test with block-constrained randomization ---
# Create permutation control object
# Permutations are restricted within each block (so that treatment labels are shuffled only within blocks)
ctrl_how <- how(nperm = 999, blocks = block, within = Within(type = "free"))

# Run permutational tests
pt_mag_total <- permutest(bd_mag_total, permutations = ctrl_how)
pt_mag_turn  <- permutest(bd_mag_turn, permutations = ctrl_how)
pt_mag_nest  <- permutest(bd_mag_nest, permutations = ctrl_how)

# Inspect results
pt_mag_total
pt_mag_turn
pt_mag_nest

# --- 4) Visualization: PCoA and boxplots of sample dispersion values ---
# bin genome information
mag_results_trans <- decostand(mag_comm, method = "hellinger")
mag_dist <-vegdist(mag_results_trans, "bray")
#permanova test the difference in compositional variance
adonis2(mag_dist ~ Group, data = metadata)

#Visualization for the overall difference by PCoA plot.
#mag
pcoa_mag_plot <- PCoA_plot_fun(mag_dist)

# Boxplot for beta dispersion
df_mag_disp <- data.frame(
  sample = names(bd_mag_total$distances),
  disp_total = bd_mag_total$distances,
  disp_turn  = bd_mag_turn$distances,
  disp_nest  = bd_mag_nest$distances,
  treatment = treatment,
  block = block
) %>% select(disp_total, disp_turn, disp_nest, treatment) %>%
  pivot_longer(-treatment, names_to = "disp", values_to = "value") %>%
  mutate(treatment = factor(treatment, levels = c("Un-collapsed", "Collapsed"))) 

# Plot dispersion (divergence) distributions
library(gghalves)
disp_index <- c("disp_total", "disp_turn", "disp_nest")
disp_mag_plot <- df_mag_disp %>%
  filter(disp == "disp_total") %>%
  mutate(disp = factor(disp, levels = disp_index)) %>%
  ggplot(aes(treatment, value)) +
  geom_half_violin(aes(fill = treatment), color = NA, 
                   position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, alpha = 0.65) +
  geom_boxplot(aes(fill = treatment), width = 0.35, 
               size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(colour = treatment), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  # ggpp::geom_text_npc(data = sig.div.labs, aes(npcx = x1, npcy = y1, 
  #                                              label = sig.labels), 
  #                     size = 2, inherit.aes = F) +
  labs(x = "Group", y = "β dispersion") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  # facet_wrap(~disp, scales = "free_y", ncol = 3,
  #            labeller = as_labeller(c(
  #              `disp_total` = "Total dispersion",
  #              `disp_turn` = "Turnover dispersion",
  #              `disp_nest` = "Nestedness dispersion"
  #            ))) +
  main_theme +
  theme(legend.position = "none")

# --- β partition analysis ---
# ================
# 1) Inputs
# ================
# comm : community table, sites (rows) x species (columns)
# meta: data.frame with columns: sample, block, treatment
library(tidyverse)
comm <- mag_comm
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

# Pre-calculate permutations
library(ecopart)
library(gtools)
SiteComb <- permutations(n=5, r=5, v=1:5) # 5 quarts in each gully and treatment

# Set number of parallel cores
n_cores <- detectCores() - 4 # here use 8 cores
cl <- makeCluster(n_cores)
cl <- makeCluster(n_cores, type = "PSOCK", port = 20000)
registerDoParallel(cl)
# Parallel computation to partition β diversity changes into 4 components
blocks <- unique(metadata$Gully_id) # set the gully id as blocks
DBeta_blocks_two.comp <- foreach(g = blocks, 
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
                                                    index = "baselga", components = "two")
                                    }, FUN.VALUE = numeric(2)) # 2 components
                                    # output
                                    result <- list(rowMeans(ResMatrix))
                                    names(result) <- as.character(g)
                                    result
                                  } %>% imap_dfr(~ data.frame(
                                    gully = .y,
                                    components = names(.x),
                                    value = as.numeric(.x)
                                  ))

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
                                      ecopart.multi(d1_blockg[SiteComb[x,], ], d2_blockg, 
                                                    index = "baselga")
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
# write.csv(DBeta_blocks_two.comp, file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_two.comp.csv")
# write.csv(DBeta_blocks_four.comp, file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_four.comp.csv")
# write.csv(DBeta_blocks_sp, file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_sp.csv")
DBeta_blocks_two.comp <- read.csv(file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_two.comp.csv",
                                  row.names = 1)
DBeta_blocks_four.comp <- read.csv(file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_four.comp.csv",
                                   row.names = 1)
DBeta_blocks_sp <- read.csv(file ="E:/thermokarst_gully/revision/result/MAGs_DBeta_blocks_sp.csv",
                            row.names = 1)
DBeta_blocks_four.comp[1:5, ]
DBeta_blocks_sp[1:5, ]

# ================
# 3) Aggregate across gullies
# ================
# Calculate mean contribution and standard error for each genera in each component
DBeta_blocks_four.comp %>%
  group_by(components) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

plot_component_df <- rbind(DBeta_blocks_four.comp %>%
                             select(components, value),
                           DBeta_blocks_four.comp %>%
                             group_by(gully) %>%
                             summarise(value = sum(value, na.rm = TRUE)) %>%
                             mutate(components = rep("Total change (Δβ)", 6)) %>%
                             select(components, value),
                           DBeta_blocks_two.comp %>%
                             select(components, value)) %>%
  group_by(components) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(se = case_when(mean < 0 ~ se * (-1),
                        TRUE ~ se),
         components = factor(components, 
                             levels = c("Total change (Δβ)",
                                        "Additive component",
                                        "Subtractive component",
                                        "Additive differentiation",
                                        "Additive homogenization",
                                        "Subtractive differentiation",
                                        "Subtractive homogenization")))
β_mag_barplot <- plot_component_df %>%
  ggplot(aes(components, mean, fill = components)) +
  geom_bar(position = position_dodge(0.5), stat = 'identity', 
           colour = 'black', width = 0.45, )+
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = .2) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.23, 0.3)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Components", y = "Changes in β diversity") +
  scale_fill_manual(values = c("#9966cc", "#79ceb8", "#e95f5c", 
                               "#007fff", "#5cc3e8", "#ffdb00",
                               "#21abcd")) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -0.23, ymax = 0.3, alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = -0.23, ymax = 0.3, alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = -0.23, ymax = 0.3, alpha = 0.1, fill = "#5cc3e8") +
  main_theme + 
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(color = "black", linewidth = 0.4))
β_mag_barplot

# Combine the PCoA plot, boxplot for beta dispersion, and beta partition
library(cowplot)
all_β_mag_plot <-  ggdraw() +
  draw_plot(pcoa_mag_plot, x = 0, y = 3/5, width = 1/2, height = 2/5) +
  draw_plot(disp_mag_plot, x = 1/2, y = 3/5, width = 1/2, height = 2/5) +
  draw_plot(β_mag_barplot, x = 0, y = 0, width = 1, height = 3/5) +
  draw_plot_label(label = c("a", "b", "c"), size = 9,
                  x = c(0, 1/2, 0), y = c(1, 1, 3/5))
all_β_mag_plot 
# ggsave("E:/thermokarst_gully/revision/result/all_β_mag_plot.pdf",
#        all_β_mag_plot, width = 135, height = 120, units = "mm", dpi = 600)


summary_DBeta_sp_mags_data <- DBeta_blocks_sp %>%
  filter(!value == 0) %>%
  group_by(taxonomy, components) %>%
  summarise(
    mean_contribution = mean(value, na.rm = TRUE),
    se_contribution = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Select top 20 genera with highest absolute mean contribution in each component
top20_mags_data <- summary_DBeta_sp_mags_data %>%
  group_by(components) %>%
  slice_max(order_by = abs(mean_contribution), n = 20) %>%
  mutate(taxonomy = fct_reorder(taxonomy, mean_contribution)) %>%
  ungroup()

library(dplyr)

cum_70_df_filtered <- summary_DBeta_sp_mags_data %>%
  group_by(components) %>%
  arrange(desc(mean_contribution)) %>%
  mutate(
    cumulative_contribution = cumsum(mean_contribution)
  ) %>%
  filter(cumulative_contribution <= 0.7) %>%
  # select(-cumulative_contribution) %>%
  ungroup()


# --- create plot ---
df_plot <- summary_DBeta_sp_mags_data %>%
  group_by(components) %>%
  mutate(
    total_contribution = sum(mean_contribution, na.rm = TRUE), 
    relative_contribution = mean_contribution / total_contribution 
  ) %>%
  arrange(desc(relative_contribution)) %>%
  mutate(
    cumulative_contribution = case_when(components %in% c("Additive homogenization", "Subtractive homogenization")  ~ cumsum(relative_contribution) * (-1),
                                                       TRUE ~ cumsum(relative_contribution)),
    x_rank = row_number() 
  ) %>%
  ungroup()

# Determine the scale factor
max_left_y <- max(df_plot$mean_contribution, na.rm = TRUE)
max_right_y <- max(df_plot$cumulative_contribution, na.rm = TRUE) 
scaling_factor <- max_right_y / max_left_y
head(df_plot)

parti_4_comp_plot <- ggplot(df_plot, aes(x = x_rank)) +
  
  # --- 3.1. 左轴：relative_contribution (柱状图) ---
  geom_col(
    aes(y = mean_contribution, fill = components), 
    position = "dodge", 
    alpha = 0.7
  ) +
  
  # --- 3.2. 右轴：cumulative_contribution (折线图和点) ---
  geom_line(
    aes(y = cumulative_contribution / scaling_factor, group = components),
    color = "red", 
    size = 0.75
  ) +
  geom_point(
    aes(y = cumulative_contribution / scaling_factor, group = components),
    color = "red", 
    size = 0.5
  ) +
  
  # --- 3.3. 设置 X 轴：使用 taxonomy 标签替换行号 ---
  scale_x_continuous(
    # 确保 breaks 和 labels 对应正确
    expand = expansion(mult = c(0, 0)),
    breaks = df_plot$x_rank[order(df_plot$components, df_plot$x_rank)],
    labels = df_plot$taxonomy[order(df_plot$components, df_plot$x_rank)],
    name = "Taxonomy"
  ) +
  
  # --- 3.4. 设置 Y 轴标签和双坐标轴 ---
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0.1)),
    name = "Mean contribution (Individual)", # 左轴名称已更新
    # labels = scales::percent, 
    sec.axis = sec_axis(
      ~ . * scaling_factor, 
      name = "Cumulative contribution", # 右轴名称已更新
      labels = scales::percent 
    )
  ) +
  
  # --- 3.5. 分组绘图 ---
  facet_wrap(~ components, scales = "free") + 
  
  # --- 3.6. 其它美化 ---
  main_theme +
  theme(axis.text.x = element_blank(), 
    legend.position = "none"
  )

# Combine all plots
library(cowplot)
disp_mags_parti_plot <- plot_grid(disp_mag_plot, β_mag_barplot,
                                  nrow = 1, align = "hv",
                                  rel_widths = c(1, 2))
mags_disp_parti_4.comp_plots <- plot_grid(disp_mags_parti_plot, parti_4_comp_plot,
                        ncol = 1, align = "h", axis = "lr",
                        rel_heights = c(1.5, 2))
mags_disp_parti_4.comp_plots

# Save the combined plot
# ggsave(file.path("E:/thermokarst_gully/revision/result/mags_disp_parti_4.comp_plots.pdf"),
#        mags_disp_parti_4.comp_plots, width = 120, height = 155, units = "mm")

# Step 3: Create plot
top20_mags_plot_data <- top20_mags_data %>% rename("ID" = "taxonomy") %>%
  left_join(mags_tax %>% select(ID, Genus), by = "ID") %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%
  mutate(taxonomy = paste0(ID, " (", Genus, ")"))

top20_4process_amgs_plot <- ggplot(top20_mags_plot_data) +
  geom_segment(aes(x = interaction(taxonomy, rank(mean_contribution)),
                   xend = interaction(taxonomy, rank(mean_contribution)), 
                   y = 0, yend = mean_contribution, colour = components), 
               linewidth = 0.25) +
  geom_point(aes(x = interaction(taxonomy, rank(mean_contribution)), 
                 y = mean_contribution, color = components), 
             size = 1.25, shape = 16) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("#007fff", "#5cc3e8", "#ffdb00",
                                "#21abcd")) +
  facet_wrap(~ components, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
  coord_flip() +
  labs(x = "Genera", 
       y = "Contribution (mean ± SE)") +
  main_theme +
  theme(strip.text = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank()
  )

top20_4process_amgs_plot

# ggsave(filename = "E:/thermokarst_gully/revision/result/top20_4process_plot.pdf",
#        top20_4process_plot, width = 89, height = 89, units = "mm", dpi = 600)



  

# top20_4process_amgs_plot <- ggplot(top20_data, aes(x = interaction(taxonomy, rank(mean_contribution)), y = mean_contribution)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_errorbar(aes(ymin = mean_contribution - se_contribution, 
#                     ymax = mean_contribution + se_contribution),
#                 width = 0.3, linewidth = 0.25) +
#   geom_point(size = 1, color = "#2C7BB6", shape = 16) +
#   facet_wrap(~ components, scales = "free") +
#   scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
#   coord_flip() +
#   labs(x = "Genera", 
#        y = "Contribution (mean ± SE)") +
#   main_theme +
#   theme(strip.text = element_text(face = "bold", size = 6),
#         axis.text = element_text(size = 6),
#         plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
#         panel.grid.minor = element_blank()
#   )
# 
# top20_4process_amgs_plot

library(ggpmisc)
summary_DBeta_sp_mags_data  %>%
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

