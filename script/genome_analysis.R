wd = "E:/thermokarst_gully/data/metagenome/MAGs"
save.dir <- "E:/thermokarst_gully/result"
mags.att.file <- "Network_Cytoscape_sif_0.37_RMT_2024_10_18_sif_default_node.txt"
abun.file <- "bin_abundance_coverm.txt"

#loading packages
library(tidyverse)
setwd(wd)
mags.att <- read.table(mags.att.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE)
abun <- read.table(abun.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)[-1, ]

dat2 <- abun %>% 
  mutate(Control = rowMeans(select(., grep('_C', colnames(abun), value = T)))) %>%
  mutate(Collapsed = rowMeans(select(., grep('_T', colnames(abun), value = T)))) %>%
  rownames_to_column("ID") %>%
  select(c(ID, Control, Collapsed)) %>%
  pivot_longer(cols = -c(ID), names_to = "Group", values_to = 'rel_abun') %>%
  mutate(Group = factor(Group, levels = c('Control', 'Collapsed')))
# write the table for itol annotation
write.csv(dat2, "E:/thermokarst_gully/result/MAGs/itol/abundance_annotation_coverm.csv")


#ZP plot
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 7),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 7),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.title = element_text(colour = 'black', size = 7),
        legend.text = element_text(colour = 'black', size = 6),
        legend.key.size = unit(0.5, 'cm'))
ZP_plot <- mags.att %>%
  mutate(`Network roles` = factor(`Network roles`, levels = c("Network hubs", "Module hubs",
                                                              "Connector hubs", "Peripheral species"))) %>%
  ggplot(aes(x = Pi, y = Zi)) +
  geom_point(shape = 19, alpha = 0.75, aes(colour = `Network roles`)) +
  scale_color_manual(values = c("#e95f5c", "#79ceb8", "#5cc3e8", "#ffdb00")) +
  scale_x_continuous(limits = c(0, 1.0)) +
  scale_y_continuous(limits = c(-2, 7)) +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5) +
  labs(x = "Among−module connectivity (Pi)",
       y = "Within−module connectivity (Zi)") +
  main_theme +
  theme(legend.position = c(0.2, 0.8))
# save the plot
# ggsave(file.path(save.dir, "./figs/MAGs/ZP_plot.pdf"),
#        ZP_plot, width = 90, height = 90, units = "mm")
ZP_plot

#difference in the core genomes relative abundance between the uncollapsed and collapsed samples
net_hubs <- mags.att %>% as_tibble() %>% filter(`Network roles` == "Network hubs") %>% pull(`shared name`)
net_cor <- abun[net_hubs, ] %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Sample_id") %>%
  mutate(Gully_id = gsub("_.+$", "", Sample_id)) %>%
  mutate(Group = case_when(grepl("_C", Sample_id) ~ "Uncollapsed",
                           grepl("_T", Sample_id) ~ "Collapsed"))

library(lme4)
library(lmerTest)
lmm_mags_modes <- lapply(net_hubs, function(x) {
  lmer(substitute(i ~ Group + (1|Gully_id), list(i = as.name(x))), data = net_cor)})
summary.model <- function(model){
  F.value <- anova(model)$'F value'
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
df <- NULL
for(i in 1:length(net_hubs)) {
  tmp <- summary.model(lmm_mags_modes[[i]])
  if (is.null(df)){
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}

mags_result_lmm <- data.frame(Variables = net_hubs, group1 = rep("Uncollapsed", length(net_hubs)),
                            group2 = rep("Collapsed", length(net_hubs)), df)
mags_result_lmm
  

library(gghalves)
annot_text <- tibble(Genome = factor(mags_result_lmm$Variables, levels = net_hubs),
                x1 = rep(0.5, length(net_hubs)),
                y1 = rep(0.95, length(net_hubs)),
                sig.labels = mags_result_lmm$sig,
                x2 = rep(0.1, length(net_hubs)),
                y2 = rep(1, length(net_hubs)),
                panel.labels = letters[as.numeric(Genome)])


mags_abun_plot <- net_cor %>% select(c("Group", net_hubs)) %>%
  gather(Genome, value, -c("Group")) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  mutate(Genome = factor(Genome, levels = net_hubs)) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = Group), shape = 21, size = 1, width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = "Relative abundance (%)") +
  # add a sigment to denote the comparison between groups
  ggpp::geom_text_npc(data = annot_text, aes(npcx = x1, npcy = y1, label = sig.labels), inherit.aes = F) +
  ggpp::geom_text_npc(data = annot_text, aes(npcx = x2, npcy = y2, label = panel.labels), inherit.aes = F) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~Genome, scales = "free_y", ncol = 5) +
  main_theme +
  theme(legend.position = "none")


# save the plot
# ggsave(file.path(save.dir, "./figs/MAGs/mags_abun_plot.pdf"),
#        mags_abun_plot, width = 7.59, height = 2.74, units = "in")
mags_abun_plot



source("e:/thermokarst_gully/script/read_data_all.R")
# bin genome information
metabolic_output_file <- file.path(wd, "METABOLIC_result.csv")
metabolic_tab <- read.csv(metabolic_output_file) %>%
  select(c('Category', 'Function', grep('Hmm.presence', colnames(.))))

target_cat <- c("Thermophilic specific", "Ethanol fermentation", "Complex carbon degradation", "Fermentation", "C1 metabolism", "Carbon fixation", "Methane metabolism", "Amino acid utilization", "Fatty acid degradation", "Aromatics degradation",  "Nitrogen cycling", "Sulfur cycling", "Hydrogenases", "Iron cycling", "Manganese cycling")

C_N_metapathway <- c("Amino acid utilization", "Amylolytic enzymes", 
                     "Cellulose degrading", "Endohemicellulases",
                     "Hemicullulose debranching", "Chitin degrading", 
                     "Phenol => Benzoyl-CoA", "acyl-CoA dehydrogenase",
                     "Other oligosaccharide degrading",
                     "Acetogenesis", "Lactate utilization", "Pyruvate oxidation",  
                     "Pyruvate <=> acetyl-CoA + formate", "Alcohol utilization",
                     "Acetate to acetyl-CoA", "Methane production",
                     "Methane oxidation - Partculate methane monooxygenase",
                     "Methane oxidation - Soluble methane monoxygenase", 
                     "Ammonia oxidation", "Nitrate reduction", "Nitric oxide reduction",
                     "Nitrite oxidation", "Nitrite reduction to ammonia", "Nitrite reduction",
                     "Nitrous oxide reduction", "N2 fixation")
other_metapathway <- c("Sulfate reduction", "Sulfide oxidation", "Sulfite reduction",
                       "Sulfur oxidation", "Sulfur reduction", "Thiosulfate disproportionation", 
                       "Thiosulfate oxidation", "Iron oxidation", "Iron reduction",
                       "Arsenite oxidation", "Arsenate reduction", "Selenate reduction", 
                       "FeFe hydrogenase", "Ni-Fe Hydrogenase")
dat1 <- metabolic_tab %>% 
  filter(Category %in% target_cat) %>%
  mutate(across(starts_with("s_"), ~ifelse(. == "Present", 1, 0))) %>%
  group_by(Category, Function) %>% 
  summarise(across(everything(), sum)) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>% 
  mutate(across(starts_with("s_"), ~ifelse(. >= 1, 1, 0))) %>%
  pivot_longer(cols = -c(Category, Function), names_to = "ID", values_to = 'presence_or_absent') %>%
  mutate(ID = gsub(".Hmm.presence", "", ID))

# merge the tax table and the metabolic table
dat2 <- dat1 %>% filter(ID %in% net_hubs) %>%
  as.data.frame() %>% 
  select(Category, Function, ID, presence_or_absent) %>%
  group_by(Category, Function, ID) %>%
  dplyr::summarize(presence = sum(presence_or_absent)) %>%
  ungroup() # "Adding missing grouping variables" message 


# plot heatmap with no clusters
pacman::p_load(pheatmap)
category <- unique(dat2$Category)

my_col <- c("#E31A1C", "#893f45", "#1F78B4", "#6A3D9A", "#b31b1b", "#5d3954", "#008b8b", "#8b008b", "#e75480", "#872657", "#00bfff",  "#c08081", "#355e3b", "#29ab87", "#f56991", "#32cd32", "#035096", "#bb3385", "#009E73")

plot_list = list()
for (i in 1:13) {
  p <- pheatmap(dat2 %>% filter(Category == category[i]) %>%
                  dplyr::select(Function, ID, presence) %>% 
                  pivot_wider(names_from = Function, values_from = presence) %>%
                  arrange(match(ID, net_hubs)) %>%
                  column_to_rownames(var = "ID"),
                cluster_rows = F,
                cluster_cols = F,
                # show_rownames = FALSE,
                show_rownames = FALSE,
                border_color = 'Black',
                cellwidth = 7, cellheight = 7,
                fontsize_col = 8, fontsize_row = 8,
                color = c(colorRampPalette(c("white", my_col[i]))(2)),
                breaks = c(0, 0.5, 1),
                legend = F, silent = T
  )
  plot_list[[i]] = p
}

p14 <- pheatmap(dat2 %>% filter(Category == category[14]) %>%
                  dplyr::select(Function, ID, presence) %>% 
                  pivot_wider(names_from = Function, values_from = presence) %>%
                  arrange(match(ID, net_hubs)) %>%
                  column_to_rownames(var = "ID"),
                cluster_rows = F,
                cluster_cols = F,
                border_color = 'Black',
                cellwidth = 7, cellheight = 7,
                fontsize_col = 8, fontsize_row = 8,
                color = c(colorRampPalette(c("white", my_col[14]))(2)),
                breaks = c(0, 0.5, 1),
                legend = F, silent = T
)
core_genomes_heatmaps <- cowplot::plot_grid(plot_list[[1]][[4]], plot_list[[2]][[4]], plot_list[[3]][[4]], plot_list[[4]][[4]], 
                               plot_list[[5]][[4]], plot_list[[6]][[4]], plot_list[[7]][[4]], plot_list[[8]][[4]],
                               plot_list[[9]][[4]], plot_list[[10]][[4]], plot_list[[11]][[4]], plot_list[[12]][[4]],
                               plot_list[[13]][[4]], p14[[4]], nrow = 1, align = "h",
                               rel_widths = c(8, 3, 5, 3, 6, 2, 1, 5, 2, 3, 1, 3, 8, 13))

# ggsave(file.path(save.dir, "./figs/MAGs/core_species_heatmap_tax_fun_cor.pdf"),
#        core_genomes_heatmaps, width = 8, height = 8.5, units = "in")
core_genomes_heatmaps



###########
abun_uncollapsed <- abun[, grep("_C", colnames(abun), value = T)]
abun_uncollapsed <- abun_uncollapsed[rowSums(abun_uncollapsed) > 0, , drop = FALSE]
nrow(abun_uncollapsed)
rownames(abun_uncollapsed)[rownames(abun_uncollapsed) %in% net_hubs]

abun_collapsed <- abun[, grep("_T", colnames(abun), value = T)]
abun_collapsed <- abun_collapsed[rowSums(abun_collapsed) > 0, , drop = FALSE]
nrow(abun_collapsed)
rownames(abun_collapsed)[rownames(abun_collapsed) %in% net_hubs]



library(rstatix)
mags.att %>% kruskal_test(Genome_size ~ `Network roles`)

genome_size_uncolla <- mags.att[rownames(mags.att) %in% rownames(abun_uncollapsed), "Genome_size"]
genome_size_colla <- mags.att[rownames(mags.att) %in% rownames(abun_uncollapsed), "Genome_size"]

genome_size_df <- rbind(data.frame(Group = rep("Uncollapsed", length(genome_size_uncolla)), 
                                   genome_size = genome_size_uncolla),
                        data.frame(Group = rep("Collapsed", length(genome_size_colla)), 
                                   genome_size = genome_size_colla))

my_comparisons_group <- list(c('Uncollapsed', 'Collapsed'))
genome_size_group <- genome_size_df %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  ggplot(aes(Group, genome_size/1000000, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.2), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.3, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = Group), shape = 21, size = 1.5, width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_group, 
                     method = "wilcox.test", p.adjust.method = "BH", 
                     bracket.size = 0.3,
                     size = 3.5, tip.length = 0.00) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00", "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78")) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 14),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        strip.text = element_text(colour = "black", size = 12),
        legend.key.size = unit(1,"line"),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "none")

genome_size_group





genome_size <- mags.att %>% select(`Network roles`, Genome_size) %>%
  mutate(`Network roles` = factor(`Network roles`, 
                                  levels = c("Network hubs", "Module hubs",
                                             "Connector hubs", "Peripheral species"))) %>%
  ggplot(aes(`Network roles`, Genome_size, fill = `Network roles`)) +
  geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.8) +
  geom_jitter(aes(fill = `Network roles`), shape = 21, size = 1.5, width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.1, npcy = 0.95, label = "Kruskal-Wallis test: P < 0.001") +
  # ggpmisc::annotate(aes(npcx = 0.1, npcy = 0.95, label = "Kruskal-Wallis test: P < 0.001"), inherit.aes = F) +
  labs(x = "Network roles", y = "Genome size") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")

genome_size


# Fitting exponential decay curve
library(stats)
model <- nls(`Node degree` ~ a * exp(b * Genome_size) + c, data = mags.att, 
             start = list(a = 200, b = -0.25, c = 6))
summary(model)

lm(Genome_size ~ predict(model), data = mags.att) %>% broom::glance()

#plot
genome_size_degree_plot <- mags.att %>%
  mutate(`Network roles` = factor(`Network roles`, 
                                  levels = c("Network hubs", "Module hubs",
                                             "Connector hubs", "Peripheral species"))) %>%
  ggplot(aes(Genome_size, `Node degree`, colour = `Network roles`)) + 
  geom_point(size = 3.5, alpha = 0.8) +
  geom_smooth(method = "nls", 
              formula = "y ~ A * exp(B * x) + C", 
              method.args = list(start = c(A = 34.4961, B = -0.4447, C = 5.5699)),
              size = 1, se = F, linetype = 1,  colour = 'black') +
  scale_color_manual(values = c("#e95f5c", "#79ceb8", "#5cc3e8", "#ffdb00")) +
  xlab("Genome size (Mb)") + ylab("Node degree") +
  ggpp::annotate("text_npc", npcx = 0.8, npcy = 0.95, label = "R^2 = 0.77; P < 0.001") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(color = "black", size = 10),
        axis.ticks.length = unit(0.2, "lines"), 
        axis.ticks = element_line(color = "black"), 
        axis.line = element_blank(), 
        axis.text.y = element_text(colour = "black", size = 8), 
        axis.text.x = element_text(colour = "black", size = 8), 
        strip.text = element_text(size = 14), legend.position = "none")
# ggsave(file.path(save.dir.revision, "FTC_curve.pdf"),
#        FTC_plot, width = 4, height = 3.5, units = "in")
genome_size_degree_plot


library(ggpubr) # 'ggplot2' Based Publication Ready Plots
library(ggpmisc) # Miscellaneous Extensions to 'ggplot2'
mags.att %>% 
  ggplot(aes(Genome_size, `Node degree`, fill = `Network roles`))+
  #point
  geom_point(size = 2)+
  #regression line
  geom_smooth(method = "loess", color = "black", fill = "grey", se = T, 
              formula = y ~ x,
              linetype = 1, alpha = 0.2)+
  # #add p values
  # stat_poly_eq(formula = y ~ x, 
  #              aes(label = paste(after_stat(adj.rr.label),..p.value.label..,sep = "~~~~")), parse = TRUE,
  #              size= 2.5, label.x = "left")+
  # #color
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  #legend
  guides(color=guide_legend(override.aes = list(size=5,alpha=1)),
         size = "none")+
  #axis labes
  labs(x = "Genome size", y = "Node degree")+
  main_theme



###########
# Flux analysis
wd = 'e:/thermokarst_gully/data'
flux.file <- "flux_data.txt"
setwd(wd)
flux <- read.table(flux.file, header = TRUE, sep = "\t",
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE)

flux_ave.time_df <- flux %>% select(c("Flux_type", "Gully_id", "Rep", "Uncollapsed", "Collapsed")) %>%
  group_by(Flux_type, Gully_id, Rep) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
flux_ave.site_df <- flux %>% select(c("Flux_type", "Gully_id", "Uncollapsed", "Collapsed")) %>%
  group_by(Flux_type, Gully_id) %>%
  get_summary_stats(c("Uncollapsed", "Collapsed"), type = "common") #or using type = "mean_sd"

core_mags_abun <- abun[net_hubs, ] %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Sample_id") %>%
  mutate(Gully_id = case_when(grepl("G1", Sample_id) ~ "EB",
                              grepl("G2", Sample_id) ~ "ML",
                              grepl("G3", Sample_id) ~ "RS",
                              grepl("G4", Sample_id) ~ "SLH",
                              grepl("G5", Sample_id) ~ "HSX",
                              grepl("G6", Sample_id) ~ "HH")) %>%
  mutate(Group = case_when(grepl("_C", Sample_id) ~ "Uncollapsed",
                           grepl("_T", Sample_id) ~ "Collapsed")) %>%
  mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
  select(Gully_id, Group, net_hubs, Total) %>%
  group_by(Gully_id, Group) %>%
  get_summary_stats(c(net_hubs, Total), type = "common") #or using type = "mean_sd"
  
# write.csv(core_mags_abun, "./core_mags_abun_total.csv")
# write.csv(flux_ave.site_df, "./flux_ave.site_df.csv")


mags_flux.file <- "core_mags_abun_total.txt"
mags_flux <- read.table(mags_flux.file, header = TRUE, sep = "\t", row.names = 1, 
           as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
           check.names = FALSE) %>% dplyr::select(Gully_id, Group, variable, mean, se) %>%
  pivot_wider(names_from = variable, values_from = c(mean, se))

# write.csv(mags_flux, "./mags_flux_plot.csv")


library(tidyverse)

# Load the data
wd = 'e:/thermokarst_gully/data'
flux_plot_file <- "mags_flux_plot_longdata.txt"
setwd(wd)
flux_plot_df <- read.table(flux_plot_file, header = TRUE, sep = "\t", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE)

#plot
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(ggpubr) # ggplot2 Based Publication Ready Plots
library(ggpmisc) # Miscellaneous Extensions to 'ggplot2'
colnames(flux_plot_df) <- make.names(colnames(flux_plot_df), unique = TRUE)
flux_plot_df %>% filter(Flux_type %in% c("CH4", "CO2")) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  ggplot(aes(Mean_abundance, Mean_flux))+
  #errorbar for x variable
  geom_errorbar(aes(xmin = Mean_abundance - Se_abundance, 
                    xmax = Mean_abundance + Se_abundance), linewidth = 0.4, width = 0)+
  #errorbar for y variable
  geom_errorbar(aes(ymin = Mean_flux - Se_flux,
                    ymax = Mean_flux + Se_flux), linewidth = 0.4, width = 0)+
  #point
  geom_point(aes(color = Group), size = 2)+
  #regression line
  geom_smooth(method = "lm", color = "black", fill = "grey", se = T, 
              formula = y ~ x,
              linetype = 1, alpha = 0.2)+
  #add p values
  # stat_cor(method = "pearson",label.x = 17, label.y = 27, size=4.5)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(adj.rr.label),..p.value.label..,sep = "~~~~")), parse = TRUE,
               size= 2.5, label.x = "left")+
  # #color
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  #legend
  guides(color=guide_legend(override.aes = list(size=5,alpha=1)),
         size = "none")+
  #axis labes
  labs(x = "Mean relative abundance (%)", y = NULL, color = NULL)+
  facet_grid(Flux_type ~ Genome_id, scales = "free") +
  main_theme


##################Max Growth Rate##################
devtools::install_github("jlw-ecoevo/gRodon2")
library(gRodon)
library(Biostrings)

path_to_metagenome <- system.file('extdata',
                                  'ERR2143764_fastp_prokka_scaffolds.ffn.gz',
                                  package = 'gRodon')
genes <- readDNAStringSet(path_to_metagenome)
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
predictGrowth(genes, highly_expressed, mode = "metagenome_v1")



path_to_coverage <- system.file('extdata',
                                'ERR2143764_fastp_map2ffn_counts.tsv',
                                package = 'gRodon')
read_depths <- read.delim(path_to_coverage,
                          stringsAsFactors = FALSE)
depths <- read_depths$meandepth
#Make sure in the correct order
names(depths) <- read_depths$X.rname
depth_of_coverage <- depths[gsub(" .*", "", names(genes))]
head(depth_of_coverage)

predictGrowth(genes, 
              highly_expressed, 
              mode = "metagenome_v2", 
              depth_of_coverage = depth_of_coverage)



####################################
####################################
####################################
library(gRodon)
library(Biostrings)

# Specify the directory containing the .ffn files
input_directory <- "E:/thermokarst_gully/data/metagenome/MAGs/prokka_ffn/"
output_file <- "E:/thermokarst_gully/data/metagenome/MAGs/growth_rate_results.csv"

# Get a list of all .ffn files in the directory
fnn_files <- list.files(input_directory, pattern = "\\.ffn$", full.names = TRUE)

# Initialize an empty data frame to store the results
results <- data.frame(
  File = character(),
  CUBHE = numeric(),
  GC = numeric(),
  GCdiv = numeric(),
  ConsistencyHE = numeric(),
  CUB = numeric(),
  CPB = numeric(),
  FilteredSequences = numeric(),
  nHE = numeric(),
  dCUB = numeric(),
  d = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each .ffn file and calculate the metrics
for (file in ffn_files) {
  cat("Processing:", file, "\n")
  tryCatch({
    # Read the DNA sequences
    genes <- readDNAStringSet(file)
    
    # Identify highly expressed genes (ribosomal proteins)
    highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)
    
    # Predict the growth metrics
    growth_rate <- predictGrowth(genes, highly_expressed, mode = "partial")
    
    # Append the metrics to the results data frame
    results <- rbind(results, data.frame(
      File = basename(file),
      CUBHE = growth_rate$CUBHE,
      GC = growth_rate$GC,
      GCdiv = growth_rate$GCdiv,
      ConsistencyHE = growth_rate$ConsistencyHE,
      CUB = growth_rate$CUB,
      CPB = growth_rate$CPB,
      FilteredSequences = growth_rate$FilteredSequences,
      nHE = growth_rate$nHE,
      dCUB = growth_rate$dCUB,
      d = growth_rate$d,
      LowerCI = growth_rate$LowerCI,
      UpperCI = growth_rate$UpperCI
    ))
  }, error = function(e) {
    cat("Error processing", file, ":", e$message, "\n")
  })
}

# Write the results to a CSV file
write.csv(results, file = output_file, row.names = FALSE)

cat("Processing complete. Results saved to:", output_file, "\n")



##########################################
############################################
# Set the genomes_dir directory
genomes_dir <- "/data01/kangluyao/thermokarst_gully/binning/prokka_fna"
genomes_files = list.files(genomes_dir, full.names = T, recursive = T, pattern = ".fna$")
out_dir = "/data01/kangluyao/thermokarst_gully/binning/microtraits"

#Determine the number of "cores"
message("Number of cores:", parallel::detectCores(), "\n")


#Here, we use ~70% of the available cores to process 100 genomes:
library("tictoc")
tictoc::tic.clearlog()

microtrait_results = extract.traits.parallel(genomes_files, out_dir = dirname(genomes_files), ncores = 30)

#Pull the paths to the corresponding "rds" files as follows:
rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file", mc.cores = 30))

#Building trait matrices from microTrait outputs
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = 30)
                                             






library(parallel)

library(tictoc)

library(microtrait)



message("Running on: ", system("cat /proc/cpuinfo | grep "model name" | uniq | sed 's/.*: //'", intern = TRUE), "\n")

message("Number of cores:", detectCores(), "\n")

tictoc::tic.clearlog()
tictoc::tic(paste0("Running microtrait for ", length(genomes_files)))

microtrait_results <- mclapply(1:length(genomes_files), function(i) {
  r = extract.traits(genomes_files[i], out_dir = "/data01/kangluyao/thermokarst_gully/binning/microtraits")
  saveRDS(r, file = file.path("/data01/kangluyao/thermokarst_gully/binning/microtraits/rds", 
                              paste0(fs::path_file(genomes_files[i]), ".microtrait.rds.1")))
  r
}, mc.cores = 30)

tictoc::toc(log = "TRUE")

rds_files = unlist(parallel::mclapply(genomeset_results, "[[", "rds_file", mc.cores = 30))
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = 1)
