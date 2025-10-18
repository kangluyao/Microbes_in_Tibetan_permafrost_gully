rm(list = ls())
# Set work directory
setwd('e:/thermokarst_gully/')
wd_16s <- file.path(getwd(),"data/16S/rdp")
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")

# Loading packages
pacman::p_load(vegan, tidyverse, ggrepel, EnhancedVolcano, edgeR, ggplot2)

# Data input
source("script/read_data.R")
# Permutational multivariate analysis of variance using distance matrices (adonis) to test the difference of functional composition
library(vegan)
# determine the dissimilarity matrix based on the bray-curties distance
fun_dist <-vegdist(t(ko_tpm_table), "bray" )
# permanova, ANOSIM and MRPP analysis
adonis2(fun_dist ~ Group, data = metadata)
mrpp(fun_dist, metadata$Group, perm = 999)
anosim(fun_dist, metadata$Group, perm = 999)

# PCoA plot
ord.fun <-  cmdscale(fun_dist,  k = 2, eig = T, add = T)
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 6),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.position = "none")
pcoa_fun_plot <- data.frame(Group = metadata$Group, scores(ord.fun)) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
  ggplot(aes(x = Dim1, y = Dim2)) + 
  geom_point(size = 1, alpha = 0.8, shape = 21, colour = "black", aes(fill = Group)) + 
  stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1, 
               show.legend = FALSE, level = 0.95) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  # scale_x_continuous(expand = c(0.03, 0.03)) +
  # scale_y_continuous(expand = c(0.03, 0.03)) +
  labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), digits = 3), "%)", sep = ""),
       y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), digits = 3), "%)", sep = "")) +
  main_theme +
  theme(legend.background = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.key.size = unit(0.4, 'cm'))
# creat a directory for kegg-gene results
if (!dir.exists(file.path(save.dir, "figs/kegg/"))) {
  dir.create(file.path(save.dir, "figs/kegg/"))
}
# ggsave(file.path(save.dir, "./figs/kegg/PCoA_fun_bray.pdf"),
#        pcoa_fun_plot, width = 3.5, height = 2.5, units = "in")
pcoa_fun_plot



# Determine the funtional homogenization
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
          'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
distance_fun_data <- lapply(vars, function(x) 
  usedist::dist_subset(fun_dist, 
                       grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
distance_fun_df <- distance_fun_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
library(lme4)
library(lmerTest)
lmm_dis_fun_mod <- lmer(distance ~ Group + Time + Slope + MAP + (1 | Gully_id),  data = distance_fun_df)
summary.model <- function(model){
  F.value <- anova(model)$'F value'
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
summary.model(lmm_dis_fun_mod)


# Create a box plot
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
library(gghalves)
dis_fun_plot <- distance_fun_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.5, npcy = 0.95, 
                 size = 2, label = "P < 0.001") +
  labs(x = "Group", y = "Dissimilarity in gene composition") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")
dis_fun_plot


# Determine the funtional homogenization
#permanova test the difference in compositional variance
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
similar_fun_data <- lapply(vars, function(x) usedist::dist_subset(fun_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

library(lme4)
library(lmerTest)
# determine the effect size of the permafrost thawing for the functional gene composition
dist_scale <- similar_fun_data %>% 
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

fm1 <- lmer(distance ~ Group + (1 | Gully_id), data = dist_scale)
presult <- car::Anova(fm1, type = 2)
coefs <- coef(summary(fm1))[, "Estimate"]  ##four coefs
names(coefs) <- paste0(names(coefs), ".mean")
SEvalues <- coef(summary(fm1))[, "Std. Error"]  ##standard errors
names(SEvalues) <- paste0(names(SEvalues), ".se")
tvalues <- coef(summary(fm1))[, "t value"]  ##t values
names(tvalues) <- paste0(names(tvalues), ".t")
chisqP <- c(presult[, 1], presult[, 3])
names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
result <- data.frame(variables = "Function", GroupCollapsed.mean = coefs[2], GroupCollapsed.se = SEvalues[2], sig = "***")


fun_dist_comparison <-  ggplot(data = result, aes(x = variables, y = GroupCollapsed.mean, colour = "#e95f5c")) +
  geom_hline(aes(yintercept = 0), size = 0.7,  colour = "gray2")+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom") +
  theme(legend.position = "none", 
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.1), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_rect(colour = "white"))

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_div_comparison, width = 2.7, height = 5, units = "in")
# fun_dist_comparison <- ggplotGrob(fun_dist_comparison)
# pcoa_fun_plot <- pcoa_fun_p1 + 
#   annotation_custom(fun_dist_comparison, 
#                     xmin = 0.06, xmax = 0.14, 
#                     ymin = -0.075, ymax = -0.03)
fun_dist_comparison

## Explore the losted genes after permafrost thawing
# library
library(ggvenn)
x <- list(
  Control = ko_tpm_table %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_C', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames(), 
  Collapsed = ko_tpm_table %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_T', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames()
)
# plot
fun.venn <- ggvenn(
  x, 
  fill_color = c("#f8766d", "#a3a500", "#00b0f6"),
  stroke_color = NA,
  set_name_size = 4,
  text_size = 4,
  show_percentage = F
)
fun.venn

#####
traits_cate_file <- 'E:/thermokarst_gully/data/metagenome/MAGs/microtraits/microtraits_cate.txt'
traits_cate_df <- read.table(traits_cate_file, header = TRUE, sep = "\t", 
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE)

ko_tpm_table %>% rownames_to_column("KO") %>%
  inner_join(traits_cate_df[, c(1,2,5)], by = "KO") %>%
  relocate(c(KO, `microtrait_hmm-name`, `microtrait_hmm-description`), .before = 1)


# loading packages
library(tidyverse)
#reading data
KO_L2_file <- 'E:/thermokarst_gully/data/metagenome/KEGG.PathwayL2.raw.txt'

KO_L2_df <- read.table(KO_L2_file, header = TRUE, sep = "\t", 
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE)

# Using dplyr to process the data
aggre_data <- KO_L2_df %>%
  # Add a column for row means (excluding the L2 column)
  mutate(RowMean = rowMeans(dplyr::select(., -PathwayL2))) %>%
  # Arrange rows by RowMean in descending order
  arrange(desc(RowMean)) %>%
  # Add a new column to classify rows as "Top 10" or "Other"
  mutate(Pathway = if_else(dplyr::row_number(.) <= 25, as.character(PathwayL2), "Other")) %>%
  # Group by "Pathway" and summarize
  group_by(Pathway) %>%
  dplyr::summarize(., across(starts_with("G"), sum)) %>% # Sum numerical columns 
  # Retain row mean for ordering purposes
  mutate(RowMean = rowMeans(select(., -Pathway))) %>%
  # Arrange by RowMean in descending order, placing "Other" last
  mutate(PathwayOrder = if_else(Pathway == "Other", -Inf, RowMean)) %>%
  arrange(desc(PathwayOrder)) %>%
  # factor(Pathway = factor(Pathway, levels = Pathway)) %>% # Order the Pathway
  select(-PathwayOrder, -RowMean)  # Drop intermediate columns

# Pull out the pathway for ordering purposes
pathway_order <- aggre_data %>% pull(Pathway)

##Stack plot using ggplot2 
pathway_plot <- aggre_data %>% 
  pivot_longer(-Pathway, names_to = 'Sample_id', values_to = 'TPM') %>%
  mutate(Pathway = factor(Pathway, ordered = T, levels = pathway_order),
         Group = case_when(grepl("_C", Sample_id) ~ "Un-collapsed",
                           grepl("_T", Sample_id) ~ "Collapsed"),
         Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(sample_id = rep(1:60, 26)) %>%
  group_by(Group, Sample_id) %>%
  mutate(prop = (TPM * 100 / sum(TPM))) %>%
  ggplot(aes(x = sample_id, y = prop, fill = Pathway)) +
  geom_area() +
  labs(x = 'Sample', y = 'Proportion (%)') +
  scale_fill_manual(values = pals::brewer.accent(26)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~ Group, scales = 'free_x', space = 'free_x') +
  main_theme +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, "lines"),
        panel.spacing = unit(0, "lines"))
pathway_plot


# Combine all plot
library(cowplot)
path_plot <- ggdraw() +
  draw_plot(fun.venn, x = 0, y = 0.5, width = 1/3, height = 1/2) +
  draw_plot(pcoa_fun_p1, x = 0, y = 0, width = 1/3, height = 1/2) +
  draw_plot(pathway_plot, x = 1/3, y = 0, width = 2/3, height = 1) +
  draw_plot_label(label = c("a", "b", 'c'), size = 8,
                  x = c(0, 0, 1/3), y = c(1, 0.5, 1))
# ggsave(file.path(save.dir, "/figs/kegg/path_plot1.pdf"),
#        path_plot, width = 5.5, height = 4, units = "in")
path_plot

# Test the effect size of permafrost collapse on the functional pathways
pathway_scale <- KO_L2_df %>% filter(PathwayL2 != "") %>% 
  column_to_rownames("PathwayL2") %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample_id") %>% 
  mutate(Gully_id = sapply(stringr::str_split(Sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(Group = case_when(grepl("_C", Sample_id) ~ "Uncollapsed",
                           grepl("_T", Sample_id) ~ "Collapsed")) %>%
  mutate(across(where(is.numeric), scale)) %>%
  relocate(where(is.character)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(-Sample_id)

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
pathway_S1 <- sapply(3:ncol(pathway_scale), function(j) {
  if (length(unique(pathway_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(pathway_scale[, j] ~ Group + (1 | Gully_id), data = pathway_scale)
    
    presult <- car::Anova(fm1, type = 2)
    coefs <- coef(summary(fm1))[, "Estimate"]  ##four coefs
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm1))[, "Std. Error"]  ##standard errors
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    tvalues <- coef(summary(fm1))[, "t value"]  ##t values
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    chisqP <- c(presult[, 1], presult[, 3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    result <- c(coefs, tvalues, SEvalues, chisqP)
  }
  result
})
colnames(pathway_S1)<-colnames(pathway_scale)[-c(1:2)]
pathway_S1[1:6, 1:6]

ko_9_levels_file <- 'E:/thermokarst_gully/data/metagenome/ko_9_levels.txt'

ko_9_levels <- read.table(ko_9_levels_file, header = TRUE, sep = "\t", 
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE)

pathway_index <- unique(ko_9_levels[, c("L1", "L2")]) %>% pull(L2)
pathway_index <- pathway_index[pathway_index != ""]

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*",  "")))}
diff_in_pathway_plot <- pathway_S1 %>%
  t() %>%
  data.frame(check.names = F) %>%
  tibble::rownames_to_column(., "variables") %>%
  filter(variables %in% pathway_index) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(pathway_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, color = colour)) +
  geom_hline(aes(yintercept =0), size=0.7,  colour="gray2")+
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se + 0.2)),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") + 
  scale_colour_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0, xmax = 7.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#ff6666") +
  annotate("rect", xmin = 7.5, xmax = 15.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#f19837") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#e56eee") +
  annotate("rect", xmin = 20.5, xmax = 23.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#5fb236") +
  annotate("rect", xmin = 23.5, xmax = 28.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#035096") +
  annotate("rect", xmin = 28.5, xmax = 38, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#8b008b") +
  main_theme +
  theme(strip.background = element_rect(fill = c("#FFF6E1")), legend.position = "none")

# ggsave(file.path(save.dir, "/figs/kegg/diff_in_pathway_plot.pdf"),
#        diff_in_pathway_plot, width = 3, height =6, units = "in")
diff_in_pathway_plot

## Select the genes involved in C, N, S, and Fe cyclings
#### Arrange the log fold change table
# loading packages
pacman::p_load(phyloseq, picante, microbiome, readxl, tidyverse, ggpubr, ggplot2)

# Read data
sel_ko <- read_csv(file.path(wd_fun, 'CN_ko_input.csv'), col_names = T)

C_N_ko_count_tab <- ko_count_table[rownames(ko_count_table) %in% sel_ko$KO, ]
C_N_ko_count_tab <- cbind(KO = rownames(C_N_ko_count_tab), C_N_ko_count_tab)

C_N_ko_tpm_tab <- ko_tpm_table[rownames(ko_tpm_table) %in% sel_ko$KO, ]
C_N_ko_tpm_tab <- cbind(KO = rownames(C_N_ko_tpm_tab), C_N_ko_tpm_tab)


## prepare the data for LMMs
C_N_path1 <- C_N_ko_tpm_tab %>% 
  inner_join(sel_ko, c("KO" = "KO")) %>% 
  select(-c(1)) %>%
  group_by(Category, Subcategory, Enzyme_protein_encoded) %>%
  summarise(across(everything(), sum)) %>% data.frame() %>%
  mutate(Subcategory = factor(Subcategory, levels = unique(sel_ko$Subcategory), ordered = T))
C_N_path2 <- C_N_path1 %>% select(c(3:63)) %>% 
  tibble::column_to_rownames('Enzyme_protein_encoded')

# Test the difference in functional genes involving in carbon, nitrogen, sulfur, and other elements using linear mixed models
count_scale <- data.frame(metadata[, c("Gully_id", "Group", "Time", "Slope", "MAP")], t(C_N_path2), check.names=FALSE) %>% 
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
count_S1 <- sapply(6:ncol(count_scale), function(j) {
  if (length(unique(count_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(count_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = count_scale)
    
    presult <- car::Anova(fm1, type = 2)
    coefs <- coef(summary(fm1))[, "Estimate"]  ##four coefs
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm1))[, "Std. Error"]  ##standard errors
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    tvalues <- coef(summary(fm1))[, "t value"]  ##t values
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    chisqP <- c(presult[, 1], presult[, 3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    result <- c(coefs, tvalues, SEvalues, chisqP)
  }
  result
})
colnames(count_S1)<-colnames(count_scale)[-c(1:5)]
count_S1[1:6, 1:6]

# prepare the names for each category
C_names <- c("Alpha-amylase", "Glucoamylase", "Pullulanase", "Isopullulanase",
             "Arabinofuranosidase", "Beta_mannanase", "Xylanase", "Xylose isomerase",
             "Beta-glucosidase",  "Endoglucanase", "Exoglucanase",
             "Acetylglucosaminidase", "Endochitinase", "Exochitinase", "Pectinase", 
             "Aryl-aldehyde oxidase", "Isocitrate lyase", "Limonene-1, 2-epoxide hydrolase", 
             "Malate synthase", "Vanillate demethylase", 
             "Glyoxal oxidase", "Phenol oxidase (tyrosinase)",
             "Methyl coenzyme M reductase", "Particulate methane monooxygenase", "Soluable methane monooxygenase",
             "Pyruvate oxidation", "Pyruvate <=> acetyl-CoA + formate", "Acetogenesis", 
             "Lactate utilization", "Acetate to acetyl-CoA", "RuBisCo")
N_names <- c("amoA",	"amoB",	"amoC",	"hao",	"nxrA",	"nxrB", "narG",	"narH",	"narI", "napA",	
             "napB", "nirK", "nirS",	"norB",	"norC",	"nosZ", "HzsC", "HzsB", "HzsA", "hdh",
             "nrfA",	"nrfH", "nirB",	"nirD",	
             "nasA",	"nasB", "narB",	"NR", "NIT-6", "nirA",	"nifD", "nifK", "nifH", "nrtA",	"nrtB", "nrtC",	
             "nrtD", "nmo", "gdh_K00261", "gdh_K00262", "gdh_K15371", "glsA", "ureA", "ureC", "glnA")
S_names <- c("sat", "cysC", "cysD", "cysNC", "csyH", "cysJ", "cysN", "aprA", "aprB", "dsrA", "dsrB",
             "asrA", "asrB", "asrC", "sir", "sor", "sreA", "sreB", "sreC", "hydA", "hydD", "hydB",
             "hydG", "psrA", "psrB", "psrC", "sqr", "fccA", "fccB", "soxD", "soxX", "soxA", "soxB",
             "soxC", "soxY", "soxZ", "ttrA", "ttrB", "ttrC", "phsA", "phsB", "phsC")
Other_names <- c("Cyc1", "Cyc2", "MtrA", "MtrB", "MtrC", "arsC (grx)", "arsC (trx)", 
                 "aioA", "arsM", "ygfM", "xdhD", "YgfK")

# Effect size plot for carbon metabolisms
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*",  "")))}
eff.size_carbon_plot <- count_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  filter(variables %in% C_names) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = C_names)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, fill = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2") +
  geom_bar(stat="identity", width = 0.6) +
  # geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), 
                cex = 0.5) +
  geom_text(aes(label = sig, x = variables, 
                y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_fill_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  scale_x_discrete(position = "bottom") +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#E31A1C") +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#893f45") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#1F78B4") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#6A3D9A") +
  annotate("rect", xmin = 20.5, xmax = 22.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 22.5, xmax = 25.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#00bfff") +
  annotate("rect", xmin = 25.5, xmax = 30.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 30.5, xmax = 31.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#8b008b") +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )
my_col <- c("#E31A1C", "#893f45", "#1F78B4", "#6A3D9A", "#b31b1b", "#5d3954", "#008b8b", "#8b008b", "#e75480", "#872657", "#00bfff",  "#c08081", "#355e3b", "#29ab87", "#f56991", "#32cd32", "#035096", "#bb3385", "#009E73")

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir, "/figs/kegg/eff.size_carbon_plot.pdf"),
#        eff.size_carbon_plot, width = 5.5, height = 2.5, units = "in")
eff.size_carbon_plot

# Effect size plot for nitrogen metabolisms
eff.size_nitrogen_plot <- count_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  filter(variables %in% N_names) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = N_names)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, 
                y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#9457eb") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 20.5, xmax = 22.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 22.5, xmax = 25.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 25.5, xmax = 27.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  annotate("rect", xmin = 27.5, xmax = 28.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#9457eb") +
  annotate("rect", xmin = 28.5, xmax = 29.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 29.5, xmax = 30.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )
my_col <- c("#E31A1C", "#893f45", "#1F78B4", "#6A3D9A", "#b31b1b", "#5d3954", "#008b8b", "#8b008b", "#e75480", "#872657", "#00bfff",  "#c08081", "#355e3b", "#29ab87", "#f56991", "#32cd32", "#035096", "#bb3385", "#009E73")

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir, "/figs/env/diff_in_env_plot.pdf"),
#        diff_in_env_plot, width = 2.5, height = 3.5, units = "in")
eff.size_nitrogen_plot

# Effect size plot for sulfur metabolisms
eff.size_sulfur_plot <- count_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  filter(variables %in% S_names) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = S_names)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, 
                y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#9457eb") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 20.5, xmax = 22.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 22.5, xmax = 25.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 25.5, xmax = 27.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  annotate("rect", xmin = 27.5, xmax = 28.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#9457eb") +
  annotate("rect", xmin = 28.5, xmax = 29.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 29.5, xmax = 30.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )
my_col <- c("#E31A1C", "#893f45", "#1F78B4", "#6A3D9A", "#b31b1b", "#5d3954", "#008b8b", "#8b008b", "#e75480", "#872657", "#00bfff",  "#c08081", "#355e3b", "#29ab87", "#f56991", "#32cd32", "#035096", "#bb3385", "#009E73")

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir, "/figs/env/diff_in_env_plot.pdf"),
#        diff_in_env_plot, width = 2.5, height = 3.5, units = "in")
eff.size_sulfur_plot

# Effect size plot for other elements metabolisms
eff.size_other_plot <- count_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  filter(variables %in% Other_names) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = Other_names)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, 
                y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 4.5, xmax = 8.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )
my_col <- c("#E31A1C", "#893f45", "#1F78B4", "#6A3D9A", "#b31b1b", "#5d3954", "#008b8b", "#8b008b", "#e75480", "#872657", "#00bfff",  "#c08081", "#355e3b", "#29ab87", "#f56991", "#32cd32", "#035096", "#bb3385", "#009E73")

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir, "/figs/env/diff_in_env_plot.pdf"),
#        diff_in_env_plot, width = 2.5, height = 3.5, units = "in")
eff.size_other_plot




# Reading the microTraits catergory table
library(data.table)
trait_cater <- fread("E:/thermokarst_gully/data/metagenome/MAGs/microtraits/microtraits_cate.txt", 
              header = TRUE, stringsAsFactors = FALSE)
setnames(trait_cater, 
         old = c("microtrait_trait.name1", "YAS", "YAS.2", "YAS.3", "microtrait_hmm.dbxref_kegg"), 
         new = c("microtrait_trait", "level1", "level2", "level3", "KO"))

# Extract a subset of the KO abundance table based on the KO IDs corresponding to specific traits
aggre_trait_data <- ko_tpm_table %>%
  rownames_to_column("KO") %>%
  inner_join(trait_cater, by = "KO") %>%
  filter(!is.na(microtrait_trait))

# fwrite(aggre_trait_data, "E:/thermokarst_gully/data/metagenome/MAGs/microtraits/aggre_trait_data.csv")
aggre_trait_group_data <- aggre_trait_data %>%
  column_to_rownames("KO") %>%
  select(level1, level2, level3, 1:60) %>%
  group_by(level1, level2, level3) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(level2)) %>%
  select(3:63) %>%
  column_to_rownames("level3") %>%
  t() %>% data.frame()

# Merge the genes with the metadata
aggre_trait_env_data <- aggre_trait_group_data %>%
  rownames_to_column("Sample_name") %>%
  left_join(metadata[, c("Sample_name", "Gully_id", "Group", 
                         "MAP", "Time", "Slope")],
            by = "Sample_name")
# Linear mixed models test the effect of collapsed and gully_id
gene_id_trait <- c("aromatic.acid.transport",
                   "biopolymer.transport", "carbohydrate.transport", 
                   "carboxylate.transport", "free.amino.acids.transport",
                   "ion.transport", "lipid.transport", "N.compound.transport", 
                   "nucleic.acid.component.transport", 
                   "organophosphorus.transport", "osmolyte.transport",
                   "other.transport", "peptide.transport", 
                   "S.compound.transport", "secondary.metabolite.transport", 
                   "vitamin.transport", "complex.carbohydrate.depolymerization",
                   "simple.compound.degradation", "C1.compounds", "N.compounds",
                   "S.compounds", "P.compounds", "Fermentation",
                   "aerobic.respiration", "anaerobic.respiration",
                   "chemolithoautotrophy", "photosystem", "pigments", 
                   "General", "high.temperature", "low.temperature",
                   "desiccation.osmotic.salt.stress", "pH.stress",
                   "oxidative.stress", "oxygen.limitation", "envelope.stress")

gene_trait_scale <- aggre_trait_env_data %>% 
  cbind(metadata[, c("Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", gene_id_trait))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
gene_trait_S1 <- sapply(6:ncol(gene_trait_scale), function(j) {
  if (length(unique(gene_trait_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(gene_trait_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), 
                data = gene_trait_scale)
    
    presult <- car::Anova(fm1, type = 2)
    coefs <- coef(summary(fm1))[, "Estimate"]  ##four coefs
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm1))[, "Std. Error"]  ##standard errors
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    tvalues <- coef(summary(fm1))[, "t value"]  ##t values
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    chisqP <- c(presult[, 1], presult[, 3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    result <- c(coefs, tvalues, SEvalues, chisqP)
  }
  result
})
colnames(gene_trait_S1)<-colnames(gene_trait_scale)[-c(1:5)]
data.frame(gene_trait_S1)

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", "")))}
# Create a plot
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
all_gene_trait_comparison <- gene_trait_S1 %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>% filter(variables %in% gene_id_trait) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(gene_id_trait))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.1),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2.5, 2.5)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0, xmax = 8.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#e56eee") +
  annotate("rect", xmin = 8.5, xmax = 14.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#ffdb00") +
  annotate("rect", xmin = 14.5, xmax = 18.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#5cc3e8") +
  annotate("rect", xmin = 18.5, xmax = 20.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#e95f5c") +
  annotate("rect", xmin = 20.5, xmax = 36.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#79ceb8") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path("E:/thermokarst_gully/result2/all_tax_div_comparison.pdf"),
#        all_tax_div_comparison, width = 2.5, height = 5, units = "in")
all_gene_trait_comparison


aggre_trait_aqui <- aggre_trait_data %>%
  filter(!is.na(level2) & level1 == "Resource Acquisition") %>%
  select(KO, 2:61) %>%
  column_to_rownames("KO")

aggre_trait_use <- aggre_trait_data %>%
  filter(!is.na(level2) & level1 == "Resource Use") %>%
  select(KO, 2:61) %>%
  column_to_rownames("KO")

aggre_trait_stress <- aggre_trait_data %>%
  filter(!is.na(level2) & level1 == "Stress Tolerance") %>%
  select(KO, 2:61) %>%
  column_to_rownames("KO")


library(vegan)
# Transform the gene abundance data
aggre_trait_aqui_trans <- decostand(t(aggre_trait_aqui), method = "hellinger")
aggre_trait_use_trans <- decostand(t(aggre_trait_use), method = "hellinger")
aggre_trait_stress_trans <- decostand(t(aggre_trait_stress), method = "hellinger")
aggre_trait_aqui_dist <-vegdist(aggre_trait_aqui_trans, "bray")
aggre_trait_use_dist <-vegdist(aggre_trait_use_trans, "bray")
aggre_trait_stress_dist <-vegdist(aggre_trait_stress_trans, "bray")

#permanova test the difference in compositional variance
adonis2(aggre_trait_aqui_dist ~ Group, data = metadata)
adonis2(aggre_trait_use_dist ~ Group, data = metadata)
adonis2(aggre_trait_stress_dist ~ Group, data = metadata)

# PCoA plot
PCoA_plot_fun <- function(dist) {
  ord.fun <-  cmdscale(dist,  k = 2, eig = T, add = T)
  pcoa.plot <- data.frame(Group = metadata$Group, scores(ord.fun)) %>%
    mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
    ggplot(aes(x = Dim1, y = Dim2)) + 
    geom_point(size = 1, alpha = 0.8, shape = 21, colour = "black", aes(fill = Group)) + 
    stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1, 
                 show.legend = FALSE, level = 0.95) +
    scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    # scale_x_continuous(expand = c(0.03, 0.03)) +
    # scale_y_continuous(expand = c(0.03, 0.03)) +
    labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), digits = 3), "%)", sep = "")) +
    main_theme +
    theme(legend.background = element_blank(),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key = element_blank(),
          legend.position = c(0.85, 0.85),
          legend.key.size = unit(0.4, 'cm'))
  return(pcoa.plot)
}
# Set the main theme for ggplot2
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 6),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.position = "none")
# Plots
trait_aqui.fun <- PCoA_plot_fun(aggre_trait_aqui_dist)
trait_use.fun <- PCoA_plot_fun(aggre_trait_use_dist)
trait_stress.fun <- PCoA_plot_fun(aggre_trait_stress_dist)


## Explore the difference in taxonomic variance between uncollapsed and collapsed soils
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
          'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_deter_fun <- function(dist) {
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


similar_aqui_data <- similar_deter_fun(aggre_trait_aqui_dist)
similar_use_data <- similar_deter_fun(aggre_trait_use_dist)
similar_stress_data <- similar_deter_fun(aggre_trait_stress_dist)
  
similar_gene_trait_data <- data.frame(Group = similar_aqui_data$Group,
                           Gully_id = similar_aqui_data$Gully_id, 
                           distance_aqui = similar_aqui_data$distance,
                           distance_use = similar_use_data$distance,
                           distance_stress = similar_stress_data$distance)


# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
similar_gene_traits_df <- similar_gene_trait_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
dis_index <- c("distance_aqui", "distance_use", "distance_stress")
lmm_fun <- function(vars, df) {
  library(lme4)
  library(lmerTest)
  lmm_dist_modes <- lapply(vars, function(x) {
    lmer(substitute(i ~ Group + Time + Slope + MAP + (1 | Gully_id), list(i = as.name(x))), 
         data = df)})
  summary.model <- function(model){
    F.value <- anova(model)$'F value'
    p.value <- anova(model)$'Pr(>F)'
    p.stars <- function(p.values) {
      unclass(symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "")))}
    sig <- p.stars(p.value)
    results<-data.frame(F.value, p.value, sig)
    return(results)
  }
  df <- NULL
  for(i in 1:length(vars)) {
    tmp <- summary.model(lmm_dist_modes[[i]])
    if (is.null(df)){
      df <- tmp
    } else {
      df <- rbind(df, tmp)
    }
  }
  result_lmm <-data.frame(dist_index = rep(vars, each = 4), 
                          variables = rep(c("Group", "Time", "Slope", "MAP"), 
                                          length(vars)), df)
  return(result_lmm)
}

lmm.results <- lmm_fun(dis_index, similar_gene_traits_df)
lmm.results


# Add the significant symbols manually
sig.dis.labs <- tibble(dis_index = factor(dis_index, levels = dis_index),
                       x1 = rep(0.5, length(dis_index)),
                       y1 = rep(0.95, length(dis_index)),
                       sig.labels = lmm.results %>% filter(variables == "Group") %>%
                         select(sig) %>% pull(),
                       x2 = rep(0.1, length(dis_index)),
                       y2 = rep(1, length(dis_index)),
                       panel.labels = letters[as.numeric(dis_index)])

# Create a box plot
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
library(gghalves)
dis_gene_traits_plot <- similar_gene_traits_df %>% 
  select(c("Group", dis_index)) %>%
  gather(dis_index, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(dis_index = factor(dis_index, 
                            levels = c("distance_aqui", "distance_use", "distance_stress"))) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::geom_text_npc(data = sig.dis.labs, aes(npcx = x1, npcy = y1, label = sig.labels), inherit.aes = F) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~dis_index, scales = "free_y", ncol = 1) +
  main_theme +
  theme(legend.position = "none")
dis_gene_traits_plot
