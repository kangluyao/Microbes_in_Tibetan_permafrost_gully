# Set work directory
setwd('e:/thermokarst_gully/')
wd_16s <- file.path(getwd(),"data/16S/rdp")
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")

# Loading packages
pacman::p_load(vegan, tidyverse, ggrepel, EnhancedVolcano, edgeR, ggplot2)

# Data input
source("script/read_data.R")
source("script/Functions.R")
# Permutational multivariate analysis of variance using distance matrices (adonis) to test the difference of functional composition
library(vegan)
# determine the dissimilarity matrix based on the bray-curties distance
ko_dist <-vegdist(t(ko_tpm_table), "bray" )
# permanova, ANOSIM and MRPP analysis
adonis2(ko_dist ~ Group, data = metadata)
mrpp(ko_dist, metadata$Group, perm = 999)
anosim(ko_dist, metadata$Group, perm = 999)

# PCoA plot
ko_pcoa_plot <- PCoA_plot_fun(ko_dist)

distance_ko_data <- similar_deter_fun(ko_dist)

# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
distance_ko_df <- distance_ko_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
var = "distance"
lmm_dis_ko_mod <- lmm_fun(var, distance_ko_df)
lmm_dis_ko_mod

# Create a box plot
library(gghalves)
dis_ko_plot <- distance_ko_df %>% 
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
dis_ko_plot

## Explore the losted genes after permafrost thawing
library(ggvenn)
x <- list(
  `Un-collapsed` = ko_tpm_table %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_C', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames(), 
  Collapsed = ko_tpm_table %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_T', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames()
)
# plot
ko.venn <- ggvenn(
  x, 
  fill_color = c("#f8766d", "#a3a500", "#00b0f6"),
  stroke_color = NA,
  set_name_size = 4,
  text_size = 4,
  show_percentage = F
)
ko.venn

# loading packages
library(tidyverse)
#reading data
ko_9_levels_file <- 'E:/thermokarst_gully/data/metagenome/ko_9_levels.txt'
ko_9_levels <- fread(ko_9_levels_file, header = TRUE, stringsAsFactors = FALSE)
YAS_sel_L3 <- fread("E:/thermokarst_gully/data/YAS.csv")

KO_L3_df <- ko_tpm_table %>%
  rownames_to_column("KO") %>%
  left_join(ko_9_levels, by = "KO") %>%
  select(KO, L3, starts_with("G")) %>%
  distinct() %>%
  filter(!is.na(L3)) %>%
  inner_join(YAS_sel_L3, by ="L3")
# Test the effect size of permafrost collapse on the YAS functional pathways


YAS_pathway_scale <- KO_L3_df %>%
  select(L3, starts_with("G")) %>%
  group_by(L3) %>%
  summarise(across(where(is.numeric), sum)) %>%
  filter(L3 != "") %>% 
  column_to_rownames("L3") %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample_id") %>% 
  mutate(Gully_id = sapply(stringr::str_split(Sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(Group = case_when(grepl("_C", Sample_id) ~ "Uncollapsed",
                           grepl("_T", Sample_id) ~ "Collapsed")) %>%
  mutate(across(where(is.numeric), scale)) %>%
  relocate(where(is.character)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(-Sample_id)

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
YAS_pathway_S1 <- sapply(3:ncol(YAS_pathway_scale), function(j) {
  if (length(unique(YAS_pathway_scale[, j])) < 3) {
    result <- rep(NA, 40)
  } else {
    fm1 <- lmer(YAS_pathway_scale[, j] ~ Group + (1 | Gully_id), data = YAS_pathway_scale)
    
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
colnames(YAS_pathway_S1)<-colnames(YAS_pathway_scale)[-c(1:2)]
YAS_pathway_S1[1:6, 1:6]

pathway_index <- YAS_sel_L3$L3
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*",  "")))}
YAS_diff_in_pathway_plot <- YAS_pathway_S1 %>%
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
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se + 0.2)),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") + 
  scale_colour_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0, xmax = 25.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#ff6666") +
  annotate("rect", xmin = 25.5, xmax = 38.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#79ceb8") +
  annotate("rect", xmin = 38.5, xmax = 56.5, ymin = -2.5, ymax = 2.5, alpha = 0.2, fill = "#e56eee") +
  main_theme +
  theme(strip.background = element_rect(fill = c("#FFF6E1")), legend.position = "none")

# ggsave(file.path("E:/thermokarst_gully/revision/result/diff_in_pathway_plot.pdf"),
#        diff_in_pathway_plot, width = 4.5, height = 16, units = "in")
YAS_diff_in_pathway_plot


# Explore the permafrost effects on the Y-A-S gene groups

aggre_trait_YAS <- KO_L3_df %>%
  filter(!is.na(L3)) %>%
  select(Strategy, starts_with("G")) %>%
  group_by(Strategy) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Strategy") %>%
  t() %>% as.data.frame() %>%
  # Please ensure your row order same
  cbind(metadata[, c("Gully_id", "Group", "MAP", "Time", "Slope")]) %>%
  as.data.frame()

yas_vars <- c("Growth yield", "Resource acquisition", "Stress tolerance")
yas_gene.abundance_results <- lmm_fun(yas_vars, aggre_trait_YAS)           


# Add the significant symbols manually
sig.yas_gene.labs <- tibble(yas_vars = factor(yas_vars, levels = yas_vars),
                            x1 = rep(0.5, length(yas_vars)),
                            y1 = rep(0.95, length(yas_vars)),
                            sig.labels = yas_gene.abundance_results %>% 
                              filter(variables == "Group") %>%
                              select(sig) %>% pull(),
                            x2 = rep(0.1, length(yas_vars)),
                            y2 = rep(1, length(yas_vars)),
                            panel.labels = letters[as.numeric(yas_vars)])

# Create a box plot
library(gghalves)
yas_gene.abundance_plot <- aggre_trait_YAS %>% 
  select(c("Group", yas_vars)) %>%
  gather(yas_vars, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(yas_vars = factor(yas_vars, levels = c("Growth yield", 
                                                "Resource acquisition", 
                                                "Stress tolerance"))) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::geom_text_npc(data = sig.yas_gene.labs, 
                      aes(npcx = x1, npcy = y1, label = sig.labels), 
                      inherit.aes = F) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~yas_vars, scales = "free_y", ncol = 1) +
  main_theme +
  theme(legend.position = "none")
yas_gene.abundance_plot

# Now exploring the permafrost collapse effects on gene compositions
aggre_trait_Y <- KO_L3_df %>%
  filter(!is.na(L3) & Strategy == "Growth yield") %>%
  select(KO, starts_with("G")) %>%
  group_by(KO) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("KO")

aggre_trait_A <- KO_L3_df %>%
  filter(!is.na(L3) & Strategy == "Resource acquisition") %>%
  select(KO, starts_with("G")) %>%
  group_by(KO) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("KO")

aggre_trait_S <- KO_L3_df %>%
  filter(!is.na(L3) & Strategy == "Stress tolerance") %>%
  select(KO, starts_with("G")) %>%
  group_by(KO) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("KO")


library(vegan)
# Transform the gene abundance data
aggre_trait_Y_trans <- decostand(t(aggre_trait_Y), method = "hellinger")
aggre_trait_A_trans <- decostand(t(aggre_trait_A), method = "hellinger")
aggre_trait_S_trans <- decostand(t(aggre_trait_S), method = "hellinger")
aggre_trait_Y_dist <-vegdist(aggre_trait_Y_trans, "bray")
aggre_trait_A_dist <-vegdist(aggre_trait_A_trans, "bray")
aggre_trait_S_dist <-vegdist(aggre_trait_S_trans, "bray")

#permanova test the difference in compositional variance
adonis2(aggre_trait_Y_dist ~ Group, data = metadata)
adonis2(aggre_trait_A_dist ~ Group, data = metadata)
adonis2(aggre_trait_S_dist ~ Group, data = metadata)

# PCoA plot
PCoA_plot_fun <- function(dist) {
  ord.fun <-  cmdscale(dist,  k = 2, eig = T, add = T)
  pcoa.plot <- data.frame(Group = metadata$Group, scores(ord.fun)) %>%
    mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
    ggplot(aes(x = Dim1, y = Dim2, colour = Group)) + 
    geom_point(size = 1, alpha = 0.8, shape = 16) + 
    stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1, 
                 show.legend = FALSE, level = 0.95) +
    scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    # scale_x_continuous(expand = c(0.03, 0.03)) +
    # scale_y_continuous(expand = c(0.03, 0.03)) +
    labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), digits = 3), "%)", sep = "")) +
    main_theme
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
trait_Y.plot <- PCoA_plot_fun(aggre_trait_Y_dist)
trait_A.plot <- PCoA_plot_fun(aggre_trait_A_dist)
trait_S.plot <- PCoA_plot_fun(aggre_trait_S_dist)
# Combine all PCoA plots
PCoA_YAS_all_plots <- cowplot::plot_grid(trait_Y.plot, trait_A.plot, 
                                         trait_S.plot, ncol = 1)
PCoA_YAS_all_plots


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


similar_Y_data <- similar_deter_fun(aggre_trait_Y_dist)
similar_A_data <- similar_deter_fun(aggre_trait_A_dist)
similar_S_data <- similar_deter_fun(aggre_trait_S_dist)

similar_YAS_trait_data <- data.frame(Group = similar_Y_data$Group,
                                      Gully_id = similar_Y_data$Gully_id, 
                                      distance_Y = similar_Y_data$distance,
                                      distance_A = similar_A_data$distance,
                                      distance_S = similar_S_data$distance)


# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
similar_YAS_traits_df <- similar_YAS_trait_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
dis_index <- c("distance_Y", "distance_A", "distance_S")
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

lmm.results <- lmm_fun(dis_index, similar_YAS_traits_df)
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
library(gghalves)
dis_YAS_traits_plot <- similar_YAS_traits_df %>% 
  select(c("Group", dis_index)) %>%
  gather(dis_index, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(dis_index = factor(dis_index, 
                            levels = c("distance_Y", "distance_A", "distance_S"))) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::geom_text_npc(data = sig.dis.labs, aes(npcx = x1, npcy = y1, label = sig.labels), inherit.aes = F) +
  labs(x = "Group", y = "Dissimilarity in gene composition") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~dis_index, scales = "free_y", ncol = 1) +
  main_theme +
  theme(legend.position = "none")
dis_YAS_traits_plot

# Combine all plots
library(cowplot)
dist_YAS.gene_trait_plot <- plot_grid(YAS_diff_in_pathway_plot, 
                                      PCoA_YAS_all_plots, 
                                      dis_YAS_traits_plot,
                                      nrow = 1, rel_widths = c(2, 1, 1))
dist_YAS.gene_trait_plot

ggsave(file.path("E:/thermokarst_gully/revision/result/dist_YAS.gene_trait_plot.pdf"),
       dist_YAS.gene_trait_plot, width = 183, height = 120, units = "mm")
