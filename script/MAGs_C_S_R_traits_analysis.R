################################################
######### Metagenomic assembled genomes ########
################################################
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
# bin genome information
mag_results_trans <- decostand(t(mags_abun_tab), method = "hellinger")
mag_dist <-vegdist(mag_results_trans, "bray")
#permanova test the difference in compositional variance
adonis2(mag_dist ~ Group, data = metadata)

#Visualization for the overall difference by PCoA plot.
#mag
pcoa_mag_plot <- PCoA_plot_fun(mag_dist)

pcoa_mag_plot


## Explore the difference in MAGs taxonomic variance between uncollapsed and collapsed soils
distance_mag_data <- similar_deter_fun(mag_dist)
# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
distance_mags_df <- distance_mag_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
lmm_dis_mags_mod <- lmm_fun("distance", distance_mags_df)
# Create a box plot
library(gghalves)
dis_mag_plot <- distance_mags_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.5, npcy = 0.95, 
                 size = 2, label = "P < 0.001") +
  labs(x = "Group", y = "Dissimilarity in MAGs composition") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")
dis_mag_plot


################################################
##### Metagenomic assembled genome traits ######
################################################
# Comparing all microbial traits
library(microtrait)
microtrait_results <- readRDS(file.path(wd_fun, "/MAGs/microtraits/thermokarst_gully.microtraitresults.rds"))

# Normalizing the traits matrices by genome length
microtrait_results_metadata_norm = microtrait_results %>% trait.normalize(normby = "genome_length")

spec_names <- data.frame(microtrait_results_metadata_norm$trait_matrixatgranularity1)$id
trait_data <- data.frame(microtrait_results_metadata_norm$trait_matrixatgranularity1)[,-1]
rownames(trait_data) <- spec_names

# Calculate CWM for each trait and each community
trait_data <- trait_data[rownames(mags_abun_tab), ]
relative_abundance <- mags_abun_tab/100
cwm_results <- data.frame(matrix(NA, nrow = ncol(relative_abundance), 
                                 ncol = ncol(trait_data)))
row.names(cwm_results) <- colnames(relative_abundance)
colnames(cwm_results) <- colnames(trait_data)

for (i in 1:ncol(trait_data)) {
  for (j in 1:ncol(relative_abundance)) {
    cwm_results[j, i] <- sum(relative_abundance[, j] * trait_data[, i], na.rm = TRUE)/sum(relative_abundance[, j])
  }
}

print(cwm_results)

## Test the overall difference in the community weighted mean traits.
library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
cwm_results_final <- cwm_results[ ,colSums(cwm_results[])>0]
cwm_results_final <- cwm_results_final[metadata$Sample_name, ]

dim(cwm_results_final)
# cwm_results_trans <- scale(cwm_results_final, center = TRUE, scale = T)
cwm_results_trans <- decostand(cwm_results_final, method = "log")
trait_dist <-vegdist(cwm_results_trans, "bray", binary = F)
#permanova test the difference in compositional variance
adonis2(trait_dist ~ Group, data = metadata)

#Visualization for the overall difference by PCoA plot.
pcoa_cwm_trait_plot <- PCoA_plot_fun(trait_dist)
pcoa_cwm_trait_plot

# Extract the distance between un-collapsed and collapsed samples
distance_trait_data <- similar_deter_fun(trait_dist)

# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
distance_trait_df <- distance_trait_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
lmm_dis_traits_mod <- lmm_fun("distance", distance_trait_df)
lmm_dis_traits_mod

# Create a box plot
library(gghalves)
dis_trait_plot <- distance_trait_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.5, npcy = 0.95, 
                 size = 2, label = "P < 0.001") +
  labs(x = "Group", y = "Dissimilarity in CWM traits") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")
dis_trait_plot

# Linear mixed models test the effect of collapsed and gully_id
trait_id <- colnames(cwm_results)
trait_scale <- cwm_results %>% 
  cbind(metadata[, c("Gully_id","Group", "Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", trait_id))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
trait_S1 <- sapply(6:ncol(trait_scale), function(j) {
  if (length(unique(trait_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(trait_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = trait_scale)
    
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
colnames(trait_S1)<-colnames(trait_scale)[-c(1:5)]
data.frame(trait_S1)

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", "")))}
# Create a plot
each_trait_comparison <- trait_S1 %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>% filter(variables %in% trait_id) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = trait_id)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, 
                y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) 
                                                                    + GroupCollapsed.se)*1.25),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2.5, 2.5)) +
  scale_x_discrete(position = "bottom") +
  annotate("rect", xmin = 0.5, xmax = 16.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 16.5, xmax = 18.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 18.5, xmax = 22.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#5cc3e8") +
  annotate("rect", xmin = 22.5, xmax = 29.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#ffdb00") +
  annotate("rect", xmin = 29.5, xmax = 37.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 37.5, xmax = 40.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#5fb236") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
each_trait_comparison

# Combine all plots
library(cowplot)
dist_gmags_trait_plot <- plot_grid(pcoa_mag_plot, dis_mag_plot,
                                  pcoa_cwm_trait_plot, dis_trait_plot,
                                  nrow = 1, align = "hv")
dist_gmags_trait_plot

mags_plots <- plot_grid(dist_gmags_trait_plot, each_trait_comparison,
                        nrow = 2, rel_heights = c(1, 2.8))
mags_plots

# Save the combined plot
# ggsave(file.path("E:/thermokarst_gully/revision/result/mags_plots.pdf"), 
#        mags_plots, width = 183, height = 155, units = "mm")












# Functional dispersion calculation
trait_fd_list <- list()

for(trait_name in colnames(trait_data)) {
  trait_fd_list[[trait_name]] <- dbFD(
    x = trait_data[, trait_name, drop = FALSE],
    a = t(mags_abun_tab),
    messages = FALSE
  )
}

# Extract FDis of each trait
fdis_per_trait <- sapply(trait_fd_list, function(x) x$FDis)

# Test the permafrsot collase effect on Functional diversity
fdis_scale <- fdis_per_trait %>% 
  cbind(metadata[, c("Gully_id","Group", "Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", trait_id))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
fdis_S1 <- sapply(6:ncol(fdis_scale), function(j) {
  if (length(unique(fdis_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(fdis_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = fdis_scale)
    
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
colnames(fdis_S1)<-colnames(fdis_scale)[-c(1:5)]

# Add significance symbol
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
all_fdis_comparison <- fdis_S1 %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>% filter(variables %in% trait_id) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(trait_id))) %>%
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
  coord_flip() + scale_x_discrete(position = "bottom") +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 3.5, xmax = 11.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 11.5, xmax = 18.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#5cc3e8") +
  annotate("rect", xmin = 18.5, xmax = 40.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#ffdb00") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )

# ggsave(file.path("E:/thermokarst_gully/result2/all_tax_div_comparison.pdf"),
#        all_tax_div_comparison, width = 2.5, height = 5, units = "in")
all_fdis_comparison
