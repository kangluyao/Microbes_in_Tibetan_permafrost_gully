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


# Determine the alpha diversity
mag_comp = as.matrix(t(mags_abun_tab))
mag_shannon  <- diversity(mag_comp, index = "shannon")
mag_richness <- specnumber(mag_comp)

alpha_div_mag <- data.frame(
  Richness = mag_richness, Shannon = mag_shannon) %>%
  mutate(Evenness = Shannon / log(Richness))

# Arrange the data for statistics
col_div_names <- c("Sample_name", 'Gully_id', 'Group', 
                   "Richness", "Shannon", "Evenness")
div_table_mag <- cbind(metadata[, c("Sample_name", 'Gully_id', 'Group')], alpha_div_mag) %>%
  rename_with(~ col_div_names) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

div_table_mag %>%
  group_by(Gully_id, Group) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = T), 
                        sd = ~sd(., na.rm = T), 
                        n = ~n())))

# Linear mixed models test the effect of collapsed and gully_id
div_index <- c("Richness", "Shannon", "Evenness")
div_scale_mag <- div_table_mag %>% 
  cbind(metadata[, c("Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", div_index))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
library(lme4)
library(lmerTest)
div_S1_mag <- sapply(6:ncol(div_scale_mag), function(j) {
  if (length(unique(div_scale_mag[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(div_scale_mag[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = div_scale_mag)
    
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
colnames(div_S1_mag)<-colnames(div_scale_mag)[-c(1:5)]
data.frame(div_S1_mag)

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", "")))}

# Add the significant symbols manually
sig.div.labs <- tibble(div_index = factor(div_index, levels = div_index),
                       x1 = rep(0.5, length(div_index)),
                       y1 = rep(0.95, length(div_index)),
                       sig.labels = c("P < 0.001", "P < 0.001", "P = 0.101"),
                       panel.labels = letters[as.numeric(div_index)])

# Create a plot
library(gghalves)
div_mag_plot <- div_table_mag %>% select(c("Group", div_index)) %>%
  gather(div_index, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(div_index = factor(div_index, 
                            levels = c("Richness", "Shannon", "Evenness"))) %>%
  ggplot(aes(Group, value)) +
  geom_half_violin(aes(fill = Group), color = NA, 
                   position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, alpha = 0.65) +
  geom_boxplot(aes(fill = Group), width = 0.35, 
               size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(colour = Group), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::geom_text_npc(data = sig.div.labs, aes(npcx = x1, npcy = y1, 
                                               label = sig.labels), 
                      size = 2, inherit.aes = F) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~div_index, scales = "free_y", ncol = 3) +
  main_theme +
  theme(legend.position = "none")
div_mag_plot
# bin genome information
mag_results_trans <- decostand(mag_comp, method = "hellinger")
mag_dist <-vegdist(mag_results_trans, "bray")
#permanova test the difference in compositional variance
set.seed(123)
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
library(lme4)
library(lmerTest)
lmm_dis_mags_mod <- lmer(distance ~ Group + Time + 
                        Slope + MAP + (1 | Gully_id), data = distance_mags_df)

car::Anova(lmm_dis_mags_mod)
# Create a box plot
library(gghalves)
dis_mag_plot <- distance_mags_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 16, size = 0.5,
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

# --- β dispersion ---
mag_comm <- t(mags_abun_tab)
treatment <- metadata$Group
block <- metadata$Gully_id

# --- 1) Beta-diversity decomposition ---
# For abundance data
library(betapart)
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
set.seed(123)
pt_mag_total <- permutest(bd_mag_total, permutations = ctrl_how)
pt_mag_turn  <- permutest(bd_mag_turn, permutations = ctrl_how)
pt_mag_nest  <- permutest(bd_mag_nest, permutations = ctrl_how)

# Inspect results
pt_mag_total
pt_mag_turn
pt_mag_nest

# --- 4) Visualization: boxplots of sample dispersion values ---
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
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~disp, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               `disp_total` = "Total dispersion",
               `disp_turn` = "Turnover dispersion",
               `disp_nest` = "Nestedness dispersion"
             ))) +
  main_theme +
  theme(legend.position = "none")
disp_mag_plot
# ggsave(file.path("E:/thermokarst_gully/revision/result/disp_mag_plot.pdf"),
#        disp_mag_plot, width = 120, height = 45, units = "mm", dpi = 600)

################################################
##### Metagenomic assembled genome traits ######
################################################
# Comparing all microbial traits
library(microtrait)
microtrait_results <- readRDS(file.path(wd_fun, "/MAGs/microtraits/thermokarst_gully.microtraitresults.rds"))

# Normalizing the traits matrices by genome length
microtrait_results_metadata_norm = microtrait_results %>% 
  trait.normalize(normby = "genome_length")

# Inputting other traits including genome size, max growth rate, and so on

mags_result_eff.size <- mags_result_lmm %>% rownames_to_column("Genome_id")
growth_rate_file <- "E:/thermokarst_gully/data/metagenome/MAGs/growth_rate_results.csv"
growth_rate <- read.table(growth_rate_file, header = TRUE, sep = ",", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE)
MAGs_all_aatr <- mags.att %>% left_join(growth_rate, by = "Genome_id") %>% 
  left_join(mags_result_eff.size, by = "Genome_id") %>%
  mutate(Genome_size = round(Genome_size/1000000, 1), `Max growth rate` = 1/d) %>%
  select(Genome_id, Genome_size, GC, `Max growth rate`)

# Combining all MAGs attributes
trait_data = microtrait::add.metadata(microtrait_results_metadata_norm, 
                                                       MAGs_all_aatr, genome_metadata_idcol = "Genome_id")
spec_names <- data.frame(trait_data$trait_matrixatgranularity1)$id
trait_data <- data.frame(trait_data$trait_matrixatgranularity1, check.names=FALSE)[,-c(1,39:41)]
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
set.seed(123)
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
lmm_dis_trait_mod <- lmer(distance ~ Group + Time + 
                           Slope + MAP + (1 | Gully_id), data = distance_trait_df)

car::Anova(lmm_dis_trait_mod)

# Create a box plot
library(gghalves)
dis_trait_plot <- distance_trait_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, 
                   color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 16, size = 0.5,
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

# --- β dispersion ---
trait_comm <- cwm_results_trans
treatment <- metadata$Group
block <- metadata$Gully_id

# --- 1) Beta-diversity decomposition ---
# For abundance data
bpair_trait_abund <- beta.pair.abund(trait_comm, index.family = "bray")

# Convert to dist objects for betadisper()
dist_trait_total <- as.dist(bpair_trait_abund$beta.bray)
dist_trait_turn  <- as.dist(bpair_trait_abund$beta.bray.bal)
dist_trait_nest  <- as.dist(bpair_trait_abund$beta.bray.gra)

# --- 2) Calculate dispersion (distance to group centroid) for each treatment ---
bd_trait_total <- betadisper(dist_trait_total, treatment)
bd_trait_turn  <- betadisper(dist_trait_turn, treatment)
bd_trait_nest  <- betadisper(dist_trait_nest, treatment)

# --- 3) Permutation test with block-constrained randomization ---
# Create permutation control object
# Permutations are restricted within each block (so that treatment labels are shuffled only within blocks)
ctrl_how <- how(nperm = 999, blocks = block, within = Within(type = "free"))

# Run permutational tests
pt_trait_total <- permutest(bd_trait_total, permutations = ctrl_how)
pt_trait_turn  <- permutest(bd_trait_turn, permutations = ctrl_how)
pt_trait_nest  <- permutest(bd_trait_nest, permutations = ctrl_how)

# Inspect results
pt_trait_total
pt_trait_turn
pt_trait_nest

# --- 4) Visualization: boxplots of sample dispersion values ---
df_trait_disp <- data.frame(
  sample = names(bd_trait_total$distances),
  disp_total = bd_trait_total$distances,
  disp_turn  = bd_trait_turn$distances,
  disp_nest  = bd_trait_nest$distances,
  treatment = treatment,
  block = block
) %>% select(disp_total, disp_turn, disp_nest, treatment) %>%
  pivot_longer(-treatment, names_to = "disp", values_to = "value") %>%
  mutate(treatment = factor(treatment, levels = c("Un-collapsed", "Collapsed"))) 

# Plot dispersion (divergence) distributions
library(gghalves)
disp_index <- c("disp_total", "disp_turn", "disp_nest")
disp_trait_plot <- df_trait_disp %>% 
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
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~disp, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               `disp_total` = "Total dispersion",
               `disp_turn` = "Turnover dispersion",
               `disp_nest` = "Nestedness dispersion"
             ))) +
  main_theme +
  theme(legend.position = "none")
disp_trait_plot
# ggsave(file.path("E:/thermokarst_gully/revision/result/disp_trait_plot.pdf"),
#        disp_trait_plot, width = 120, height = 45, units = "mm", dpi = 600)



# Linear mixed models test the effect of collapsed and gully_id
strategy_vars <- c("Genome_size", "Max growth rate", 
                   "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization",
                   "Resource Acquisition:Substrate degradation:simple compound degradation",
                   "Resource Acquisition:Substrate uptake:carbohydrate transport",
                   "Resource Acquisition:Substrate uptake:N compound transport",
                   "Resource Acquisition:Substrate uptake:organophosphorus transport",
                   "Resource Acquisition:Substrate assimilation:C1 compounds",
                   "Resource Acquisition:Substrate assimilation:N compounds",
                   "Resource Acquisition:Substrate assimilation:P compounds",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration",
                   "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation",
                   "Resource Use:Chemotrophy:chemolithoautotrophy",
                   "Resource Use:Phototrophy:photosystem",
                   "Stress Tolerance:Specific:temperature:low",
                   "Stress Tolerance:Specific:pH stress",
                   "Stress Tolerance:Specific:desiccation/osmotic/salt stress",
                   "Stress Tolerance:Specific:oxidative stress",
                   "Stress Tolerance:Specific:envelope stress")

# writing a function to change the facet labels
strategy_names <- list(
  "Genome_size" = "Genome size (Mb)",
  "Max growth rate" = "Max growth rate",
  "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization" = "Complex carbohydrate depolymerization",
  "Resource Acquisition:Substrate degradation:simple compound degradation" = "Simple compound degradation",
  "Resource Acquisition:Substrate uptake:carbohydrate transport" = "Carbohydrate transport",
  "Resource Acquisition:Substrate uptake:N compound transport" = "N compound transport",
  "Resource Acquisition.Substrate uptake:organophosphorus transport" = "organophosphorus transport",
  "Resource Acquisition:Substrate assimilation:C1 compounds" = "C1 compounds assimilation",
  "Resource Acquisition:Substrate assimilation:N compounds" = "N compounds assimilation",
  "Resource Acquisition.Substrate assimilation:P compounds" = "P compounds assimilation",
  "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration" = "Aerobic respiration", 
  "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration" = "Anaerobic respiration", 
  "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation" = "Fermentation",
  "Resource Use:Chemotrophy:chemolithoautotrophy" = "Chemolithoautotrophy",
  "Resource Use:Phototrophy:photosystem" = "Photosystem",
  "Stress Tolerance:Specific:temperature:low" = "Low temperature stress",
  "Stress Tolerance:Specific:pH stress" = "pH stress",
  "Stress Tolerance:Specific:desiccation/osmotic/salt stress" = "Desiccation/osmotic/salt stress", 
  "Stress Tolerance:Specific:oxidative stress" = "Oxidative stress",
  "Stress Tolerance:Specific:envelope stress" = "Envelope stress")

trait_id <- strategy_vars
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
  scale_x_discrete(position = "bottom", labels = strategy_names) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#79ceb8") +
  annotate("rect", xmin = 2.5, xmax = 10.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#e95f5c") +
  annotate("rect", xmin = 10.5, xmax = 15.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#5cc3e8") +
  annotate("rect", xmin = 15.5, xmax = 20.5, ymin = -2.5, ymax = 2.5, alpha = 0.1, fill = "#ffdb00") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
each_trait_comparison

# Combine all plots
library(cowplot)
dist_mags_trait_plot <- plot_grid(pcoa_mag_plot, dis_mag_plot,
                                  pcoa_cwm_trait_plot, dis_trait_plot,
                                  nrow = 1, align = "hv")
dist_mags_trait_plot

mags_plots <- plot_grid(dist_mags_trait_plot, each_trait_comparison,
                        nrow = 2, rel_heights = c(1, 2.8))
mags_plots

# Save the combined plot
ggsave(file.path("E:/thermokarst_gully/revision/result/mags_plots.pdf"),
       mags_plots, width = 183, height = 155, units = "mm")












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
