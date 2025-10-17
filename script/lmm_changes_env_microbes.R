rm(list = ls())
# Data input 
## Set working and saving directory
setwd('e:/thermokarst_gully/')
save.dir <- file.path(getwd(),"result")

# Loading packages and read data
pacman::p_load(phyloseq, ape, vegan, Biostrings, microbiome, tidytable, tidyverse, rstatix)
source("script/read_data.R")

# Test the difference in the environmental variables between the uncollapsed and collapsed sites
env_vars <- c("Plant_richness", "AGB", "pH", "Soil_moisture", "Clay_Silt", "SOC",
              "NH4_N", "NO3_N", "AP")
env_stats <- metadata %>% group_by(Group) %>%
  get_summary_stats(env_vars, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  arrange(variable, Group)
env_stats # the discriptive statistics 

# LMM
library(lme4)
library(lmerTest)
lmm_env_modes <- lapply(env_vars, function(x) {
  lmer(substitute(i ~ Group + Time + Slope + MAP + (1 | Gully_id), list(i = as.name(x))), data = metadata)})
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
for(i in 1:length(env_vars)) {
  tmp <- summary.model(lmm_env_modes[[i]])
  if (is.null(df)){
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}

env_result_lmm <-data.frame(Variables = rep(env_vars, each = 4),
                            fixed_factors = rep(c("Group", "Time", "Slope", "MAP"), length(env_vars)),
                            group1 = rep("Un-collapsed", length(env_vars)),
                            group2 = rep("Collapsed", length(env_vars)), df)
env_result_lmm


# Bar plot
rep_str1 = list("Plant_richness" = "Plant richness",
                "AGB" = expression(paste("AGB ", "(g/", "m"^2, ")")),
                "pH" = "pH",
                "Soil_moisture" = "Moisture (%)",
                "Clay_Silt" = "Clay+Silt (%)",
                "SOC" = "SOC (g/kg)",
                "NH4_N" = expression(paste("NH"[4]^"+", "-N", " (mg/kg)")),
                "NO3_N" = expression(paste("NO"[3]^"-", "-N", " (mg/kg)")),
                "AP" = "AP (mg/kg)"
)
# write a function to change the facet labels
facet_labeller <- function(variable,value){
  return(rep_str1[value])
}

# Create a bar plot
env_plot <- env_stats %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(variable = factor(variable, levels = env_vars)) %>%
  ggplot(aes(Group, mean, fill = Group)) +
  geom_bar(position = position_dodge(0.8), stat = 'identity', 
           colour = 'black', width = 0.6, )+
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = .2) +
  # geom_text(aes(y = 1.05*(mean + se), label = sig.lab), 
  #           position = position_dodge(width = .75), size = 2.5, vjust = 0.01) + # add lowercase letters to indicate the significance
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~variable, scales = "free", ncol = 3, strip.position = "left", labeller = facet_labeller) +
  theme_classic() + 
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(colour = NA),
        strip.text = element_text(colour = 'black', size = 7, margin = margin()),
        axis.title = element_text(color = 'black',size = 7),
        axis.text = element_text(colour = 'black', size = 7),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        legend.title = element_text(colour = 'black', size = 7),
        legend.text = element_text(colour = 'black', size = 7),
        legend.key.size = unit(0.5, 'cm'))

lines <- tibble(variable = factor(c("Plant_richness", "AGB", "pH", "Soil_moisture", "Clay_Silt", "SOC", "NH4_N", "NO3_N", "AP"), levels = env_vars),
                x = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
                xend = c(2, 2, 2, 2, 2, 2, 2, 2, 2),
                y = c(17, 270, 8, 220, 72, 185, 19, 40, 2.4),
                yend = c(17, 270, 8, 220, 72, 185, 19, 40, 2.4)
)
stars <- tibble(variable = factor(c("Plant_richness", "AGB", "pH", "Soil_moisture", "Clay_Silt", "SOC", "NH4_N", "NO3_N", "AP"), levels = env_vars),
                x1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                y1 = c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95),
                sig.labels = env_result_lmm %>% filter(fixed_factors == "Group") %>%
                  select(sig) %>% pull,
                x2 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                y2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
                panel.labels = letters[as.numeric(variable)]
)

env_plot <- env_plot +
  # geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE) + # add a sigment to denote the comparison between groups
  ggpp::geom_text_npc(data = stars, aes(npcx = x1, npcy = y1, label = sig.labels), 
                      size = 3, inherit.aes = F) +
  ggpp::geom_text_npc(data = stars, aes(npcx = x2, npcy = y2, label = panel.labels), 
                      size = 3, inherit.aes = F)

# save the plot
# ggsave(file.path(save.dir, "./figs/env_effect/env_comparison_plot.pdf"),
#        env_plot, width = 130, height = 90, units = "mm")
env_plot

# Explore the effect of permafrost collapse on the environmental heterogeneity.
env_df <- metadata %>% select(env_vars)
env.dist <- vegdist(scale(env_df), 'euclidean', na.rm = T)
# difference in environmental variance among uncollapsed and collapsed sample at plot level
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# extract the values within site (at plot level) of uncollapsed and collapsed samples, respectively.
env_dis_plot_table <- lapply(vars, function(x) usedist::dist_subset(env.dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
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

# merge similar_all_data with meta_unique
env_dis_plot_df <- env_dis_plot_table %>%
  left_join(meta_unique, by = "Gully_id")


# Linear mixed model test the difference of environmental heterogeneity among uncollapsed and collapsed samples
library(lme4)
library(lmerTest)
summary.model <- function(model){
  F.value <- anova(model)$`F value`
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
env.dist_mode <- lmer(distance ~ Group + Time + Slope + MAP + (1 | Gully_id), 
                      data = env_dis_plot_df)
env_lmm_model_results <- summary.model(env.dist_mode)
env_lmm_model_results

main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 6),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.title = element_text(colour = 'black', size = 6),
        legend.text = element_text(colour = 'black', size = 5),
        legend.key.size = unit(0.5, 'cm'))
# Create a box plot
library(gghalves)
env_dist_plot <- env_dis_plot_df %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = Group), shape = 21, size = 1.5, width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  ggpp::geom_text_npc(aes(npcx = 0.5, npcy = 0.95, label = "**")) +
  # geom_text(aes(x = 1.5, y = Inf, label = "**"), inherit.aes = F) +
  labs(x = NULL, y = "Environmental heterogeneity") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")

# save the plot
# ggsave(file.path(save.dir, "./figs/env_effect/environmental_heterogeneity_plot.pdf"),
#        env_dist_plot, width = 45, height = 45, units = "mm")
env_dist_plot

# Test the collinearity of the environmental variables using corrplot
# Load required libraries
library(corrplot)
library(dplyr)

# Define environmental variables
env_vars <- c("MAT", "MAP", "Time", "Slope", "Plant_richness", "AGB", "pH", "Soil_moisture", "Clay_Silt", 
              "SOC", "NH4_N", "NO3_N", "AP")

# Extract the environmental variables from metadata
env_data <- metadata[, env_vars]

# Calculate correlation matrix
cor_matrix <- cor(env_data, method = "spearman")

# Print correlation matrix
print(cor_matrix)

# Create correlation plot
corrplot(cor_matrix,
         type = "lower",             # Display type: "full", "lower", "upper"
         order = "original",           # Ordering method: "original", "AOE", "FPC", "hclust", "alphabet"
         tl.cex = 0.6,              # Text label size
         tl.col = "black",          # Text label color
         tl.srt = 45,               # Text label rotation angle
         col = COL2('PiYG'),
         addCoef.col = "black",     # Add correlation coefficients
         number.cex = 0.6,          # Size of correlation coefficients
         diag = FALSE)              # Remove diagonal


# Test the environmental effects on the taxonomic composition.
library(vegan)
# Determine the microbial dissimilarity matrix based on the bray-curties distance
tax_16s_hell_trans <- decostand(t(otu_16s), method = "hellinger")
tax_16s_dist <-vegdist(tax_16s_hell_trans, "bray")
tax_its_hell_trans <- decostand(t(otu_its), method = "hellinger")
tax_its_dist <-vegdist(tax_its_hell_trans, "bray")

fun_dist <-vegdist(t(ko_tpm_table), "bray" )

mag_hell_trans <- decostand(t(abundance_tab), method = "hellinger")
mag_dist <-vegdist(mag_hell_trans, "bray")
cwm_results_trans <- decostand(cwm_results_final, method = "log")
trait_dist <-vegdist(cwm_results_trans, "bray")



# Determine the environmental dissimilarity matrix based on the euclidean distance
sel.vars <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP") # Select the environmental variables without collinearity
dist.env.list <- list()
for (i in sel.vars) {
  df <- data.frame(metadata[,i])
  rownames(df) <- rownames(metadata)
  env_dist <- vegdist(scale(df), "euclidean", na.rm = T)
  dist.env.list[[i]] <- env_dist
}

# Determine the taxonomic difference between the collapsed and un-collapsed areas within each site based on the bray-curties distance
gp <- paste(metadata$Gully_id, metadata$Group, sep = "_")
extract_fun <- function(df) {
  simi_group_df <- usedist::dist_groups(df, gp) # transform the distance matrices into similarity matrices
  between_df <- simi_group_df %>% filter(Label %in% grep("^Between", simi_group_df$Label, value = T))
  simi_between_group_df <- rbind(between_df %>% filter(str_detect(Item1, 'G1')) %>% filter(str_detect(Item2, 'G1')), 
                                 between_df %>% filter(str_detect(Item1, 'G2')) %>% filter(str_detect(Item2, 'G2')), 
                                 between_df %>% filter(str_detect(Item1, 'G3')) %>% filter(str_detect(Item2, 'G3')), 
                                 between_df %>% filter(str_detect(Item1, 'G4')) %>% filter(str_detect(Item2, 'G4')), 
                                 between_df %>% filter(str_detect(Item1, 'G5')) %>% filter(str_detect(Item2, 'G5')), 
                                 between_df %>% filter(str_detect(Item1, 'G6')) %>% filter(str_detect(Item2, 'G6'))) %>%
    as_tibble() %>% 
    select(Label, Distance) %>%
  return(simi_between_group_df)
}

# Bacteria
dis_between_16s_df <- extract_fun(tax_16s_dist)
# Fungi
dis_between_its_df <- extract_fun(tax_its_dist)
# Function
dis_between_fun_df <- extract_fun(fun_dist)
# MAGs
dis_between_mag_df <- extract_fun(mag_dist)
# Traits
dis_between_trait_df <- extract_fun(trait_dist)

# For each environmental variable, and change the distance column name to selected variable name
dis_between_env_list <- list()
for (i in sel.vars){
  dis_between_env_list[[i]] <- extract_fun(dist.env.list[[i]])
  colnames(dis_between_env_list[[i]])[2] <- i
}
# Merge the all data frame using cbind according to the Label column
dis_between_env_table <- cbind(dis_between_env_list[[sel.vars[1]]], 
                               dis_between_env_list[[sel.vars[2]]][,2],
                               dis_between_env_list[[sel.vars[3]]][,2],
                               dis_between_env_list[[sel.vars[4]]][,2],
                               dis_between_env_list[[sel.vars[5]]][,2],
                               dis_between_env_list[[sel.vars[6]]][,2],
                               dis_between_env_list[[sel.vars[7]]][,2],
                               distance_16s = dis_between_16s_df$Distance,
                               distance_its = dis_between_its_df$Distance,
                               distance_fun = dis_between_fun_df$Distance,
                               distance_mag = dis_between_mag_df$Distance,
                               distance_trait = dis_between_trait_df$Distance) %>%
  as.data.frame() %>%
  mutate(Gully_id = case_when(grepl("EB_", Label) ~ "EB", # add the gully id
                           grepl("ML_", Label) ~ "ML",
                           grepl("RS_", Label) ~ "RS",
                           grepl("SLH_", Label) ~ "SLH",
                           grepl("HSX_", Label) ~ "HSX",
                           grepl("HH_", Label) ~ "HH"),
         Time = case_when(Gully_id == "EB" ~ 52, # add the collapsed time
                          Gully_id == "ML" ~ 39,
                          Gully_id == "RS" ~ 54,
                          Gully_id == "SLH" ~ 24,
                          Gully_id == "HSX" ~ 42,
                          Gully_id == "HH" ~ 25),
         Slope = case_when(Gully_id == "EB" ~ 2.2, # add the slopes of the gullies
                           Gully_id == "ML" ~ 4.8,
                           Gully_id == "RS" ~ 11,
                           Gully_id == "SLH" ~ 4.8,
                           Gully_id == "HSX" ~ 8.6,
                           Gully_id == "HH" ~ 3),
         MAP = case_when(Gully_id == "EB" ~ 367, # add the mean annual precipitation of the gullies
                           Gully_id == "ML" ~ 411,
                           Gully_id == "RS" ~ 493,
                           Gully_id == "SLH" ~ 402,
                           Gully_id == "HSX" ~ 353,
                           Gully_id == "HH" ~ 436))

# Load the required packages for modeling and plot
library(lme4)
library(lmerTest)
library(ggplot2)
library(tidyr)
library(dplyr)

# Define predictors
predictors <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP")

# Create models for both distance indices
create_models <- function(response_var) {
  models <- lapply(predictors, function(x) {
    lmer(reformulate(c(x, "Time", "Slope", "MAP", "(1 | Gully_id)"), response_var), 
         data = dis_between_env_table)
  })
  names(models) <- paste0("mod", 1:7)
  return(models)
}

# Get models for both distance indices
models_16s <- create_models("distance_16s")
models_its <- create_models("distance_its")
models_fun <- create_models("distance_fun")
models_mag <- create_models("distance_mag")
models_trait <- create_models("distance_trait")

# Load additional required package for R-squared calculation
library(MuMIn)  # For r.squaredGLMM function

# Function to extract comprehensive model results including R-squared
extract_lmm_results <- function(models, distance_type, predictors) {
  map2_dfr(models, predictors, ~{
    model <- .x
    predictor <- .y
    
    # Get model summary
    model_summary <- summary(model)
    
    # Extract fixed effects
    fixed_coeff <- model_summary$coefficients
    
    # Extract random effects variance
    random_var <- as.data.frame(VarCorr(model))
    gully_var <- random_var$vcov[random_var$grp == "Gully_id"]
    residual_var <- random_var$vcov[random_var$grp == "Residual"]
    
    # Calculate ICC (Intraclass Correlation Coefficient)
    icc <- gully_var / (gully_var + residual_var)
    
    # Calculate R-squared values using MuMIn package
    r_squared <- r.squaredGLMM(model)
    marginal_r2 <- r_squared[1, "R2m"]    # R²m: variance explained by fixed effects only
    conditional_r2 <- r_squared[1, "R2c"] # R²c: variance explained by fixed + random effects
    
    # Alternative manual calculation of marginal R-squared
    # Get fitted values with only fixed effects
    fixed_fitted <- predict(model, re.form = NA)
    observed <- model.response(model.frame(model))
    
    # Manual marginal R-squared calculation
    var_fixed <- var(fixed_fitted)
    var_observed <- var(observed)
    marginal_r2_manual <- var_fixed / var_observed
    
    # Manual conditional R-squared calculation
    fitted_full <- fitted(model)
    var_fitted_full <- var(fitted_full)
    conditional_r2_manual <- var_fitted_full / var_observed
    
    # Extract model fit statistics
    loglik <- logLik(model)
    aic <- AIC(model)
    bic <- BIC(model)
    
    # Get ANOVA results for the main predictor - INCLUDING DEGREES OF FREEDOM
    anova_result <- anova(model)
    f_value <- anova_result$`F value`[1]
    p_value <- anova_result$`Pr(>F)`[1]
    
    # Extract degrees of freedom
    num_df <- anova_result$NumDF[1]  # Numerator degrees of freedom
    den_df <- anova_result$DenDF[1]  # Denominator degrees of freedom
    
    # Extract coefficients for all fixed effects
    predictor_row <- which(rownames(fixed_coeff) == predictor)
    
    # Main predictor results
    if(length(predictor_row) > 0) {
      predictor_est <- fixed_coeff[predictor_row, "Estimate"]
      predictor_se <- fixed_coeff[predictor_row, "Std. Error"]
      predictor_t <- fixed_coeff[predictor_row, "t value"]
      predictor_p <- fixed_coeff[predictor_row, "Pr(>|t|)"]
      predictor_df <- fixed_coeff[predictor_row, "df"]  # t-test degrees of freedom
    } else {
      predictor_est <- predictor_se <- predictor_t <- predictor_p <- predictor_df <- NA
    }
    
    # Covariate results
    time_row <- which(rownames(fixed_coeff) == "Time")
    slope_row <- which(rownames(fixed_coeff) == "Slope") 
    map_row <- which(rownames(fixed_coeff) == "MAP")
    
    time_est <- if(length(time_row) > 0) fixed_coeff[time_row, "Estimate"] else NA
    time_p <- if(length(time_row) > 0) fixed_coeff[time_row, "Pr(>|t|)"] else NA
    
    slope_est <- if(length(slope_row) > 0) fixed_coeff[slope_row, "Estimate"] else NA
    slope_p <- if(length(slope_row) > 0) fixed_coeff[slope_row, "Pr(>|t|)"] else NA
    
    map_est <- if(length(map_row) > 0) fixed_coeff[map_row, "Estimate"] else NA
    map_p <- if(length(map_row) > 0) fixed_coeff[map_row, "Pr(>|t|)"] else NA
    
    # Calculate standardized coefficient for main predictor
    if(!is.na(predictor_est)) {
      pred_sd <- sd(dis_between_env_table[[predictor]], na.rm = TRUE)
      resp_sd <- sd(dis_between_env_table[[distance_type]], na.rm = TRUE)
      std_coeff <- predictor_est * (pred_sd / resp_sd)
    } else {
      std_coeff <- NA
    }
    
    tibble(
      Response = distance_type,
      Predictor = predictor,
      # Main predictor
      Estimate = predictor_est,
      SE = predictor_se,
      t_value = predictor_t,
      t_df = predictor_df,  # Added t-test degrees of freedom
      p_value = predictor_p,
      Std_Coeff = std_coeff,
      # ANOVA results with degrees of freedom
      F_value = f_value,
      NumDF = num_df,      # Added numerator degrees of freedom
      DenDF = den_df,      # Added denominator degrees of freedom
      F_p_value = p_value,
      # R-squared values
      Marginal_R2 = marginal_r2,      # R²m: fixed effects only
      Conditional_R2 = conditional_r2, # R²c: fixed + random effects
      Marginal_R2_manual = marginal_r2_manual,    # Manual calculation
      Conditional_R2_manual = conditional_r2_manual, # Manual calculation
      # Covariates
      Time_Est = time_est,
      Time_p = time_p,
      Slope_Est = slope_est,
      Slope_p = slope_p,
      MAP_Est = map_est,
      MAP_p = map_p,
      # Random effects
      Gully_Var = gully_var,
      Residual_Var = residual_var,
      ICC = icc,
      # Model fit
      LogLik = as.numeric(loglik),
      AIC = aic,
      BIC = bic,
      N_obs = nobs(model)
    )
  })
}

# Extract results for both distance indices
results_16s <- extract_lmm_results(models_16s, "distance_16s", predictors)
results_its <- extract_lmm_results(models_its, "distance_its", predictors)
results_fun <- extract_lmm_results(models_fun, "distance_fun", predictors)
results_mag <- extract_lmm_results(models_mag, "distance_mag", predictors)
results_trait <- extract_lmm_results(models_trait, "distance_trait", predictors)

all_results <- bind_rows(results_16s, results_its, results_fun, results_mag, results_trait)

# Add formatted labels and significance
all_results <- all_results %>%
  mutate(
    Response_Label = case_when(
      Response == "distance_16s" ~ "Bacteria",
      Response == "distance_its" ~ "Fungi",
      Response == "distance_fun" ~ "Functional gene",
      Response == "distance_mag" ~ "Metagenome assembled genomes",
      Response == "distance_trait" ~ "CWM traits"
    ),
    Predictor_Label = case_when(
      Predictor == "Plant_richness" ~ "Plant richness",
      Predictor == "AGB" ~ "Aboveground biomass",
      Predictor == "pH" ~ "Soil pH",
      Predictor == "SOC" ~ "Soil organic carbon",
      Predictor == "NH4_N" ~ "NH₄⁺-N",
      Predictor == "NO3_N" ~ "NO₃⁻-N", 
      Predictor == "AP" ~ "Available phosphorus"
    ),
    # Significance stars
    Sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    F_Sig = case_when(
      F_p_value < 0.001 ~ "***",
      F_p_value < 0.01 ~ "**",
      F_p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Format estimates with significance
    Est_Formatted = paste0(sprintf("%.3f", Estimate), Sig),
    Std_Coeff_Formatted = paste0(sprintf("%.3f", Std_Coeff), Sig),
    # Format F-statistic with degrees of freedom
    F_Formatted = paste0("F(", NumDF, ",", sprintf("%.1f", DenDF), ") = ", sprintf("%.2f", F_value))
  )

# Create comprehensive results table with degrees of freedom
results_table <- all_results %>%
  select(
    `Response Variable` = Response_Label,
    `Environmental Predictor` = Predictor_Label, 
    `n` = N_obs,
    `Estimate` = Estimate,
    `Standardized β` = Std_Coeff,
    `t-value` = t_value,
    `NumDF` = NumDF, # Added numerator degrees of freedom
    `DenDF` = DenDF,  # Added denominator degrees of freedom
    `F-statistic` = F_value,
    `r_squared` = `Marginal_R2`,
    `p-value` = p_value,
    `Significance` = Sig
  ) %>%
  mutate(
    `t-value` = round(`t-value`, 2),
    `DenDF` = round(`DenDF`, 1),
    `Standardized β` = round(`Standardized β`, 2),
    `F-statistic` = round(`F-statistic`, 2),
    `p-value` = case_when(
      `p-value` < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.3f", `p-value`)
    )
  )

# Display the results table
print(results_table)

# 1. Heatmap for all predictors
rep_str2 = c("Plant_richness" = "Plant richness",
             "AGB" = "AGB",
             "pH" = "pH",
             "SOC" = "SOC",
             "NH4_N" = expression(paste("NH"[4]^"+", "-N")),
             "NO3_N" = expression(paste("NO"[3]^"-", "-N")),
             "AP" = "AP"
)
rep_str3 = c("distance_16s" = "Bacteria",
             "distance_its" = "Fungi",
             "distance_fun" = "Functional gene",
             "distance_mag" ~ "Metagenome assembled genomes",
             "distance_trait" ~ "CWM traits")
# Prepare the heatmap data
all_coeff <- all_results %>%
  select(Response, Predictor, Std_Coeff, Marginal_R2, p_value, Sig) %>%
  mutate(Response = factor(Response, levels = c("distance_16s", "distance_its",
                                                "distance_fun", "distance_mag", 
                                                "distance_trait")),
         Predictor = factor(Predictor, levels = rev(predictors)))

heatmap_plot <- all_coeff %>%
  ggplot(aes(Predictor, Response, fill = Std_Coeff)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Sig), 
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#79ceb8", mid = "white", high = "#e95f5c", 
                       name = "Standardized\nCoefficient") +
  coord_flip() + 
  scale_x_discrete(labels = rep_str2) +
  scale_y_discrete(labels = rep_str3) +
  labs(
    x = "Changes in environmental variables",
    y = "Changes in microbial community composition"
  ) +
  theme_bw() +
  main_theme +
  theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'))

# 2. Regression plots for significant predictors only
# Filter significant predictors
sig_coeff <- filter(all_coeff, Predictor %in% c("pH", "SOC", "NO3_N", "AP"))

if (nrow(sig_coeff) > 0) {
  # Prepare data for significant predictors
  sig_data <- dis_between_env_table %>%
    select(all_of(c(as.vector(unique(sig_coeff$Predictor)), "distance_16s", 
                    "distance_its", "distance_fun",
                    "distance_mag", "distance_trait"))) %>%
    pivot_longer(c(distance_16s, distance_its, distance_fun, 
                   distance_mag, distance_trait), 
                 names_to = "Response", values_to = "distance_val") %>%
    pivot_longer(all_of(unique(sig_coeff$Predictor)), names_to = "Predictor",
                 values_to = "pred_val") %>%
    inner_join(sig_coeff, by = c("Response", "Predictor")) %>%
    mutate(Response = factor(Response, levels = c("distance_16s", "distance_its",
                                                  "distance_fun", "distance_mag", 
                                                  "distance_trait"))) %>%
    mutate(line_type = ifelse(p_value <= 0.05, "dashed", "solid")) %>%
    mutate(
      annotation = paste0("R² = ", round(Marginal_R2, 2), "\n",
                          if_else(p_value < 0.001, "p < 0.001", 
                                  paste("p =", round(p_value, 3))))
    )
}



# Create regression plot
regression_plot <- sig_data %>%
  mutate(Predictor = factor(Predictor, 
                            levels = c("pH", "SOC", "NO3_N", "AP"))) %>%
  ggplot(aes(pred_val, distance_val)) +
  geom_point(alpha = 0.6, size = 1.25, color = "#79ceb8") +
  geom_smooth(aes(linetype = line_type), method = "lm", 
              se = TRUE, color = "black", linewidth = 0.9) +
  ggpp::geom_text_npc(aes(npcx = 0.90, npcy = 0.20, label = annotation), 
                      hjust = 1.1, vjust = 1.1, size = 1.75, stat = "unique") +
  facet_grid(Response ~ Predictor, scales = "free") +
  labs(x = "Changes in environmental variables",
       y = "Changes in microbial community composition") +
  theme_bw() +
  main_theme +
  theme(legend.position = "none")

# Combine the two plots
library(cowplot)
combined_plot <-  ggdraw() +
  draw_plot(heatmap_plot, x = 0, y = 0, width = 0.35, height = 1) +
  draw_plot(regression_plot, x = 0.35, y = 0, width = 0.65, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 9,
                  x = c(0, 0.35), y = c(1, 1))

# Display the combined plot
combined_plot
# Save the combined plot
ggsave(file.path("E:/thermokarst_gully/revision/result/lmm_changes_env_microbes.pdf"), combined_plot,
       width = 8.8, height = 6.4)
