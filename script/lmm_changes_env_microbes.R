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

mag_hell_trans <- decostand(t(mags_abun_tab), method = "hellinger")
mag_dist <-vegdist(mag_hell_trans, "bray")

cwm_results_trans <- decostand(cwm_results_final, method = "log")
trait_dist <-vegdist(cwm_results_trans, "bray")



# Determine the environmental dissimilarity matrix based on the euclidean distance
# Select the environmental variables without collinearity
sel.vars <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP")
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
dist.env.list <- list()
for (i in sel.vars) {
  df <- data.frame(metadata[,i])
  rownames(df) <- rownames(metadata)
  env_dist <- vegdist(scale(df), "euclidean", na.rm = T)
  dist.env.list[[i]] <- env_dist
}

dist.gene.list <- list()
for (i in gene_id_trait) {
  df <- data.frame(aggre_trait_group_data[,i])
  rownames(df) <- rownames(aggre_trait_group_data)
  gene_dist <- vegdist(scale(df), "euclidean", na.rm = T)
  dist.gene.list[[i]] <- gene_dist
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

dis_between_gene_list <- list()
for (i in gene_id_trait){
  dis_between_gene_list[[i]] <- extract_fun(dist.gene.list[[i]])
  colnames(dis_between_gene_list[[i]])[2] <- i
}

# Merge the all data frame using cbind according to the Label column
dis_between_all_table <- cbind(
  dis_between_env_list[[sel.vars[1]]],
  do.call(cbind, lapply(sel.vars[-1], function(v) 
    dis_between_env_list[[v]][, 2, drop = FALSE])),
  do.call(cbind, lapply(gene_id_trait, function(v) 
    dis_between_gene_list[[v]][, 2, drop = FALSE])),
  
  data.frame(
    distance_16s = dis_between_16s_df$Distance,
    distance_its = dis_between_its_df$Distance,
    distance_fun = dis_between_fun_df$Distance,
    distance_mag = dis_between_mag_df$Distance,
    distance_trait = dis_between_trait_df$Distance
  )
) %>%
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
# Define predictors and response variables
predictors <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP")
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
# Create models for both distance indices
create_models <- function(response_var) {
  models <- lapply(predictors, function(x) {
    lmer(reformulate(c(x, "Time", "Slope", "MAP", "(1 | Gully_id)"), response_var), 
         data = dis_between_env_table)
  })
  names(models) <- paste0("mod", 1:7)
  return(models)
}


library(purrr)
library(tibble)

# Create all conbination
response_vars <- c("distance_16s", "distance_its", "distance_fun", 
                   "distance_mag", "distance_trait", gene_id_trait)
model_combinations <- expand_grid(
  response = response_vars,
  predictor = predictors
)

# Build the models
model_results <- model_combinations %>%
  mutate(
    formula = map2(predictor, response, ~reformulate(
      termlabels = c(.x, "Time", "Slope", "MAP", "(1 | Gully_id)"),
      response = .y
    )),
    model = map(formula, ~lmer(.x, data = dis_between_all_table))
  )

# Results
model_results

# 
model_results %>%
  filter(response == "aromatic.acid.transport", 
         predictor == "Plant_richness") %>%
  pull(model) %>%
  .[[1]]


# Get models for both distance indices
models_16s <- create_models("distance_16s")
models_its <- create_models("distance_its")
models_fun <- create_models("distance_fun")
models_mag <- create_models("distance_mag")
models_trait <- create_models("distance_trait")

# Load additional required package for R-squared calculation
extract_model_stats <- function(model, response, predictor, data) {
  if(is.null(model)) {
    return(tibble(Response = response, Predictor = predictor, 
                  across(everything(), ~NA)))
  }
  
  tryCatch({
    # basic stastics
    coef_table <- summary(model)$coefficients
    pred_idx <- which(rownames(coef_table) == predictor)
    
    # fixed effects
    if(length(pred_idx) > 0) {
      est <- coef_table[pred_idx, "Estimate"]
      se <- coef_table[pred_idx, "Std. Error"]
      t_val <- coef_table[pred_idx, "t value"]
      p_val <- coef_table[pred_idx, "Pr(>|t|)"]
      t_df <- coef_table[pred_idx, "df"]
      std_coef <- est * sd(data[[predictor]], na.rm = TRUE) / 
        sd(data[[response]], na.rm = TRUE)
    } else {
      est <- se <- t_val <- p_val <- t_df <- std_coef <- NA
    }
    
    # ANOVA
    anova_res <- anova(model)
    f_val <- anova_res$`F value`[1]
    f_p <- anova_res$`Pr(>F)`[1]
    num_df <- anova_res$NumDF[1]
    den_df <- anova_res$DenDF[1]
    
    # R²
    r2 <- r.squaredGLMM(model)
    
    # covariable
    get_coef <- function(var) {
      idx <- which(rownames(coef_table) == var)
      if(length(idx) > 0) c(coef_table[idx, "Estimate"], coef_table[idx, "Pr(>|t|)"])
      else c(NA, NA)
    }
    time <- get_coef("Time")
    slope <- get_coef("Slope")
    map <- get_coef("MAP")
    
    # random effects
    var_df <- as.data.frame(VarCorr(model))
    gully_var <- var_df$vcov[var_df$grp == "Gully_id"]
    resid_var <- var_df$vcov[var_df$grp == "Residual"]
    
    tibble(
      Response = response,
      Predictor = predictor,
      Estimate = est,
      SE = se,
      t_value = t_val,
      t_df = t_df,
      p_value = p_val,
      Std_Coeff = std_coef,
      F_value = f_val,
      NumDF = num_df,
      DenDF = den_df,
      F_p_value = f_p,
      Marginal_R2 = r2[1, "R2m"],
      Conditional_R2 = r2[1, "R2c"],
      Time_Est = time[1],
      Time_p = time[2],
      Slope_Est = slope[1],
      Slope_p = slope[2],
      MAP_Est = map[1],
      MAP_p = map[2],
      Gully_Var = gully_var,
      Residual_Var = resid_var,
      ICC = gully_var / (gully_var + resid_var),
      LogLik = as.numeric(logLik(model)),
      AIC = AIC(model),
      BIC = BIC(model),
      N_obs = nobs(model)
    )
  }, error = function(e) {
    tibble(Response = response, Predictor = predictor, 
           across(everything(), ~NA))
  })
}

# Extract the model results
all_model_stats <- model_results %>%
  mutate(stats = pmap(list(model, response, predictor), 
                      ~extract_model_stats(..1, ..2, ..3, dis_between_all_table))) %>%
  unnest(stats)


# 1.Heatmap
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
all_model_stats %>%
  mutate(Sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""),
    Predictor = factor(Predictor, levels = c("Plant_richness", "AGB", "pH",
                                             "SOC", "NH4_N", "NO3_N", "AP")),
    Response = factor(Response, levels = c(unique(all_model_stats$Response)))) %>% 
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
  theme(strip.background = element_rect(fill = c("#FFF6E1")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.key.size = unit(0.25, 'cm'), #change legend key size
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
# ggsave(file.path("E:/thermokarst_gully/revision/result/lmm_changes_env_microbes.pdf"), combined_plot,
#        width = 8.8, height = 6.4)

# Explore the environmental effects on the C-S-R gene abundance and CWM traits
sel.vars <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP")
aggre_trait_group_env <- aggre_trait_group_data %>%
  rownames_to_column("Sample_name") %>%
  left_join(metadata[, c("Sample_name", "Gully_id", sel.vars)], by = "Sample_name")

#linner mixed model for matrix to matrix
lmm.mat.cal <- function(y, x, data, method){
  y <- as.matrix(y)
  x <- as.matrix(x)
  df<-NULL
  for(i in colnames(y)){
    for(j in colnames(x)){
      a <- y[, i, drop = F]
      b <- x[, j, drop = F]
      mode <- lmerTest::lmer(a ~ b + (1 | Gully_id), data, na.action=na.omit)
      coeff <- summary(mode)$coefficients[2,1]
      r.square <- round(MuMIn::r.squaredGLMM(mode)[1], 2)
      AIC <- round(AIC(mode), 2)
      p.value <- round(anova(mode)$Pr[1], 3)
      if (coeff>0) r = sqrt(r.square)
      else r = (-1) * sqrt(r.square)
      tmp <- c(i, j, r, r.square, AIC, p.value)
      if(is.null(df)){
        df <- tmp  
      }
      else{
        df <- rbind(df, tmp)
      }    
    }
  }
  df<-data.frame(row.names=NULL,df)
  colnames(df)<-c("dependent.variables","predictor.variables","Correlation","r.square", "AIC", "Pvalue")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$AdjPvalue<-rep(0,dim(df)[1])
  df$Correlation<-as.numeric(as.character(df$Correlation))
  #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
  # 1 -> donot adjust
  # 2 -> adjust predictor.variables + Type (column on the correlation plot)
  # 3 -> adjust dependent.variables + Type (row on the correlation plot for each type)
  # 4 -> adjust dependent.variables (row on the correlation plot)
  # 5 -> adjust predictor.variables (panel on the correlation plot)
  adjustment_label<-c("NoAdj","Adjpredictor.variablesAndType","Adjdependent.variablesAndType","Adjdependent.variables","Adjpredictor.variables")
  adjustment<-5
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$predictor.variables)){
      for(j in unique(df$Type)){
        sel<-df$predictor.variables==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$dependent.variables)){
      for(j in unique(df$Type)){
        sel<-df$dependent.variables==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$dependent.variables)){
      sel<-df$dependent.variables==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  } else if (adjustment==5){
    for(i in unique(df$predictor.variables)){
      sel<-df$predictor.variables==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  df$dependent.variables <-factor(df$dependent.variables, ordered = T, levels = rev(colnames(y)))
  df$predictor.variables <-factor(df$predictor.variables, ordered = T, levels = colnames(x))
  return(df)
}

# Define predictors
predictors <- c("Plant_richness", "AGB", "pH", "SOC", "NH4_N", "NO3_N", "AP")
gene_id_trait

x <- aggre_trait_group_env[, predictors]
y <- aggre_trait_group_env[, gene_id_trait]
corr_results <- lmm.mat.cal(y, x, aggre_trait_group_env, "spearman")

heatmap_lmm_plot <- corr_results %>%
  ggplot(aes(predictor.variables, dependent.variables, fill = round(Correlation, 2))) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Significance), 
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#79ceb8", mid = "white", high = "#e95f5c", 
                       name = "Standardized\nCoefficient") +
  coord_flip() + 
  scale_x_discrete(labels = rep_str2) +
  # scale_y_discrete(labels = rep_str3) +
  labs(
    x = "Changes in environmental variables",
    y = "Changes in microbial community composition"
  ) +
  main_theme +
  theme(strip.background = element_rect(fill = c("#FFF6E1")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.key.size = unit(0.25, 'cm'), #change legend key size
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'))
