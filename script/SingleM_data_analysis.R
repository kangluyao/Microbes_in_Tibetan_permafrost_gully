library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(stringr)
# Community composition
singlem_phyla <- fread("E:/thermokarst_gully/data/metagenome/singlem/gully-phylum.tsv") %>%
  filter(taxonomy != "unassigned") %>%
  mutate(
    taxonomy = str_extract(taxonomy, "p__[^;]+"),
  ) %>%
  mutate(across(c(taxonomy),  ~str_replace_all(., ".*__", ""))) %>%
  column_to_rownames("taxonomy")

singlem_phyla_df <- singlem_phyla %>%
  rownames_to_column("taxonomy") %>%
  pivot_longer(-taxonomy, names_to = "Sample", values_to = "Relative_abundance") %>%
  mutate(Group = case_when(grepl("_C", Sample) ~ "Un-collapsed",
                           grepl("_T", Sample) ~ "Collapsed"),
         Group = factor(Group, levels = c("Un-collapsed", "Collapsed")))

# 假设 df 包含列：Sample, Group, Phylum, Relative_abundance
# 计算每组平均相对丰度
plot_df <- singlem_phyla_df %>%
  group_by(Group, taxonomy) %>%
  summarise(mean_abund = mean(Relative_abundance), .groups = "drop")

# 排序 & 只保留前10类群
top_taxa <- plot_df %>% 
  group_by(taxonomy) %>%
  summarise(total = sum(mean_abund)) %>%
  # arrange(desc(total)) %>%
  slice_max(total, n = 10) %>%
  pull(taxonomy)

plot_df <- plot_df %>%
  mutate(taxonomy = ifelse(taxonomy %in% top_taxa, taxonomy, "Others")) %>%
  mutate(taxonomy = factor(taxonomy, levels = c(top_taxa, "Others")))

# 绘图
composition_plot <- ggplot(plot_df, aes(x = Group, y = mean_abund, fill = taxonomy)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_d3("category20") +
  labs(y = "Relative abundance (%)", x = NULL) +
  main_theme +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'))

## Diversity analysis
singlem_genu <- fread("E:/thermokarst_gully/data/metagenome/singlem/gully-genus.tsv") %>%
  filter(taxonomy != "unassigned") %>%
  mutate(
    taxonomy = str_extract(taxonomy, "g__[^;]+"),
  ) %>%
  mutate(across(c(taxonomy),  ~str_replace_all(., ".*__", ""))) %>%
  column_to_rownames("taxonomy")
singlem_genu = as.matrix(singlem_genu)

shannon  <- diversity(singlem_genu, index = "shannon")
richness <- specnumber(t(singlem_genu))

alpha_div <- data.frame(
  Richness = richness,
  Shannon  = shannon,
  Evenness = shannon / log(richness),
  stringsAsFactors = FALSE)

# Arrange the data for statistics
col_div_names <- c("Sample_name", 'Gully_id', 'Group', 
                   "Richness", "Shannon", "Evenness")
div_table <- cbind(metadata[, c("Sample_name", 'Gully_id', 'Group')], alpha_div) %>%
  rename_with(~ col_div_names) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

# Linear mixed models test the effect of collapsed and gully_id
div_index <- c("Richness", "Shannon", "Evenness")
div_scale <- div_table %>% 
  cbind(metadata[, c("Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", div_index))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
library(lme4)
library(lmerTest)
div_S1 <- sapply(6:ncol(div_scale), function(j) {
  if (length(unique(div_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(div_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = div_scale)
    
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
colnames(div_S1)<-colnames(div_scale)[-c(1:5)]
data.frame(div_S1)

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", "")))}

# Add the significant symbols manually
sig.div.labs <- tibble(div_index = factor(div_index, levels = div_index),
                       x1 = rep(0.5, length(div_index)),
                       y1 = rep(0.95, length(div_index)),
                       sig.labels = c("P = 0.003", "P < 0.001", "P < 0.001"),
                       panel.labels = letters[as.numeric(dis_index)])

# Create a plot
library(gghalves)
div_plot <- div_table %>% select(c("Group", div_index)) %>%
  gather(div_index, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(div_index = factor(div_index, 
                            levels = c("Richness", "Shannon", "Evenness"))) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
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
div_plot

#determine the dissimilarity matrix based on the bray-curties distance
library(vegan)
tax_singlem_dist <-vegdist(t(singlem_genu), "bray")
#permanova test the difference in compositional variance
adonis2(tax_singlem_dist ~ Group, data = metadata)

#Visualization for the overall difference by PCoA plot.
#singlem
pcoa_singlem_plot <- PCoA_plot_fun(tax_singlem_dist)

pcoa_singlem_plot

## Explore the difference in taxonomic variance between uncollapsed and collapsed soils
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
          'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_singlem_data <- lapply(vars, function(x) 
  usedist::dist_subset(tax_singlem_dist, 
                       grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_data <- data.frame(Group = similar_singlem_data$Group,
                           Gully_id = similar_singlem_data$Gully_id, 
                           distance_singlem = similar_singlem_data$distance)


# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
similar_df <- similar_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
library(lme4)
library(lmerTest)
lmm_dis_modes <- lmer(distance_singlem ~ Group + Time + 
                        Slope + MAP + (1 | Gully_id), data = similar_df)
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
summary.model(lmm_dis_modes )
# Create a box plot
library(gghalves)
dis_singlem_plot <- similar_df %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  ggplot(aes(Group, distance_singlem, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", 
                   width = 0.5, color = NA, alpha = 0.65) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.5,
              width = 0.15, alpha = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggpp::annotate("text_npc", npcx = 0.5, npcy = 0.95, 
                 size = 2, label = "P < 0.001") +
  labs(x = "Group", y = "Dissimilarity in microbial communities") +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(legend.position = "none")
dis_singlem_plot


# Combine all plots
library(cowplot)
combined_singlem_plot <-  ggdraw() +
  draw_plot(composition_plot, x = 0, y = 0, width = 2/5, height = 1) +
  draw_plot(div_plot, x = 2/5, y = 0.5, width = 3/5, height = 0.5) +
  draw_plot(pcoa_singlem_plot, x = 2/5, y = 0, width = 1.5/5, height = 0.5) +
  draw_plot(dis_singlem_plot, x = 3.5/5, y = 0, width = 1.5/5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 9,
                  x = c(0, 2/5, 2/5, 3.5/5), y = c(1, 1, 0.5, 0.5))

combined_singlem_plot

# save the plot
ggsave(file.path("E:/thermokarst_gully/revision/result/combined_singlem_plot.pdf"),
       combined_singlem_plot, width = 183, 
       height = 90, units = "mm", dpi = 600)
