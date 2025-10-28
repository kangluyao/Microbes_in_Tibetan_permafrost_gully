data(dune)
data(dune.env)
(sim_phyla <- with(metadata, simper(t(singlem_phyla), Group, permutations = 999)))
## IGNORE_RDIFF_BEGIN
summary(sim_phyla)
## IGNORE_RDIFF_END
(sim_genu <- with(metadata, simper(t(singlem_genu), Group, permutations = 999)))
## IGNORE_RDIFF_BEGIN
summary(sim_genu)
## IGNORE_RDIFF_END

library(ggplot2)
library(dplyr)

sim_result <- summary(sim_genu)$`Un-collapsed_Collapsed`

simper_df <- sim_result %>% 
  rownames_to_column("taxonomy") %>%
  mutate(taxonomy = factor(taxonomy, levels = taxonomy))

genera_id <- simper_df %>% pull(taxonomy)

cutoff_pct <- 0.70
cutoff_index <- which(simper_df$cumsum >= cutoff_pct)[1]
cutoff_species <- simper_df$taxonomy[cutoff_index]
high_var_species <- c(
  simper_df$taxonomy[simper_df$cumsum < cutoff_pct], 
  cutoff_species
)

scale_factor <- max(simper_df$average)

simper_contri_plot <- ggplot(simper_df, aes(x = taxonomy)) +
  geom_bar(
    aes(y = average), 
    stat = "identity", 
    fill = "#00BFC4", 
    alpha = 0.7
  ) +
  geom_line(
    aes(y = cumsum * scale_factor), 
    group = 1, 
    color = "red", 
    linewidth = 1
  ) +
  geom_point(
    aes(y = cumsum * scale_factor), 
    color = "red", 
    size = 1
  ) +
  geom_hline(
    yintercept = cutoff_pct * scale_factor, 
    linetype = "dashed", 
    color = "red"
  ) +
  geom_vline(
    xintercept = cutoff_index, 
    linetype = "dotted", 
    color = "darkred", 
    linewidth = 1
  ) +
  scale_y_continuous(
    name = "Individual Contribution (%)",
    sec.axis = sec_axis(
      ~ . / scale_factor, 
      name = "Cumulative Contribution",
      labels = scales::percent
    )
  ) +
  labs(
    title = "SIMPER Species Contribution Analysis",
    subtitle = paste("Cutoff species (70% threshold):", cutoff_species)
  ) +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.ticks.x = element_blank()
  )


simper_contri_plot

# merge similar_all_data with meta_unique
singlem_genera_df <- data.frame(metadata[, c("Sample_name", "Gully_id", 
                                             "Group", "Time", "Slope", "MAP")]) %>%
  left_join(singlem_genu %>% t() %>%
              as.data.frame() %>% rownames_to_column("Sample_name"),
            by = "Sample_name")

genera_abun_stats <- singlem_genera_df %>% 
  get_summary_stats(genera_id, type = "common") %>% #or using type = "mean_sd"
  arrange(desc(mean))
genera_abun_stats


# Test the effects of permafrost collapsed on community structure of each specific group
library(lme4)
library(lmerTest)
singlem_genera_scale <- singlem_genera_df %>% 
  select(c("Group", "Gully_id", "Time", "Slope", "MAP", genera_id)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to 
# wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
result_list <- lapply(6:ncol(singlem_genera_scale), function(j) {
  message("Now j=", j, " in ", ncol(singlem_genera_scale), ". ", date())
  
  if (length(unique(singlem_genera_scale[, j])) < 3) {
    return(NULL)
  }
  
  fm1 <- lmer(singlem_genera_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), 
              data = singlem_genera_scale)
  
  presult <- car::Anova(fm1, type = 2)
  coef_summary <- coef(summary(fm1))
  
  coefs <- setNames(coef_summary[, "Estimate"], paste0(rownames(coef_summary), ".mean"))
  se_values <- setNames(coef_summary[, "Std. Error"], paste0(rownames(coef_summary), ".se"))
  t_values <- setNames(coef_summary[, "t value"], paste0(rownames(coef_summary), ".t"))
  chisq_p <- setNames(c(presult[, 1], presult[, 3]), 
                      c(paste0(rownames(presult), ".chisq"), paste0(rownames(presult), ".P")))
  
  data.frame(t(c(coefs, t_values, se_values, chisq_p)))
})

singlem_genera_S1 <- do.call(rbind, result_list[!sapply(result_list, is.null)])
rownames(singlem_genera_S1) <- colnames(singlem_genera_scale)[6:ncol(singlem_genera_scale)][!sapply(result_list, is.null)]

#plot
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_genu_comparison <- singlem_genera_S1 %>%
  rownames_to_column("variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = genera_id)) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), linewidth = 0.7,  colour = "gray2")+
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  # geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
  #           position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + scale_x_discrete(position = "bottom") +
  main_theme + theme(legend.position = "none",
                     axis.ticks.x = element_blank())

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_dis_comparison, width = 2.7, height = 5, units = "in")
single_genu_comparison

cowplot::plot_grid(simper_contri_plot, single_genu_comparison,
                   ncol = 1, align = "hv")



# Microbial traits
df_cor <- singlem_genera_S1 %>%
  rownames_to_column("variable") %>%
  select(variable, GroupCollapsed.mean, Group.P) %>%
  # mutate(GroupCollapsed.mean = abs(GroupCollapsed.mean)) %>%
  inner_join(sim_result %>%
               rownames_to_column("variable"),
             by = "variable") %>%
  rename(gtdb_taxonomy = variable)


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

gtdb_trait_genera_ave <- gtdb_trait %>%
  group_by(gtdb_taxonomy) %>%
  summarise(across(c(gc_percentage, genome_size, 
                     ncbi_rrna_count), mean, na.rm = TRUE)) %>%
  inner_join(df_cor, by = "gtdb_taxonomy") %>%
  mutate(high_contribution = if_else(gtdb_taxonomy %in% 
                                       high_var_species, "True", "False"),
         eff.group = if_else(GroupCollapsed.mean >= 0, "Positive", "Negative"))

head(gtdb_trait_genera_ave)


traits_id <- c("gc_percentage", "genome_size", "ncbi_rrna_count")
my_comparisons <- list(c('Negative', 'Positive'))
singlem_traits_comparison <- gtdb_trait_genera_ave %>% 
  filter(Group.P < 0.05) %>%
  # filter(high_contribution == "True") %>%
  select(c("eff.group", all_of(traits_id))) %>%
  filter(eff.group %in% c('Negative', 'Positive')) %>%
  gather(traits, value, -c("eff.group")) %>%
  mutate(eff.group = factor(eff.group, levels = c('Negative', 'Positive'))) %>%
  mutate(traits = factor(traits, levels = traits_id)) %>%
  ggplot(aes(eff.group, value, fill = eff.group)) +
  geom_half_violin(position = position_nudge(x = 0.25), 
                   side = "r", width = 0.6, color = NA, alpha = 0.75) +
  # geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = eff.group), size = 1.5, 
              width = 0.15, stroke = 0, pch = 21, alpha = 0.75) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.35, linewidth = 0.25) +
  stat_compare_means(comparisons = my_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", 
                     bracket.size = 0.5, bracket.width = 0.1,  size = 3.5,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~traits, scales = "free_y", ncol = 6) + #
  main_theme +
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))

singlem_traits_comparison






plot_correlation(
  data = gtdb_trait_genera_ave, # %>% 
    # filter(high_contribution == "True"),
  x_var = "gc_percentage",
  y_var = "average",
  method = "pearson")


plot_correlation(
  data = gtdb_trait_genera_ave %>% 
    filter(high_contribution == "True"),
  x_var = "genome_size",
  y_var = "average",
  method = "pearson")

plot_correlation(
  data = gtdb_trait_genera_ave %>% 
    filter(high_contribution == "True"),
  x_var = "ncbi_rrna_count",
  y_var = "average",
  method = "pearson")


plot_correlation(
  data = gtdb_trait_genera_ave %>% 
    filter( Group.P < 0.05 & GroupCollapsed.mean < 0 & ncbi_rrna_count < 10),
  x_var = "ncbi_rrna_count",
  y_var = "GroupCollapsed.mean",
  method = "spearman")

plot_correlation(
  data = df_cor %>% filter(Group.P < 0.05),
  x_var = "GroupCollapsed.mean",
  y_var = "average",
  method = "pearson")






gtdb_taxonomy %in% high_var_species &



mag_dist <-vegdist(mag_results_trans, "bray")

sim <- with(metadata, simper(t(mags_abun_tab), Group, permutations = 999))
sim_result <- summary(sim)
cutoff_index <- which(simper_df$cumsum >= 0.7)[1]

sim_contri_mags <- data.frame(sim_result$`Un-collapsed_Collapsed`) %>%
  # mutate(is_cutoff = cumsum >= 0.7) %>%
  # slice_head(n = which(.$is_cutoff)[1]) %>%
  rownames_to_column("Genome_id") %>%
  select(Genome_id, average) %>%
  inner_join(trait_data %>%
            rownames_to_column("Genome_id"),
            by = "Genome_id")

plot_correlation(
  data = sim_contri_mags,
  x_var = "Resource.Acquisition.Substrate.uptake.carbohydrate.transport",
  y_var = "average",
  method = "pearson" 
)
