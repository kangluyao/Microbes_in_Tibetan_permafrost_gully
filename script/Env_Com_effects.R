# Data input 
## Set working and saving directory
setwd('e:/thermokarst_gully/')
save.dir <- file.path(getwd(),"result")

# Loading packages and read data
pacman::p_load(phyloseq, ape, vegan, Biostrings, microbiome, tidytable, tidyverse, rstatix)
source("script/read_data.R")

# Test the difference in the environmental variables between the uncollapsed and collapsed sites
env_vars <- c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt", "WHC", "SOC",
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
                "BGB" = expression(paste("BGB ", "(g/", "m"^2, ")")),
                "pH" = "pH",
                "Soil_moisture" = "Moisture (%)",
                "Clay_Silt" = "Clay+Silt (%)",
                "WHC" = "WHC (%)", 
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
  facet_wrap(~variable, scales = "free", ncol = 4, strip.position = "left", labeller = facet_labeller) +
  theme_classic() + 
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(colour = NA),
        strip.text = element_text(colour = 'black', size = 6, margin = margin()),
        axis.title = element_text(color = 'black',size = 6),
        axis.text = element_text(colour = 'black', size = 5),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        legend.title = element_text(colour = 'black', size = 6),
        legend.text = element_text(colour = 'black', size = 5),
        legend.key.size = unit(0.5, 'cm'))

lines <- tibble(variable = factor(c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt", "WHC", "SOC", "NH4_N", "NO3_N", "AP"), levels = env_vars),
                x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                xend = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                y = c(17, 270, 550, 8, 220, 72, 260, 185, 19, 40, 2.4),
                yend = c(17, 270, 550, 8, 220, 72, 260, 185, 19, 40, 2.4)
)
stars <- tibble(variable = factor(c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt", "WHC", "SOC", "NH4_N", "NO3_N", "AP"), levels = env_vars),
                x1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                y1 = c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95),
                sig.labels = env_result_lmm %>% filter(fixed_factors == "Group") %>%
                  select(sig) %>% pull,
                x2 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                y2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                panel.labels = letters[as.numeric(variable)]
)

env_plot <- env_plot +
  # geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE) + # add a sigment to denote the comparison between groups
  ggpp::geom_text_npc(data = stars, aes(npcx = x1, npcy = y1, label = sig.labels), inherit.aes = F) +
  ggpp::geom_text_npc(data = stars, aes(npcx = x2, npcy = y2, label = panel.labels), inherit.aes = F)

# save the plot
# ggsave(file.path(save.dir, "./figs/env_effect/env_comparison_plot.pdf"),
#        env_plot, width = 145, height = 90, units = "mm")
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
        strip.text = element_text(colour = 'black', size = 7),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 7),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.title = element_text(colour = 'black', size = 7),
        legend.text = element_text(colour = 'black', size = 6),
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

# Select the variables of low collinearity
sel_variables <- c("Gully_id",	"Group", "MAP",	"Time", "Slope",	
                   "Plant_richness", "AGB",	"BGB", "pH", "SOC", "NH4_N",	"NO3_N",	"AP")

# Construct a function for partial mantel test
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  #env.table <- env.table[complete.cases(env.table), ]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt", "WHC", "SOC", "NH4_N", "NO3_N", "AP")
  for (x in vars) {
    y.dist <- vegdist(scale(env.table[,x]), 'euclidean', na.rm = T)
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean', na.rm = T)
    mode <- mantel.partial(otu_table_hel_dist, y.dist, z.dist, 
                           method = "spearman", permutations = 999, na.rm = T)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}

# Extract the samples of each group
phylo_16s_uncollapsed <- subset_samples(phylo_16s, Group == 'Un-collapsed')
# Remove the ghost ASV with zero reads
phylo_16s_uncollapsed <- prune_taxa(taxa_sums(phylo_16s_uncollapsed)>=1, phylo_16s_uncollapsed)

phylo_its_uncollapsed <- subset_samples(phylo_its, Group == 'Un-collapsed')
phylo_its_uncollapsed <- prune_taxa(taxa_sums(phylo_its_uncollapsed)>=1, phylo_its_uncollapsed)

phylo_16s_collapsed <- subset_samples(phylo_16s, Group == 'Collapsed')
phylo_16s_collapsed <- prune_taxa(taxa_sums(phylo_16s_collapsed)>=1, phylo_16s_collapsed)

phylo_its_collapsed <- subset_samples(phylo_its, Group == 'Collapsed')
phylo_its_collapsed <- prune_taxa(taxa_sums(phylo_its_collapsed)>=1, phylo_its_collapsed)

# Run partial mantel tests
set.seed(123)
par.mant.16s.uncollapsed <- partial.mantel.fun(phylo_16s_uncollapsed)
set.seed(123)
par.mant.16s.collapsed <- partial.mantel.fun(phylo_16s_collapsed)

set.seed(123)
par.mant.its.uncollapsed <- partial.mantel.fun(phylo_its_uncollapsed)
set.seed(123)
par.mant.its.collapsed <- partial.mantel.fun(phylo_its_collapsed)

# Arrange the un-collapsed group table for plot
par.man.tibble.uncollapsed <- tibble(spec = c(rep('Bacteria', nrow(par.mant.16s.uncollapsed)),
                                              rep('Fungi', nrow(par.mant.its.uncollapsed))),
                                     rbind(par.mant.16s.uncollapsed, par.mant.its.uncollapsed))
par.man.tibble.uncollapsed

par.man.tibble.collapsed <- tibble(spec = c(rep('Bacteria', nrow(par.mant.16s.collapsed)),
                                              rep('Fungi', nrow(par.mant.its.collapsed))),
                                     rbind(par.mant.16s.collapsed, par.mant.its.collapsed))
par.man.tibble.collapsed

# Using linkET package to plot the associations between microbial community structure and environmental factors
library(linkET)
vars <- c("Plant_richness", "AGB", "BGB", "pH", "Soil_moisture", "Clay_Silt", "WHC", "SOC", "NH4_N", "NO3_N", "AP")
env.table <- metadata[ , vars]
mantel_uncollapsed <- par.man.tibble.uncollapsed %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
                 labels = c("<0.20", "0.20-0.4", ">0.40"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
set_corrplot_style()
p_uncollapsed <-  qcorrplot(correlate(env.table, engine = "Hmisc", method = "spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = p.value, size = r), 
              data = mantel_uncollapsed, 
              curvature = nice_curvature()) +
  scale_size_manual(values = c(0.3, 1, 1.75)) +
  scale_colour_manual(values = c("#e95f5c", "#79ceb8", '#3C5488FF', 'grey')) +
  # scale_x_discrete(labels = rep_str) +
  # scale_y_discrete(labels = rep_str) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = expression(paste("Spearman's ", rho), order = 3))) +
  theme(panel.grid=element_blank(), 
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"))

# ggsave(file.path(save.dir, './figs/env_effect/partail_matel_plot_for_uncollapsed.pdf'), p_uncollapsed, width = 5, height = 2.5, units = "in")
p_uncollapsed


# Next is the ploting for collapsed samples
mantel_collapsed <- par.man.tibble.collapsed %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
                 labels = c("<0.20", "0.20-0.4", ">0.40"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
set_corrplot_style()
p_collapsed <-  qcorrplot(correlate(env.table, engine = "Hmisc", method = "spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = p.value, size = r), 
              data = mantel_collapsed, 
              curvature = nice_curvature()) +
  scale_size_manual(values = c(0.3, 1, 1.75)) +
  scale_colour_manual(values = c("#e95f5c", "#79ceb8", '#3C5488FF', 'grey')) +
  # scale_x_discrete(labels = rep_str) +
  # scale_y_discrete(labels = rep_str) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = expression(paste("Spearman's ", rho), order = 3))) +
  theme(panel.grid=element_blank(), 
        axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5, hjust = 1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.8, "lines"))

# ggsave(file.path(save.dir, './figs/env_effect/partail_matel_plot_for_collapsed.pdf'), p_collapsed, width = 5, height = 2.5, units = "in")
p_collapsed

# Using the Null model to infer the ecological processes underlying the community assembly.
#Read in data
processimportance.file <- "E:/thermokarst_gully/data/null_model/ProcessImportance_EachGroup.txt"
iCAMP.compare.file <- "E:/thermokarst_gully/data/null_model/iCAMP.Compare.Group.txt"
processimportance.df <-  read.table(processimportance.file, header = TRUE, sep = "\t", 
                                    as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                                    check.names = FALSE)
iCAMP.compare.df <-  read.table(iCAMP.compare.file, header = TRUE, sep = "\t", 
                                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                                check.names = FALSE)


# plot
donut.df <- processimportance.df %>% 
  pivot_longer(cols = -c(Taxa, Group), names_to = "Process", 
               values_to = "Importance") %>%
  mutate(Taxa = factor(Taxa,  ordered = T,
                       levels = c("Bacteria", "Fungi")),
         Process = factor(Process, ordered = T,
                          levels = c('HeS', 'HoS', "DL", 'HD', 'DR')),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed"))) %>%
  mutate(Importance = round(Importance * 100, 1))


mycols <- c("#5cc3e8", "#ffdb00",  "#e95f5c", "#79ceb8", "#1F77B4")

# Plot
donut.bac.df <- donut.df %>% filter(Taxa == "Bacteria") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))
donut.fung.df <- donut.df %>% filter(Taxa == "Fungi") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))


donut_plot_fun <- function(df) {
  p <- ggplot(df) +
    # Outer donut, positioned further out with width for donut effect
    geom_bar(data = subset(df, Group == "Collapsed"), 
             aes(x = 3, y = Importance, fill = Process), 
             stat = "identity", color = "white", width = 0.5) +
    # Inner donut, positioned closer in with width for donut effect
    geom_bar(data = subset(df, Group == "Uncollapsed"), 
             aes(x = 2.4, y = Importance, fill = Process), 
             stat = "identity", color = "white", width = 0.5) +
    scale_fill_manual(values = mycols) +
    # Percentage labels for outer donut
    geom_text(data = subset(df, Group == "Collapsed"), 
              aes(x = 3, y = cumsum(Importance) - Importance / 2, 
                  label = paste0(Process,": ", round(Importance, 1), "%")), 
              size = 2, color = "black") +
    # Percentage labels for inner donut
    geom_text(data = subset(df, Group == "Uncollapsed"), 
              aes(x = 2.4, y = cumsum(Importance) - Importance / 2, 
                  label = paste0(Process,": ", round(Importance, 1), "%")), 
              size = 2, color = "black") +
    coord_polar(theta = "y") +
    # facet_wrap(~ Taxa, ncol = 1) +
    xlim(1, 3.5) +  # Adjust limits to create spacing and prevent clipping
    theme_void() +
    theme(legend.position = "left",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    labs(fill = NULL)
  return(p)
} 

donut_plot <- cowplot::plot_grid(donut_plot_fun(donut.bac.df), 
                                 donut_plot_fun(donut.fung.df), ncol = 2)
ggsave(file.path(save.dir, './figs/null_model/donut_plot1.pdf'), donut_plot, width = 3, height = 3)
donut_plot

######## cohen_d_plot
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}

bar_df <- iCAMP.compare.df %>% 
  mutate(Taxa = factor(Taxa,  ordered = T,
                       levels = c("Bacteria", "Fungi")),
         Process = factor(Process, ordered = T, 
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         sig = as.vector(unlist(lapply(P.value, p.stars))))


cohen_d_plot <- ggplot(data = bar_df, aes(x = Process, y = Cohen.d, fill = Process)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = mycols) +
  # scale_x_discrete(limits=rev) +
  scale_y_continuous(limits = c(-7.5, 7.5)) +
  geom_text(aes(label = sig, x = Process, y = Cohen.d * 1.1, vjust = 0.55)) +
  # labs(x = NULL, y = NULL) +
  facet_wrap(~ Taxa, ncol = 2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size = 6),
        strip.text = element_text(color='black',size = 6),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size = 6),
        axis.text.x = element_text(colour='black', size = 6),
        legend.title=element_text(size = 7),
        legend.text=element_text(size = 6),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
ggsave(file.path(save.dir, './figs/null_model/cohen_d_plot1.pdf'), cohen_d_plot, width = 3, height = 1.65)
iCAMP_plot <- cowplot::plot_grid(donut_plot, cohen_d_plot, ncol = 2)
iCAMP_plot























#Read in data
processimportance.file <- "E:/thermokarst_gully/data/null_model/ProcessImportance_EachGroup.txt"
iCAMP.compare.file <- "E:/thermokarst_gully/data/null_model/iCAMP.Compare.Group.txt"
processimportance.df <-  read.table(processimportance.file, header = TRUE, sep = "\t", 
                                    as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                                    check.names = FALSE)
iCAMP.compare.df <-  read.table(iCAMP.compare.file, header = TRUE, sep = "\t", 
                                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                                check.names = FALSE)


# plot
donut.df <- processimportance.df %>% 
  pivot_longer(cols = -c(Taxa, Group), names_to = "Process", 
               values_to = "Importance") %>%
  mutate(Taxa = factor(Taxa,  ordered = T,
                       levels = c("Bacteria", "Fungi", "Protist", "Animal")),
         Process = factor(Process, ordered = T,
                          levels = c('HeS', 'HoS', "DL", 'HD', 'DR')),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed"))) %>%
  mutate(Importance = round(Importance * 100, 1))


mycols <- c("#5cc3e8", "#ffdb00",  "#e95f5c", "#79ceb8", "#1F77B4")

# Plot
donut.bac.df <- donut.df %>% filter(Taxa == "Bacteria") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))
donut.fung.df <- donut.df %>% filter(Taxa == "Fungi") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))
donut.pro.df <- donut.df %>% filter(Taxa == "Protist") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))
donut.anim.df <- donut.df %>% filter(Taxa == "Animal") %>%
  mutate(Process = factor(Process, ordered = T,
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         Group = factor(Group, ordered = T,
                        levels = c("Uncollapsed", "Collapsed")))

donut_plot_fun <- function(df) {
  p <- ggplot(df) +
    # Outer donut, positioned further out with width for donut effect
    geom_bar(data = subset(df, Group == "Collapsed"), 
             aes(x = 3, y = Importance, fill = Process), 
             stat = "identity", color = "white", width = 0.5) +
    # Inner donut, positioned closer in with width for donut effect
    geom_bar(data = subset(df, Group == "Uncollapsed"), 
             aes(x = 2.4, y = Importance, fill = Process), 
             stat = "identity", color = "white", width = 0.5) +
    scale_fill_manual(values = mycols) +
    # Percentage labels for outer donut
    geom_text(data = subset(df, Group == "Collapsed"), 
              aes(x = 3, y = cumsum(Importance) - Importance / 2, 
                  label = paste0(Process,": ", round(Importance, 1), "%")), 
              size = 2, color = "black") +
    # Percentage labels for inner donut
    geom_text(data = subset(df, Group == "Uncollapsed"), 
              aes(x = 2.4, y = cumsum(Importance) - Importance / 2, 
                  label = paste0(Process,": ", round(Importance, 1), "%")), 
              size = 2, color = "black") +
    coord_polar(theta = "y") +
    # facet_wrap(~ Taxa, ncol = 1) +
    xlim(1, 3.5) +  # Adjust limits to create spacing and prevent clipping
    theme_void() +
    theme(legend.position = "left",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    labs(fill = NULL)
  return(p)
} 

donut_plot <- cowplot::plot_grid(donut_plot_fun(donut.bac.df), donut_plot_fun(donut.fung.df), 
                                 donut_plot_fun(donut.pro.df), donut_plot_fun(donut.anim.df), 
                                 ncol = 2)
# ggsave(file.path(save.dir, './figs/null_model/donut_plot.pdf'), donut_plot, width = 3, height = 3)
donut_plot
```
cohen_d_plot
```{r}
# all_donut_charts <- ggplot(donut.df, aes(x = 2, y = Importance, fill = Process)) +
#   geom_bar(position = 'fill', stat = 'identity') +
#   facet_grid(Taxa ~ Group, switch = "y") +
#   xlim(0.5, 2.5) +
#   scale_fill_manual(values = mycols) +
#   coord_polar(theta = 'y') +
#   labs(x = NULL, y = NULL) +
#   # geom_text(aes(y = lab.ypos, label = Importance), color = "black") +
#   theme_bw() +
#   theme(legend.title = element_text(size = 7),
#         legend.text = element_text(size = 6),
#         legend.position = "left",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank())


######## cohen_d_plot
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}

bar_df <- iCAMP.compare.df %>% 
  mutate(Taxa = factor(Taxa,  ordered = T,
                       levels = c("Bacteria", "Fungi", "Protist", "Animal")),
         Process = factor(Process, ordered = T, 
                          levels = rev(c('HeS', 'HoS', "DL", 'HD', 'DR'))),
         sig = as.vector(unlist(lapply(P.value, p.stars))))


cohen_d_plot <- ggplot(data = bar_df, aes(x = Process, y = Cohen.d, fill = Process)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = mycols) +
  # scale_x_discrete(limits=rev) +
  scale_y_continuous(limits = c(-7.5, 7.5)) +
  geom_text(aes(label = sig, x = Process, y = Cohen.d * 1.1, vjust = 0.55)) +
  # labs(x = NULL, y = NULL) +
  facet_wrap(~ Taxa, ncol = 2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size = 6),
        strip.text = element_text(color='black',size = 6),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size = 6),
        axis.text.x = element_text(colour='black', size = 6),
        legend.title=element_text(size = 7),
        legend.text=element_text(size = 6),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
# ggsave(file.path(save.dir, './figs/null_model/cohen_d_plot.pdf'), cohen_d_plot, width = 3, height = 3)
iCAMP_plot <- cowplot::plot_grid(donut_plot, cohen_d_plot, ncol = 2)
iCAMP_plot

