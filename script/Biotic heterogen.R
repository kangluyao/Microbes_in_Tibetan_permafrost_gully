# Data input
# Set working and saving directory
setwd('e:/thermokarst_gully/') 
save.dir <- file.path(getwd(),"result")
# Load packages
pacman::p_load(phyloseq, ape, vegan, Biostrings, microbiome, tidyverse, rstatix)
# read data
source("script/read_data.R")
## Test the overall difference in the taxonomic and functional composition.
library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_16s_dist <-vegdist(t(otu_16s), "bray")
tax_its_dist <-vegdist(t(otu_its), "bray")
#permanova test the difference in compositional variance
adonis2(tax_16s_dist ~ Group, data = metadata)
adonis2(tax_its_dist ~ Group, data = metadata)


#Visualization for the overall difference by PCoA plot.
#16S
ord.16s <-  cmdscale(tax_16s_dist,  k = 2, eig = T, add = T)
pcoa_16s_plot <- data.frame(Group = metadata$Group, scores(ord.16s)) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
  ggplot(aes(x = Dim1, y = Dim2, shape = Group, color = Group)) + 
  geom_point(size = 1, alpha = 0.8) + 
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.2, 
               show.legend = FALSE, level = 0.95) +
  scale_colour_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  labs(x=paste("PCoA1 (", format(100 * ord.16s$eig[1] / sum(ord.16s$eig), digits = 3), "%)", sep = ""),
       y=paste("PCoA2 (", format(100 * ord.16s$eig[2] / sum(ord.16s$eig), digits = 3), "%)", sep = "")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
#ITS
ord.its <-  cmdscale(tax_its_dist,  k = 2, eig = T, add = T)
pcoa_its_plot <- data.frame(Group = metadata$Group, scores(ord.its)) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
  ggplot(aes(x = Dim1, y = Dim2, shape = Group, color = Group)) + 
  geom_point(size = 1, alpha = 0.8) + 
  stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.2, 
               show.legend = FALSE, level = 0.95) +
  scale_colour_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  labs(x=paste("PCoA1 (", format(100 * ord.its$eig[1] / sum(ord.its$eig), digits = 3), "%)", sep = ""),
       y=paste("PCoA2 (", format(100 * ord.its$eig[2] / sum(ord.its$eig), digits = 3), "%)", sep = "")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
pcoa_tax <- cowplot::plot_grid(pcoa_16s_plot, pcoa_its_plot, ncol = 2)
pcoa_tax

## Explore the difference in taxonomic variance between uncollapsed and collapsed soils
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
          'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_16s_data <- lapply(vars, function(x) 
  usedist::dist_subset(tax_16s_dist, 
                       grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_data <- data.frame(Group = similar_16s_data$Group,
                           Gully_id = similar_16s_data$Gully_id, 
                           distance_16s = similar_16s_data$distance,
                           distance_its = similar_its_data$distance)


# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_data with meta_unique
similar_df <- similar_data %>%
  left_join(meta_unique, by = "Gully_id")

## Linear mixed models test the effect of permafrost thawing on microbial diversity
library(lme4)
library(lmerTest)
dis_index <- c("distance_16s", "distance_its")
lmm_dis_modes <- lapply(dis_index, function(x) {
  lmer(substitute(i ~ Group + Time + Slope + MAP + (1 | Gully_id), list(i = as.name(x))), 
       data = similar_df)})
summary.model <- function(model){
  F.value <- anova(model)$'F value'
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
df <- NULL
for(i in 1:length(dis_index)) {
  tmp <- summary.model(lmm_dis_modes[[i]])
  if (is.null(df)){
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}

dis_result_lmm <-data.frame(dis_index = rep(dis_index, each = 4), group1 = rep("Un-collapsed", length(dis_index)),
                            variables = rep(c("Group", "Time", "Slope", "MAP"), 2),
                            group2 = rep("Collapsed", length(dis_index)), df)
dis_result_lmm


# Add the significant symbols manually
sig.dis.labs <- tibble(dis_index = factor(dis_index, levels = dis_index),
                       x1 = rep(0.5, length(dis_index)),
                       y1 = rep(0.95, length(dis_index)),
                       sig.labels = dis_result_lmm %>% filter(variables == "Group") %>%
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
dis_plot <- similar_df %>% select(c("Group", dis_index)) %>%
  gather(dis_index, value, -c("Group")) %>% 
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  mutate(dis_index = factor(dis_index, levels = c("distance_16s", "distance_its"))) %>%
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

# save the plot
# ggsave(file.path("E:/thermokarst_gully/result2/dis_groups_plot.pdf"),
#        dis_plot, width = 1.85, height = 3.75, units = "in")
dis_plot




## Identify the generalists, specialists, and opportunists
#determine the generalist and specialist species
library(EcolUtils)
spec_gen_16s <- spec.gen(phylo_16s@otu_table %>% t(), 
                         niche.width.method = "levins", perm.method = "quasiswap",
                         n = 1000, probs = c(0.05, 0.95))
# write.table(spec_gen_16s, file.path(save.dir, "/tables/spe_gen/spec_gen_16s_95.csv"), sep = ",")
spec_gen_its <- spec.gen(phylo_its@otu_table %>% t(), 
                         niche.width.method = "levins", perm.method = "quasiswap",
                         n = 1000, probs = c(0.05, 0.95))

## Explore the effect size of permafrost collapse on the microbial structure of generalists, specialists, and opportunists
##Bacteria
spec_gen_16s <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_16s_95.csv"), sep = ",", header = T)
spec_id_16s <- spec_gen_16s %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_16s <- spec_gen_16s %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_16s <- spec_gen_16s %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% spec_id_16s, ]
gen_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% gen_id_16s, ]
neu_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% neu_id_16s, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_16s_dist <-vegdist(t(spec_16s_comm), "bray")
tax_gen_16s_dist <-vegdist(t(gen_16s_comm), "bray")
tax_neu_16s_dist <-vegdist(t(neu_16s_comm), "bray")
# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
similar_spec_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_gen_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_neu_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))
#Combine all group datasets
similar_16s_data <- data.frame(Group = similar_spec_16s_data$Group,
                               Gully_id = similar_spec_16s_data$Gully_id, 
                               distance_16s_spec = similar_spec_16s_data$distance,
                               distance_16s_gen = similar_gen_16s_data$distance,
                               distance_16s_neu = similar_neu_16s_data$distance)

##Fungi##
spec_gen_its <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_its_95.csv"), sep = ",", header = T)
spec_id_its <- spec_gen_its %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_its <- spec_gen_its %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_its <- spec_gen_its %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% spec_id_its, ]
gen_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% gen_id_its, ]
neu_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% neu_id_its, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_its_dist <-vegdist(t(spec_its_comm), "bray")
tax_gen_its_dist <-vegdist(t(gen_its_comm), "bray")
tax_neu_its_dist <-vegdist(t(neu_its_comm), "bray")


# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
similar_spec_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_gen_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_neu_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))

similar_its_data <- data.frame(Group = similar_spec_its_data$Group,
                               Gully_id = similar_spec_its_data$Gully_id, 
                               distance_its_spec = similar_spec_its_data$distance,
                               distance_its_gen = similar_gen_its_data$distance,
                               distance_its_neu = similar_neu_its_data$distance)


# Combine all data table
similar_all_data <- data.frame(similar_data, similar_16s_data[, c(3:5)],
                               similar_its_data[, c(3:5)])

distance_index <- c("distance_16s_gen", "distance_16s_spec", "distance_16s_neu", 
                    "distance_its_gen", "distance_its_spec", "distance_its_neu")
distance_stats <- similar_all_data %>% group_by(Group) %>%
  get_summary_stats(distance_index, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  arrange(variable, Group)
distance_stats

# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_all_data with meta_unique
similar_all_df <- similar_all_data %>%
  left_join(meta_unique, by = "Gully_id")

# Test the effects of permafrost collapsed on community structure of each specific group
library(lme4)
library(lmerTest)
distance_index <- unique(distance_stats$variable)
dist_scale <- similar_all_df %>% 
  select(c("Group", "Gully_id", "Time", "Slope", "MAP", distance_index)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
dist_S1 <- sapply(6:ncol(dist_scale), function(j) {
  message("Now j=", j, " in ", ncol(dist_scale), ". ", date())
  if (length(unique(dist_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(dist_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = dist_scale)
    
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
colnames(dist_S1)<-colnames(dist_scale)[-c(1:5)]
data.frame(dist_S1)


#plot
rep_str = c("distance_16s_gen" = "Generalists",
            "distance_16s_spec" = "Specialists",
            "distance_16s_neu" = "Opportunists",
            "distance_its_gen" = "Generalists",
            "distance_its_spec" = "Specialists",
            "distance_its_neu" = "Opportunists"
)

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
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_dist_comparison <- dist_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(Taxa = factor(c(rep("Bacteria", 3), rep("Fungi", 3)), 
                       levels = c("Bacteria", "Fungi"))) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(distance_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), linewidth = 0.7,  colour = "gray2")+
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("grey", "#e95f5c")) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom", labels = rep_str) +
  main_theme + theme(legend.position = "none")

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_dis_comparison, width = 2.7, height = 5, units = "in")
single_dist_comparison


# Resistance analysis
gp <- paste(metadata$Gully_id, metadata$Group, sep = "_")
extract_fun <- function(df) {
  simi_group_df <- usedist::dist_groups(1-df, gp) # transform the distance matrices into similarity matrices
  between_df <- simi_group_df %>% filter(Label %in% grep("^Between", simi_group_df$Label, value = T))
  simi_between_group_df <- rbind(between_df %>% filter(str_detect(Item1, 'G1')) %>% filter(str_detect(Item2, 'G1')), 
                                 between_df %>% filter(str_detect(Item1, 'G2')) %>% filter(str_detect(Item2, 'G2')), 
                                 between_df %>% filter(str_detect(Item1, 'G3')) %>% filter(str_detect(Item2, 'G3')), 
                                 between_df %>% filter(str_detect(Item1, 'G4')) %>% filter(str_detect(Item2, 'G4')), 
                                 between_df %>% filter(str_detect(Item1, 'G5')) %>% filter(str_detect(Item2, 'G5')), 
                                 between_df %>% filter(str_detect(Item1, 'G6')) %>% filter(str_detect(Item2, 'G6'))) %>%
    as_tibble() %>% 
    select(Label, Distance) %>% rename(Resistance = Distance)
  return(simi_between_group_df)
}
              
#Bacteria
simi_between_16s.gen_df <- extract_fun(tax_gen_16s_dist)
simi_between_16s.spec_df <- extract_fun(tax_spec_16s_dist)
simi_between_16s.neu_df <- extract_fun(tax_neu_16s_dist)
#Fungi
simi_between_its.gen_df <- extract_fun(tax_gen_its_dist)
simi_between_its.spec_df <- extract_fun(tax_spec_its_dist)
simi_between_its.neu_df <- extract_fun(tax_neu_its_dist)
#merge the tables
resitance_table <- rbind(data.frame(Taxa = rep("Bacteria", 150), Group = rep("Generalists", 150),
                                 Gully_id = rep(c("EB", "HH", "HSX", "ML", "RS", "SLH"), each =25), 
                                 simi_between_16s.gen_df),
                      data.frame(Taxa = rep("Bacteria", 150), Group = rep("Specialists", 150),
                                 Gully_id = rep(c("EB", "HH", "HSX", "ML", "RS", "SLH"), each =25), 
                                 simi_between_16s.spec_df),
                      data.frame(Taxa = rep("Bacteria", 150), Group = rep("Opportunists", 150),
                                 Gully_id = rep(c("EB", "HH", "HSX", "ML", "RS", "SLH"), each =25), 
                                 simi_between_16s.neu_df),
                      data.frame(Taxa = rep("Fungi", 150), Group = rep("Generalists", 150),
                                 Gully_id = rep(c("EB", "HH", "HSX", "ML", "RS", "SLH"), each = 25), 
                                 simi_between_its.gen_df),
                      data.frame(Taxa = rep("Fungi", 150), Group = rep("Specialists", 150),
                                 Gully_id = c("EB", "HH", "HSX", "ML", "RS", "SLH"), 
                                 simi_between_its.spec_df),
                      data.frame(Taxa = rep("Fungi", 150), Group = rep("Opportunists", 150),
                                 Gully_id = rep(c("EB", "HH", "HSX", "ML", "RS", "SLH"), each =25), 
                                 simi_between_its.neu_df))

resitance_table

# Extract the unique Gully_id and corresponding Time, Slope, MAP from metadata
meta_unique <- metadata[, c("Gully_id", "Time", "Slope", "MAP")] %>%
  distinct(Gully_id, Time, Slope, MAP)

# merge similar_all_data with meta_unique
resitance_df <- resitance_table %>%
  left_join(meta_unique, by = "Gully_id")


# Test the difference in resistance among groups using linear mixed models
library(lme4)
library(lmerTest)
lmm_resistance_mode_bac <- lmer(Resistance ~ Group + Time + Slope + MAP + (1 | Gully_id), 
                                data = resitance_df %>% filter(Taxa == "Bacteria"))
lmm_resistance_mode_fungi <- lmer(Resistance ~ Group + Time + Slope + MAP + (1 | Gully_id), 
                                  data = resitance_df %>% filter(Taxa == "Fungi"))

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

summary.model(lmm_resistance_mode_bac)
summary.model(lmm_resistance_mode_fungi)

# Calculate the effect size
resistance_stats <- resitance_df %>% group_by(Taxa, Group) %>%
  get_summary_stats(Resistance, type = "common") %>% #or using type = "mean_sd"
  mutate(Taxa = factor(Taxa, levels = c("Bacteria", "Fungi"))) %>%
  mutate(Group = factor(Group, levels = rev(c("Generalists", "Specialists", "Opportunists"))))

# Create a box plot
# plot horizontal half box plots with resistance 
library(ggplot2)
library(gghalves)
resistance_comparison_plot <- resitance_df %>%
  mutate(Taxa = factor(Taxa, levels = c("Bacteria", "Fungi"))) %>%
  mutate(Group = factor(Group, levels = rev(c("Generalists", "Specialists", "Opportunists")))) %>%
  ggplot(aes(Group, Resistance, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.45, color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.2, size = 0.35, outlier.color = NA, alpha = 0.5,) +
  geom_jitter(aes(fill = Group, colour = Group), shape = 21, size = 0.25,
              width = 0.15, alpha = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  # ggpp::geom_text_npc(data = sig.dis.labs, aes(npcx = x1, npcy = y1, label = sig.labels), inherit.aes = F) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#5cc3e8", "#e95f5c", "#79ceb8", "#ffdb00")) +
  scale_color_manual(values = c("#5cc3e8", "#e95f5c", "#79ceb8", "#ffdb00")) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 2) +
  coord_flip()+
  main_theme +
  theme(legend.position = "none")
# ggsave(file.path("E:/thermokarst_gully/result2/resistance_plot.pdf"),
#        resistance_plot, width = 2.7, height = 3.5, units = "in")
resistance_comparison_plot


# Combine the plots
compo_venn_eff.size_resis_plot <- cowplot::plot_grid(single_dist_comparison, 
                                                     resistance_comparison_plot, 
                                                     ncol = 2, align = "h")
# ggsave(file.path("E:/thermokarst_gully/result2/spe_gen_comparison1.pdf"), 
#        compo_venn_eff.size_resis_plot, width = 7.5, height = 2)
compo_venn_eff.size_resis_plot




