# Data input
# Set working and saving directory
setwd('e:/thermokarst_gully/') 
save.dir <- file.path(getwd(),"result")
# Load packages
pacman::p_load(phyloseq, ape, vegan, Biostrings, 
               microbiome, tidytable, tidyverse, rstatix, networkD3)
# read data
source("script/read_data.R")

###### Composition analysis #######
# extract the taxa at phylum level for bacteria
subphylo_16s <- tax_glom(phylo_16s, 'Phylum', NArm = F)
subphylo_16s.rel  = transform_sample_counts(subphylo_16s, function(x) 100*x / sum(x))
ntaxa(subphylo_16s.rel)
ra.tab_16s <- otu_table(subphylo_16s.rel)
sum(ra.tab_16s[, 1])
subtaxa_tab_16s <- tax_table(subphylo_16s.rel)[, 2]

# extract the taxa at class level
subphylo_class_16s <- tax_glom(phylo_16s, 'Class', NArm = F)
subphylo.class.rel_16s  = transform_sample_counts(subphylo_class_16s, function(x) 100*x / sum(x))
ntaxa(subphylo.class.rel_16s)
ra.class.tab_16s <- otu_table(subphylo.class.rel_16s)
sum(ra.class.tab_16s[, 1])
subtaxa_class_tab_16s <- tax_table(subphylo.class.rel_16s)[, 3]

#ecombine the class within proteobacteria and other phyla into one abundance table
otu_final_tab_16s <- rbind(data.frame(subtaxa_class_tab_16s, ra.class.tab_16s) %>%
                             filter(Class %in% c("Alphaproteobacteria", "Betaproteobacteria",
                                                 "Gammaproteobacteria", "Deltaproteobacteria")) %>%
                             rename(Phylum = Class), data.frame(subtaxa_tab_16s, ra.tab_16s) %>%
                             filter(!Phylum %in% c("Proteobacteria", "Unassigned")))
otu_final_tab_16s[1:5, 1:5]
# extrac the phyla with the relative abundance higher than 0.1%
otu_rel_abun_0.1perc_16s <- otu_final_tab_16s %>%
  filter(!Phylum %in% c("Unassigned")) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab_16s)))) %>%
  filter(MRA >= 0.1) %>%
  arrange(desc(MRA)) %>%
  select(., -c('MRA'))
# write.table(otu_rel_abun_1perc_16s, file.path(save.dir, 'tables/otu_rel_abun_1perc_16s.csv'),
#             sep=',',  col.names = T, row.names = F, quote = FALSE)
names_16s <- otu_rel_abun_0.1perc_16s$Phylum
names_16s

## ITS composition
# extract the taxa at phylum level
subphylo_its <- tax_glom(phylo_its, 'Phylum', NArm = F)
subphylo_its.rel  = transform_sample_counts(subphylo_its, function(x) 100*x / sum(x))
ntaxa(subphylo_its.rel)
ra.tab_its <- otu_table(subphylo_its.rel)
sum(ra.tab_its[, 1])
subtaxa_tab_its <- tax_table(subphylo_its.rel)[, 2]

#ecombine the class within proteobacteria and other phyla into one abundance table
otu_final_tab_its <- data.frame(subtaxa_tab_its, ra.tab_its)

# extrac the phyla with the relative abundance higher than 0.1%
otu_rel_abun_0.1perc_its <- otu_final_tab_its %>%
  filter(!Phylum %in% c("Unassigned")) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab_its)))) %>%
  filter(MRA >= 0.1) %>%
  arrange(desc(MRA)) %>%
  select(., -c('MRA'))
# write.table(otu_rel_abun_1perc_its, file.path(save.dir, 'tables/otu_rel_abun_1perc_its.csv'),
#             sep=',',  col.names = T, row.names = F, quote = FALSE)
names_its <- otu_rel_abun_0.1perc_its$Phylum
names_its

# Combine the bacterial and fungal taxa table
compo_table <- rbind(otu_rel_abun_0.1perc_16s %>%
                    pivot_longer(-Phylum, names_to = "Sample_id", values_to = "Rela_abun") %>%
                    mutate(Group = rep("Bacteria", nrow(.))),
                  otu_rel_abun_0.1perc_its %>%
                    pivot_longer(-Phylum, names_to = "Sample_id", values_to = "Rela_abun") %>%
                    mutate(Group = rep("Fungi", nrow(.)))) %>%
  dplyr::rename(c(Taxa = Phylum)) %>%
  mutate(Type = ifelse(grepl("_C", Sample_id), "Un-collapsed",
                       ifelse(grepl("_T", Sample_id), "Collapsed", ""))) %>%
  select(-Sample_id) %>%
  summarise(across(everything(), mean), .by = c(Type, Group, Taxa)) %>%
  as.data.frame()

compo_df <- compo_table %>% mutate(group = case_when(Type == "Un-collapsed" ~ "a",
                                                     Type == "Collapsed" ~ "b"))
compo_df[1:5, ]


# prepare the table for sankey plot
sankey <- data.table::rbindlist(list(compo_df[c("Type", "Group", "Rela_abun", "group")],
                                     compo_df[c("Group", "Taxa", "Rela_abun", "group")]), 
                                use.names = FALSE)

names(sankey) <- c('source', 'target', 'value', "group")
knitr::kable(sankey[1:5, ])

# Make a connection data frame
links <- sankey 
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name = c("Un-collapsed", "Collapsed", "Bacteria", "Fungi", 
           otu_rel_abun_0.1perc_16s$Phylum,
           otu_rel_abun_0.1perc_its$Phylum)
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Add a 'group' column to each node. Here I decide to put all of them in the same group to make them grey
nodes$group <- as.factor(nodes$name)

# Ensure the links and nodes are correctly formatted
links_datset$IDsource <- as.integer(links_datset$IDsource)
links_datset$IDtarget <- as.integer(links_datset$IDtarget)
links_datset$value <- as.numeric(links_datset$value)

# Color palette
my_colors <- c("#e95f5c", "#79ceb8", "#1F77B4", "#AEC7E8", "#FF7F0E", 
               "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", 
               "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", 
               "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", 
               "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", 
               "#9EDAE5", "#1F77B4", "#AEC7E8")

# Create color JS string
colorJS <- paste0('d3.scaleOrdinal(["', 
                  paste(my_colors, collapse = '", "'), 
                  '"])')

# Create Sankey plot
sankey_plot <- sankeyNetwork(
  Links = links_datset, 
  Nodes = nodes, 
  Source = "IDsource", 
  Target = "IDtarget", 
  Value = "value", 
  NodeID = "name",
  NodeGroup = "group",
  LinkGroup = "group",
  colourScale = colorJS,
  iterations = 0, 
  fontFamily = 'Arial', 
  fontSize = 10, 
  nodeWidth = 10, 
  nodePadding = 15, 
  height = 475, 
  width = 250,
  sinksRight = FALSE
)

# Display or save the plot
sankey_plot

## 
### determine the average relative abundance of each phylum within each group
otu_rel_abun_0.1perc_all <- rbind(otu_rel_abun_0.1perc_16s, otu_rel_abun_0.1perc_its) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Phylum") %>%
  t() %>% data.frame() %>%
  mutate(Gully_id = metadata$Gully_id, Group = metadata$Group, .before = 1)

# determine the mean, standard deviations, standard errors for each phylum
avg_abun_stastics <- otu_rel_abun_0.1perc_all %>% select(-c("Gully_id")) %>%
  group_by(Group) %>% 
  summarise(across(where(is.numeric), list(mean = ~mean(., na.rm = T), 
                                           sd = ~sd(., na.rm = T), 
                                           n = ~n()))) %>% 
  pivot_longer(cols = -c(Group), 
               names_to = c('Phylum', '.value'), 
               names_sep = '_') %>% 
  mutate(se = sd/sqrt(n)) %>%
  arrange(factor(Phylum, levels = c(names_16s, names_its)))

avg_abun_stastics


# Test the difference in abundance of main taxa using liner models
# First, we need to test the collinearity of the data
# Determine the correlation matrix
cor_matrix <- cor(metadata[, -c(1:4)])
print(cor_matrix)
# visualize the correlation matrix
library(corrplot)
corrplot(cor_matrix, method = 'circle', type = 'lower', insig='blank', 
         tl.cex = 0.7, tl.col = "black", order = "original", col = COL2('PiYG'),
         addCoef.col ='black', number.cex = 0.6, diag=FALSE)

# Select the variables of low collinearity
sel_variables <- c("Gully_id",	"Group", "MAP",	"Time", "Slope",	
                   "Plant_richness", "AGB",	"BGB", "pH", "SOC", "NH4_N",	"NO3_N",	"AP")

# Loading the lme4 package
library(lme4)
# determine the effect size of the permafrost collapse on each phylum relative abundance
# In this analysis, we set the permafrost collapse as fixed factor, and collapsed time, gully slope, and MAP as covariates, 
# while gully_id as random effect.
otu_rel_abun_0.1perc_scale <- metadata[, c("MAP", "Slope", "Time")] %>% 
  cbind(., otu_rel_abun_0.1perc_all) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed")))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
abun_S1 <- sapply(6:ncol(otu_rel_abun_0.1perc_scale), function(j) {
  message("Now j=", j, " in ", ncol(otu_rel_abun_0.1perc_scale), ". ", date())
  if (length(unique(otu_rel_abun_0.1perc_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(otu_rel_abun_0.1perc_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), 
                data = otu_rel_abun_0.1perc_scale)
    
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
colnames(abun_S1) <- colnames(otu_rel_abun_0.1perc_scale)[-c(1:5)]
data.frame(abun_S1)

### effect size plot
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


all_tax_abun_comparison <- abun_S1 %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(c(names_16s, names_its)))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0.5, xmax = 8.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 8.5, xmax = 22.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path("E:/thermokarst_gully/result2/all_tax_abun_comparison.pdf"),
#        all_tax_abun_comparison, width = 2.5, height = 4, units = "in")
all_tax_abun_comparison

# Unique otus profile among three layers
# library
library(ggvenn)
# Make the venn plot for 16S
x_16s = list(
  Uncollapsed = otu_table(phylo_16s) %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_C', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames(),
  Collapsed = otu_table(phylo_16s) %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_T', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames()
)

p1 <- ggvenn(
  x_16s, 
  fill_color = c("#79ceb8", "#e95f5c"),
  stroke_size = ,
  stroke_color = NA,
  set_name_size = 3,
  text_size = 2.5,
  show_percentage = T,
  auto_scale = F
)
# Make the venn plot for ITS
x_its = list(
  Uncollapsed = otu_table(phylo_its) %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_C', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames(),
  Collapsed = otu_table(phylo_its) %>% data.frame() %>%
    mutate(rowsum = rowSums(select(., grep('_T', metadata$Sample_name, value = T)))) %>%
    filter(rowsum > 0) %>%
    rownames()
)

p2 <- ggvenn(
  x_its, 
  fill_color = c("#79ceb8", "#e95f5c"),
  stroke_color = NA,
  set_name_size = 3,
  text_size = 2.5,
  show_percentage = T,
  auto_scale = F
)
library(cowplot)
venn_plots <- plot_grid(p1, p2, ncol = 2, align = "v")
venn_plots
