# seting work directory
setwd('e:/thermokarst_gully')
wd_16s <- file.path(getwd(),"data/16S/rdp")
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")

# Loading packages
pacman::p_load(phyloseq, ape, vegan, Biostrings, 
               microbiome, tidytable, tidyverse, rstatix, networkD3)

# bin genome information
abundance_tab.file <- file.path(wd_fun, "MAGs/bin_abundance_coverm.txt")
tax.file <- file.path(wd_fun, "MAGs/annotation.txt")
mags.att.file <- file.path(wd_fun, "MAGs/MAGs_attributes.txt")

# Reading data
abundance_tab <- read.delim(abundance_tab.file, header = TRUE, sep = "\t", row.names = 1,
                            as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                            check.names = FALSE)[-1, ]
tax_bin <- read.delim(tax.file, header = TRUE, sep = "\t", as.is = TRUE, 
                      stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)

mags.att <- read.delim(mags.att.file, header = TRUE, sep = "\t", row.names = 1,
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE)
abundance_tab[1:5, 1:5]

mags_att <- abundance_tab %>% as.data.frame() %>%
  mutate(across(where(is.numeric), ~ as.numeric(. > 0))) %>% 
  mutate(rowsum = rowSums(.)) %>%
  select(rowsum) %>%
  rownames_to_column("Genome_id") %>%
  left_join(mags.att, by = "Genome_id")
  
seq_depth <- read.delim("E:/thermokarst_gully/revision/data/sequencing depth.csv", 
                        header = TRUE, sep = ",", as.is = TRUE, 
                        stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)


# Display the sequencing depth, mapping rate
library(tidyverse)
library(gghalves)

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
        legend.text = element_text(colour = 'black', size = 6),
        legend.key.size = unit(0.5, 'cm'))

seq_depth %>%
  pivot_longer(cols = -Sample, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c(
    "16S_clean_tags_per_sample", 
    "ITS_clean_tags_per_sample",
    "Depth (GB)", 
    "contig_mapping (%)",
    "MAG_mapping rate"
  ))) %>%
  group_by(variable) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE), # Calculate mean of mpg
    sd_value = sd(value, na.rm = TRUE),     # Calculate standard deviation of mpg
    n_value = n(),                        # Count the number of observations
    se_value = sd_value / sqrt(n_value),      # Calculate standard error
    .groups = "drop"                    # Ungroup the data frame after summarising
  )

seq_depth %>%
  pivot_longer(cols = -Sample, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c(
    "16S_clean_tags_per_sample", 
    "ITS_clean_tags_per_sample",
    "Depth (GB)", 
    "contig_mapping (%)",
    "MAG_mapping rate"
  ))) %>%
  ggplot(aes(x = 1, y = value)) + 
  geom_half_violin(
    fill = "#79ceb8", 
    position = position_nudge(x = 0.01), 
    side = "r", 
    width = 0.05, 
    color = NA
  ) + 
  geom_boxplot(
    fill = "#e95f5c", 
    width = 0.01, 
    size = 1.2, 
    outlier.color = NA, 
    position = position_nudge(x = 0.01)
  ) + 
  geom_jitter(
    fill = "#5cc3e8",
    shape = 21, 
    size = 3, 
    stroke = 0.5,
    width = 0.01, 
    alpha = 0.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  facet_wrap(. ~ variable, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Value") +
  theme_bw() +
  main_theme

# Explore the MAGs completeness and contamination distributions
library(gghalves)
library(ggpubr)

mags_att %>%
  summarise(
    mean_value = mean(Completeness, na.rm = TRUE), # Calculate mean of mpg
    sd_value = sd(Completeness, na.rm = TRUE),     # Calculate standard deviation of mpg
    n_value = n(),                        # Count the number of observations
    se_value = sd_value / sqrt(n_value),      # Calculate standard error
    .groups = "drop"                    # Ungroup the data frame after summarising
  )

mags_att %>%
  summarise(
    mean_value = mean(Contamination, na.rm = TRUE), # Calculate mean of mpg
    sd_value = sd(Contamination, na.rm = TRUE),     # Calculate standard deviation of mpg
    n_value = n(),                        # Count the number of observations
    se_value = sd_value / sqrt(n_value),      # Calculate standard error
    .groups = "drop"                    # Ungroup the data frame after summarising
  )


p1 <- ggplot(mags_att, aes(Completeness, Contamination)) + 
  geom_point(shape = 21, size = 3, alpha = 0.5) + 
  theme_bw() + 
  labs(x = "Completeness (%)", y = "Contamination (%)") +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.position = "none")
p1

library(gghalves)
p2 <- ggplot(mags_att, aes(1, Contamination)) + 
  geom_half_violin(fill = "#79ceb8", 
                   position = position_nudge(x = 0.26), 
                   side = "r", width = 0.5, color = NA) + 
  geom_boxplot(fill = "#e95f5c", width = 0.1, size = 1.2, 
               outlier.color = NA, position = position_nudge(x = 0.2)) + 
  geom_jitter(fill = "#5cc3e8",
              shape = 21, size = 3, width = 0.12, alpha = 0.5) + 
  theme_void() + 
  theme(legend.position = "none")
p2

p3 <- ggplot(mags_att, aes(1, Completeness)) + 
  geom_half_violin(fill = "#79ceb8", 
                   position = position_nudge(x = 0.26), 
                   side = "r", width = 0.6, color = NA, alpha = 0.5) + 
  geom_boxplot(fill = "#e95f5c", width = 0.1, size = 1.2, 
               outlier.color = NA, position = position_nudge(x = 0.2)) + 
  geom_jitter(fill = "#5cc3e8",
              shape = 21, size = 3, width = 0.12, alpha = 0.5) + 
  theme_void() + theme(legend.position = "none") + coord_flip()
p3


library(aplot) # Decorate a 'ggplot' with Associated Information
p1 %>% insert_top(p3, height = 0.4)%>%
  insert_right(p2, width = 0.4)
