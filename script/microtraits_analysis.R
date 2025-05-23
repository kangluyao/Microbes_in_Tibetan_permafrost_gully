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

# data preparation
mag_num_in_phylum <- data.frame(table(tax_bin[,"Phylum"]))
mag_num_in_phylum <- mag_num_in_phylum %>% arrange(desc(Freq))
mag_num_in_phylum <- rbind(mag_num_in_phylum[1:11, ], data.frame(Var1 = c('Others'),
                                                                 Freq = sum(mag_num_in_phylum[-c(1:11), 2])))

mag_num_in_phylum <- data.frame(Phylum = mag_num_in_phylum$Var1, mag_num = mag_num_in_phylum$Freq ,
                                prop = mag_num_in_phylum$Freq/sum(mag_num_in_phylum$Freq)*100)
mag_count.data <- mag_num_in_phylum %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
mag_count.data$Phylum <- factor(mag_count.data$Phylum, ordered = T, levels = mag_num_in_phylum$Phylum)

# pie plot 
# Define the colors you want
mycols <- c("#89c5da", "#ffc15c", "#74d944", "#CE50CA", "#5e738f", "#C0717C", "#CBD5ec", "#5F7FC7", 
            "#00718b", "#00b0f6", "#a3a500", "#f8766d", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")

pie_for_mag_num_phylum <- ggplot(mag_count.data, aes(x = "", y = prop, 
                                                     fill = reorder(Phylum, -lab.ypos))) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start = 0) +
  geom_text(aes(x = 1.35, y = lab.ypos, 
                label = paste0(mag_num, ' (', round(prop, 1), '%', ')', sep = '')),
            color = "black", size = 3) +
  scale_fill_manual('Phylum', values = mycols) +
  guides(fill = guide_legend(reverse = T)) +
  theme_void() +
  theme(legend.position = "left",
        legend.key.size = unit(0.25, 'cm'))

# Count the number of MAGs of bacteria and fungi, respectively
tax_bin %>% dplyr::select(c("Domain")) %>%
  mutate(number = c(rep(1, nrow(tax_bin)))) %>%
  group_by(Domain) %>%
  summarise(across(everything(), sum)) %>%
  arrange(desc(number)) %>%
  mutate(Domain = gsub("d__", "", Domain)) %>%
  mutate(Domain = factor(Domain, levels = Domain))
# Domain   number
# Bacteria  306
# Archaea   24

# Count the number of MAGs of each Order
mags_num_order_plot <- tax_bin %>% dplyr::select(c("Order")) %>%
  mutate(number = c(rep(1, nrow(tax_bin)))) %>%
  group_by(Order) %>%
  summarise(across(everything(), sum)) %>%
  arrange(desc(number)) %>%
  mutate(Order = gsub("o__", "", Order)) %>%
  mutate(Order = factor(Order, levels = Order)) %>% 
  filter(!Order %in% "") %>% # remove one row with the Order of empty
  slice_head(n = 15) %>% # select the top 15 orders
  ggplot(aes(x = Order, y = number)) +
  geom_bar(stat = "identity", fill = "#ffc15c", width = 0.65) +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  labs(x = "Order", y = "Number of MAGs") +
  theme_classic() +
  theme(axis.title = element_text(color = "black", size = 7),
        axis.text.x = element_text(colour = "black", angle = 90, 
                                   vjust = 0.5, hjust=1, size = 6),
        axis.text.y = element_text(colour = "black", size = 6))
# combine the two plots
mags_num_plot <- cowplot::plot_grid(pie_for_mag_num_phylum, mags_num_order_plot,
                 nrow = 1, labels = c("a", "b"), label_size = 7,
                 rel_widths = c(0.3, 0.7), align = "hv")

# save the plot
# ggsave(file.path(save.dir, "./figs/MAGs/mags_num_plot.pdf"),
#        mags_num_plot, width = 8.9, height = 4.5, units = "in")

# determine the effect size of the permafrost thawing on each MAG
mags_abun <- abundance_tab %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample_name") %>%
  mutate(Gully_id = case_when(grepl("G1", Sample_name) ~ "EB",
                              grepl("G2", Sample_name) ~ "ML",
                              grepl("G3", Sample_name) ~ "RS",
                              grepl("G4", Sample_name) ~ "SLH",
                              grepl("G5", Sample_name) ~ "HSX",
                              grepl("G6", Sample_name) ~ "HH")) %>%
  mutate(Group = case_when(grepl("_C", Sample_name) ~ "Un-collapsed",
                           grepl("_T", Sample_name) ~ "Collapsed")) %>%
  left_join(metadata %>%
              select(Sample_name, Time, Slope, MAP), by = 'Sample_name') %>%
  relocate(c(Group, Time, Slope, MAP, Gully_id), .before = 1) %>%
  select(-Sample_name) %>% data.frame()

mags_scale <- mags_abun %>% 
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size
library(lme4)
library(lmerTest)
#Function for add a significant symbols according to the P values
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
#Function for linear mixed model
summary.model <- function(model){
  presult <- car::Anova(model, type = 2)
  coefs <- coef(summary(model))[, "Estimate"]  ##four coefs
  names(coefs) <- paste0(names(coefs), ".mean")
  
  SEvalues <- coef(summary(model))[, "Std. Error"]  ##standard errors
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  
  tvalues <- coef(summary(model))[, "t value"]  ##t values
  names(tvalues) <- paste0(names(tvalues), ".t")
  
  chisqP <- c(presult[, 1], presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
  
  results <- c(coefs, SEvalues, tvalues, chisqP)
  return(results)
}

lmm.fun <- function(data) {
  df <- NULL
  for(i in 6:ncol(data)) {
    fm <- lmer(data[, i] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = data)
    tmp <- summary.model(fm)
    if (is.null(df)){
      df <- tmp
    } else {
      df <- rbind(df, tmp)
    }
  }
  rownames(df)<-colnames(data)[-c(1:5)]
  df <- df %>% data.frame(check.names = F) %>%
    mutate(sig = p.stars(Group.P), 
           eff.group = case_when(GroupCollapsed.mean > 0 & Group.P < 0.05 ~ "Positive",
                                 GroupCollapsed.mean < 0 & Group.P < 0.05 ~ "Negative",
                                 Group.P > 0.05 ~ "Neutral")) %>%
    mutate(itol_label = case_when(eff.group == "Negative" ~ "TL|0|10|#79ceb8|Negative",
                                  eff.group == "Positive" ~ "TR|0|10|#e95f5c|Positive",
                                  eff.group == "Neutral" ~ ""))
  return(df)
}

mags_result_lmm <- lmm.fun(mags_scale)
# write.table(mags_result_lmm, file = "E:/thermokarst_gully/result/tables/MAGs/diff_in_mags.csv", 
#             sep = ",", col.names = NA, qmethod = "double")


# Exploring the composition of the MAGs within significant positive and negative groups
mags_tax_group_df <- mags_result_lmm %>% rownames_to_column("ID") %>% 
  left_join(tax_bin, by = "ID") %>%
  select(Class, eff.group) %>% 
  filter(eff.group %in% c("Positive", "Negative")) %>%
  mutate(Number = rep(1, nrow(.))) %>%
  summarize(across(everything(), sum), .by = c(Class, eff.group)) %>%
  mutate(Class = sapply(Class, function(x) gsub("c__", '', x, perl = TRUE))) %>%
  mutate(Class = factor(Class),
         Class = fct_reorder(Class, Number, .desc = T))

# loading the appropriate libraries
library(ggplot2)
library(webr)
library(dplyr)
library(RColorBrewer)
PieDonut(mags_tax_group_df, aes(eff.group, Class, count = Number), 
         r0 = 0.45, r1 = 0.9, explode = 1)


# Comparing all microbial traits
library(microtrait)
microtrait_results <- readRDS(file.path(wd_fun, "/MAGs/microtraits/thermokarst_gully.microtraitresults.rds"))

# Normalizing the traits matrices by genome length
microtrait_results_metadata_norm = microtrait_results %>% trait.normalize(normby = "genome_length")

# Inputting other traits including genome size, max growth rate, and so on
mags_result_eff.size <- mags_result_lmm %>% rownames_to_column("Genome_id")
growth_rate_file <- "E:/thermokarst_gully/data/metagenome/MAGs/growth_rate_results.csv"
growth_rate <- read.table(growth_rate_file, header = TRUE, sep = ",", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE)
MAGs_all_aatr <- mags.att %>% left_join(growth_rate, by = "Genome_id") %>% left_join(mags_result_eff.size, by = "Genome_id") %>%
  mutate(Genome_size = round(Genome_size/1000000, 1), `Max growth rate` = 1/d) %>%
  select(Genome_id, Genome_size, GC, `Max growth rate`, GroupCollapsed.mean, eff.group)

# Combining all MAGs attributes
microtrait_results_metadata = microtrait::add.metadata(microtrait_results_metadata_norm, 
                                                       MAGs_all_aatr, genome_metadata_idcol = "Genome_id")
# saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))


strategy_vars <- c("Genome_size", "Max growth rate", "GC",  
                   "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization",
                   "Resource Acquisition:Substrate degradation:simple compound degradation",
                   "Resource Acquisition:Substrate uptake:carbohydrate transport",
                   "Resource Acquisition:Substrate uptake:N compound transport",
                   "Resource Acquisition:Substrate assimilation:C1 compounds",
                   "Resource Acquisition:Substrate assimilation:N compounds",
                   "Resource Use:Chemotrophy:chemolithoautotrophy",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration",
                   "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration",
                   "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation",
                   "Stress Tolerance:General",
                   "Stress Tolerance:Specific:pH stress",
                   "Stress Tolerance:Specific:temperature:low",
                   "Stress Tolerance:Specific:desiccation/osmotic/salt stress",
                   "Stress Tolerance:Specific:oxidative stress")

# writing a function to change the facet labels
strategy_names <- list(
  "Genome_size" = "Genome size (Mb)",
  "Max growth rate" = expression(paste("Max growth rate", " (h"^-1, ")")),
  "GC" = "GC content (%)",
  "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization" = "Complex carbohydrate depolymerization",
  "Resource Acquisition:Substrate degradation:simple compound degradation" = "Simple compound degradation",
  "Resource Acquisition:Substrate uptake:carbohydrate transport" = "Carbohydrate transport",
  "Resource Acquisition:Substrate uptake:N compound transport" = "N compound transport",
  "Resource Acquisition:Substrate assimilation:C1 compounds" = "C1 compounds assimilation",
  "Resource Acquisition:Substrate assimilation:N compounds" = "N compounds assimilation",
  "Resource Use:Chemotrophy:chemolithoautotrophy" = "Chemolithoautotrophy",
  "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration" = "Aerobic respiration", 
  "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration" = "Anaerobic respiration", 
  "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation" = "Fermentation",
  "Stress Tolerance:General" = "General stress tolerance", 
  "Stress Tolerance:Specific:pH stress" = "pH stress",
  "Stress Tolerance:Specific:temperature:low" = "Low temperature stress",
  "Stress Tolerance:Specific:desiccation/osmotic/salt stress" = "Desiccation/osmotic/salt stress", 
  "Stress Tolerance:Specific:oxidative stress" = "Oxidative stress")

facet_labeller <- function(variable,value){
  return(strategy_names[value])
}

# Plot the distribution of the microbial traits
library(gghalves)
library(ggpubr)
# seting the main theme for ploting
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


my_comparisons <- list(c('Negative', 'Positive'))
microtraits_comparison <- microtrait_results_metadata$trait_matrixatgranularity1 %>% 
  select(c("eff.group", all_of(strategy_vars))) %>%
  filter(eff.group %in% c('Negative', 'Positive')) %>%
  gather(mag.attrs, value, -c("eff.group")) %>%
  mutate(eff.group = factor(eff.group, levels = c('Negative', 'Positive'))) %>%
  mutate(mag.attrs = factor(mag.attrs, levels = strategy_vars)) %>%
  ggplot(aes(eff.group, value, fill = eff.group)) +
  geom_half_violin(position = position_nudge(x = 0.25), 
                   side = "r", width = 0.6, color = NA, alpha = 0.75) +
  # geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = eff.group), size = 1.5, 
              width = 0.15, stroke = 0, pch = 21, alpha = 0.75) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.25) +
  stat_compare_means(comparisons = my_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", 
                     bracket.size = 0.5, bracket.width = 0.1,  size = 3.5,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~mag.attrs, scales = "free_y", ncol = 6, labeller = facet_labeller) + #
  main_theme +
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
# ggsave(file.path(save.dir, "./figs/MAGs/microtraits_comparison.pdf"),
#        microtraits_comparison, width = 8.9, height = 4.5, units = "in")

microtraits_comparison

####################################
#######   Max Growth Rate   ########
####################################
library(gRodon)
library(Biostrings)

# Specify the directory containing the .ffn files
input_directory <- "E:/thermokarst_gully/data/metagenome/MAGs/prokka_ffn/"
output_file <- "E:/thermokarst_gully/data/metagenome/MAGs/growth_rate_results.csv"

# Get a list of all .ffn files in the directory
fnn_files <- list.files(input_directory, pattern = "\\.ffn$", full.names = TRUE)

# Initialize an empty data frame to store the results
results <- data.frame(
  File = character(),
  CUBHE = numeric(),
  GC = numeric(),
  GCdiv = numeric(),
  ConsistencyHE = numeric(),
  CUB = numeric(),
  CPB = numeric(),
  FilteredSequences = numeric(),
  nHE = numeric(),
  dCUB = numeric(),
  d = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each .ffn file and calculate the metrics
for (file in ffn_files) {
  cat("Processing:", file, "\n")
  tryCatch({
    # Read the DNA sequences
    genes <- readDNAStringSet(file)
    
    # Identify highly expressed genes (ribosomal proteins)
    highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)
    
    # Predict the growth metrics
    growth_rate <- predictGrowth(genes, highly_expressed, mode = "partial")
    
    # Append the metrics to the results data frame
    results <- rbind(results, data.frame(
      File = basename(file),
      CUBHE = growth_rate$CUBHE,
      GC = growth_rate$GC,
      GCdiv = growth_rate$GCdiv,
      ConsistencyHE = growth_rate$ConsistencyHE,
      CUB = growth_rate$CUB,
      CPB = growth_rate$CPB,
      FilteredSequences = growth_rate$FilteredSequences,
      nHE = growth_rate$nHE,
      dCUB = growth_rate$dCUB,
      d = growth_rate$d,
      LowerCI = growth_rate$LowerCI,
      UpperCI = growth_rate$UpperCI
    ))
  }, error = function(e) {
    cat("Error processing", file, ":", e$message, "\n")
  })
}

# Write the results to a CSV file
write.csv(results, file = output_file, row.names = FALSE)


##########################################
########      Microbial Traits      ######
##########################################
# Set the genomes_dir directory
genomes_dir <- "/data01/kangluyao/thermokarst_gully/binning/prokka_fna"
genomes_files = list.files(genomes_dir, full.names = T, recursive = T, pattern = ".fna$")
out_dir = "/data01/kangluyao/thermokarst_gully/binning/microtraits"

#Determine the number of "cores"
message("Number of cores:", parallel::detectCores(), "\n")


#Here, we use 30 cores to process 330 genomes:
library("tictoc")
tictoc::tic.clearlog()

microtrait_results = extract.traits.parallel(genomes_files, out_dir = dirname(genomes_files), ncores = 30)

#Pull the paths to the corresponding "rds" files as follows:
rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file", mc.cores = 30))

#Building trait matrices from microTrait outputs
# A file named thermokarst_gully.microtraitresults.rds will be generated in your working directory
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = 30)


