# set work directory
setwd('e:/thermokarst_gully')
# set sub directory
wd_16s <- "E:/thermokarst_gully/data/16S/rdp"

# loading packages
library(tidyverse)
library(phyloseq)
library(ape)
library(Biostrings)

# read 16S data
## metadata
metadata <- read.delim(file.path(wd_16s, "./metadata.txt"), header = T, sep = "\t")

rownames(metadata) <- (metadata$Sample_name)
## otu table
otu_16s <- read.delim(file.path(wd_16s, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu_16s <- otu_16s[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_16s)]]
## rarefy out table
otu_16s_rare <- read.delim(file.path(wd_16s, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_16s_rare <- otu_16s_rare[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_16s_rare)]]
## tax
tax_16s <- read.delim(file.path(wd_16s, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax_16s <- as.matrix(tax_16s)
# 16s phyloseq object
otu_16s <- otu_table(otu_16s, taxa_are_rows = TRUE)
otu_16s_rare <- otu_table(otu_16s_rare, taxa_are_rows = TRUE)
tax_16s <- tax_table(tax_16s)
meta_dat <- sample_data(metadata)
phylo_16s <- phyloseq(otu_16s, tax_16s, meta_dat)
phylo_16s_rare <- phyloseq(otu_16s_rare, tax_16s, meta_dat)
sample_names(phylo_16s) <- metadata$Sample_name
sample_names(phylo_16s_rare) <- metadata$Sample_name

# Reading functional gene table
wd_fun <- file.path(getwd(),"data/metagenome")
# if (!dir.exists(wd_fun)) {
#   dir.create(wd_fun)
# }
ko_tpm_table <- read.delim(file.path(wd_fun, "./eggnog.KEGG_ko.raw.txt"), 
                           header = T, row.names = 1, sep = "\t")
ko_tpm_table <- ko_tpm_table[, metadata$Sample_name[metadata$Sample_name %in% colnames(ko_tpm_table)]]

# Reading genome data
# abundance_tab.file <- file.path(wd_fun, "MAGs/bin_abundance_coverm.txt")
abundance_tab.file <- "E:/thermokarst_gully/result/MAGs/coverM_75_abundance.tsv"
tax.file <- file.path(wd_fun, "MAGs/tax.txt")
mags.att.file <- file.path(wd_fun, "MAGs/MAGs_attributes.txt")
mags_abun_tab <- read.delim(abundance_tab.file, header = TRUE, sep = "\t", row.names = 1,
                            as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                            check.names = FALSE)
mags_abun_tab <- mags_abun_tab[, metadata$Sample_name]
mags_tax <- read.delim(tax.file, header = TRUE, sep = "\t", as.is = TRUE, 
                      stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)

mags.att <- read.delim(mags.att.file, header = TRUE, sep = "\t", row.names = 1,
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE)

