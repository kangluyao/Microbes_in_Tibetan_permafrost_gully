# set work directory
setwd('e:/thermokarst_gully')
wd_its <- file.path(getwd(),"data/ITS")
if (!dir.exists(wd_its)) {
  dir.create(wd_its)
}
save.dir <- file.path(getwd(),"result")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}
# loading packages
library(phyloseq)
library(ape)
library(Biostrings)
# read data

## metadata
metadata <- read.delim(file.path(wd_its, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_id)

## tree
# tree <- read_tree(file.path(wd_its, "./otus.nwk"))
## otu table
otu <- read.delim(file.path(wd_its, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu <- otu[, metadata$Sample_id[metadata$Sample_id %in% colnames(otu)]]
## rarefy out table
otu_rare <- read.delim(file.path(wd_its, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_rare <- otu_rare[, metadata$Sample_id[metadata$Sample_id %in% colnames(otu_rare)]]
## tax
tax <- read.delim(file.path(wd_its, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax <- as.matrix(tax)
# phyloseq object
otu <- otu_table(otu, taxa_are_rows = TRUE)
otu_rare <- otu_table(otu_rare, taxa_are_rows = TRUE)
tax <- tax_table(tax)
meta_dat <- sample_data(metadata)
phylo <- phyloseq(otu, tax, meta_dat)
phylo_rare <- phyloseq(otu_rare, tax, meta_dat)
sample_names(phylo) <- metadata$Sample_name
sample_names(phylo_rare) <- metadata$Sample_name
otu <- otu_table(phylo)