# set work directory
setwd('e:/thermokarst_gully')
# set sub directory
wd_16s <- "E:/thermokarst_gully/data/16S/rdp"
wd_its <- "E:/thermokarst_gully/data/ITS"
wd_18s_protist <- "E:/thermokarst_gully/data/18S/protist"
wd_18s_animal <- "E:/thermokarst_gully/data/18S/animal"

# loading packages
library(phyloseq)
library(ape)
library(Biostrings)

# read 16S data
## metadata
metadata <- read.delim(file.path(wd_16s, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_name)
## otu table
otu_16s <- read.table(file.path(wd_16s, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu_16s <- otu_16s[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_16s)]]
## rarefy out table
otu_16s_rare <- read.delim(file.path(wd_16s, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_16s_rare <- otu_16s_rare[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_16s_rare)]]
## tax
tax_16s <- read.table(file.path(wd_16s, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax_16s <- as.matrix(tax_16s)
## tree
tree_16s <- read_tree(file.path(wd_16s, "./otus.nwk"))
## Represent DNA sequence
ref.seqs_16s <- readDNAStringSet(file.path(wd_16s, "./otus.fa"),
                                 format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# 16s phyloseq object
otu_16s <- otu_table(otu_16s, taxa_are_rows = TRUE)
otu_16s_rare <- otu_table(otu_16s_rare, taxa_are_rows = TRUE)
tax_16s <- tax_table(tax_16s)
meta_dat <- sample_data(metadata)
phylo_16s <- phyloseq(otu_16s, tax_16s, meta_dat, tree_16s, ref.seqs_16s)
phylo_16s_rare <- phyloseq(otu_16s_rare, tax_16s, meta_dat, tree_16s, ref.seqs_16s)
sample_names(phylo_16s) <- metadata$Sample_name
sample_names(phylo_16s_rare) <- metadata$Sample_name

# read ITS data
## metadata
metadata <- read.delim(file.path(wd_16s, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_name)
## otu table
otu_its <- read.table(file.path(wd_its, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu_its <- otu_its[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_its)]]
## rarefy out table
otu_its_rare <- read.delim(file.path(wd_its, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_its_rare <- otu_its_rare[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_its_rare)]]
## tax
tax_its <- read.table(file.path(wd_its, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax_its <- as.matrix(tax_its)
## tree
tree_its <- read_tree(file.path(wd_its, "./otus.nwk"))

## Represent DNA sequence
ref.seqs_its <- readDNAStringSet(file.path(wd_its, "./otus.fa"),
                                 format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# its phyloseq object
otu_its <- otu_table(otu_its, taxa_are_rows = TRUE)
otu_its_rare <- otu_table(otu_its_rare, taxa_are_rows = TRUE)
tax_its <- tax_table(tax_its)
# meta_dat <- sample_data(metadata)
phylo_its <- phyloseq(otu_its, tax_its, meta_dat, tree_its, ref.seqs_its)
phylo_its_rare <- phyloseq(otu_its_rare, tax_its, meta_dat, tree_its, ref.seqs_its)
sample_names(phylo_its) <- metadata$Sample_name
sample_names(phylo_its_rare) <- metadata$Sample_name

# read protist data
## metadata
metadata <- read.delim(file.path(wd_16s, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_name)
## otu table
otu_pro <- read.table(file.path(wd_18s_protist, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu_pro <- otu_pro[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_pro)]]
## rarefy out table
otu_pro_rare <- read.delim(file.path(wd_18s_protist, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_pro_rare <- otu_pro_rare[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_pro_rare)]]
## tax
tax_pro <- read.table(file.path(wd_18s_protist, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax_pro <- as.matrix(tax_pro)
## tree
tree_pro <- read_tree(file.path(wd_18s_protist, "./otus.nwk"))

## Represent DNA sequence
ref.seqs_pro <- readDNAStringSet(file.path(wd_18s_protist, "./otus.fa"),
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# Protist phyloseq object
otu_pro <- otu_table(otu_pro, taxa_are_rows = TRUE)
otu_pro_rare <- otu_table(otu_pro_rare, taxa_are_rows = TRUE)
tax_pro <- tax_table(tax_pro)
meta_dat <- sample_data(metadata)
phylo_protist <- phyloseq(otu_pro, tax_pro, meta_dat, tree_pro, ref.seqs_pro)
phylo_protist_rare <- phyloseq(otu_pro_rare, tax_pro, meta_dat, tree_pro, ref.seqs_pro)
sample_names(phylo_protist) <- metadata$Sample_name
sample_names(phylo_protist_rare) <- metadata$Sample_name


# read animal data
## otu table
otu_anim <- read.table(file.path(wd_18s_animal, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu_anim <- otu_anim[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_anim)]]
## rarefy out table
otu_anim_rare <- read.delim(file.path(wd_18s_animal, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_anim_rare <- otu_anim_rare[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu_anim_rare)]]
## tax
tax_anim <- read.table(file.path(wd_18s_animal, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax_anim <- as.matrix(tax_anim)
## tree
tree_anim <- read_tree(file.path(wd_18s_animal, "./otus.nwk"))

## Represent DNA sequence
ref.seqs_anim <- readDNAStringSet(file.path(wd_18s_animal, "./otus.fa"),
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# Protist phyloseq object
otu_anim <- otu_table(otu_anim, taxa_are_rows = TRUE)
otu_anim_rare <- otu_table(otu_anim_rare, taxa_are_rows = TRUE)
meta_dat <- sample_data(metadata)
tax_anim <- tax_table(tax_anim)
phylo_animal <- phyloseq(otu_anim, tax_anim, meta_dat, tree_anim, ref.seqs_anim)
phylo_animal_rare <- phyloseq(otu_anim_rare, tax_anim, meta_dat, tree_anim, ref.seqs_anim)
sample_names(phylo_animal) <- metadata$Sample_name
sample_names(phylo_animal_rare) <- metadata$Sample_name

# read functional table
wd_fun <- file.path(getwd(),"data/metagenome")
# if (!dir.exists(wd_fun)) {
#   dir.create(wd_fun)
# }
ko_tpm_table <- read.delim(file.path(wd_fun, "./eggnog.KEGG_ko.raw.txt"), 
                           header = T, row.names = 1, sep = "\t")
ko_tpm_table <- ko_tpm_table[, metadata$Sample_name[metadata$Sample_name %in% colnames(ko_tpm_table)]]
