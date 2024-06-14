# set work directory
wd_18s_protist <- "E:/thermokarst_gully_18s/result/protist"
wd_18s_animal <- "E:/thermokarst_gully_18s/result/animal"

# set save directory
save.dir_protist <- "E:/thermokarst_gully/data/18S/protist"
if (!dir.exists(save.dir_protist)) {
  dir.create(save.dir_protist, recursive = TRUE)
}
save.dir_animal <- "E:/thermokarst_gully/data/18S/animal"
if (!dir.exists(save.dir_animal)) {
  dir.create(save.dir_animal)
}

# loading packages
library(phyloseq)
library(ape)
library(Biostrings)

# read protist data
## metadata
metadata <- read.table(file.path(wd_18s_protist, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_name)
## otu table
otu <- read.table(file.path(wd_18s_protist, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu <- otu[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu)]]
## tax
tax <- read.table(file.path(wd_18s_protist, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax <- as.matrix(tax)
## tree
tree <- read_tree(file.path(wd_18s_protist, "./otus.nwk"))

## Represent DNA sequence
ref.seqs <- readDNAStringSet(file.path(wd_18s_protist, "./otus.fa"),
                                  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# phyloseq object
otu <- otu_table(otu, taxa_are_rows = TRUE)
tax <- tax_table(tax)
meta_dat <- sample_data(metadata)
phylo <- phyloseq(otu, tax, meta_dat, tree, ref.seqs)
sample_names(phylo) <- metadata$Sample_name
otu <- otu_table(phylo)

phylo_even = rarefy_even_depth(phylo, sample.size = min(colSums(otu)), rngseed = 1234, replace = TRUE)
otu_rare <-otu_table(phylo_even)

# write files
write.table(otu, file.path(save.dir_protist, 'otutab.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.table(otu_rare, file.path(save.dir_protist, 'otutab_rare.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.table(tax, file.path(save.dir_protist, 'taxonomy.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.tree(tree, file.path(save.dir_protist, 'otus.nwk'))
writeXStringSet(ref.seqs, file.path(save.dir_protist, 'otus.fa'), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")



# read animal data
## metadata
metadata <- read.table(file.path(wd_18s_animal, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$Sample_name)
## otu table
otu <- read.table(file.path(wd_18s_animal, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu <- otu[, metadata$Sample_name[metadata$Sample_name %in% colnames(otu)]]
## tax
tax <- read.table(file.path(wd_18s_animal, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax <- as.matrix(tax)
## tree
tree <- read_tree(file.path(wd_18s_animal, "./otus.nwk"))

## Represent DNA sequence
ref.seqs <- readDNAStringSet(file.path(wd_18s_animal, "./otus.fa"),
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
# phyloseq object
otu <- otu_table(otu, taxa_are_rows = TRUE)
tax <- tax_table(tax)
meta_dat <- sample_data(metadata)
phylo <- phyloseq(otu, tax, meta_dat, tree, ref.seqs)
sample_names(phylo) <- metadata$Sample_name
otu <- otu_table(phylo)

phylo_even = rarefy_even_depth(phylo, sample.size = min(colSums(otu)), rngseed = 1234, replace = TRUE)
otu_rare <-otu_table(phylo_even)

# write files
write.table(otu, file.path(save.dir_animal, 'otutab.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.table(otu_rare, file.path(save.dir_animal, 'otutab_rare.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.table(tax, file.path(save.dir_animal, 'taxonomy.txt'),
            sep='\t',  col.names = NA, row.names = T , quote=FALSE)
write.tree(tree, file.path(save.dir_animal, 'otus.nwk'))
writeXStringSet(ref.seqs, file.path(save.dir_animal, 'otus.fa'), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
