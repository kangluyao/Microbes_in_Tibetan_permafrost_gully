#Set working and saving directory
setwd('e:/thermokarst_gully/')
save.dir <- file.path(getwd(),"result")

#Loading packages and reading data
pacman::p_load(phyloseq, ape, vegan, Biostrings, microbiome, tidyverse, rstatix)
#reading data
source("script/read_data_all.R")

# extract the community and tree files from the phylo projects
comm_16s <- t(otu_table(phylo_16s))
tree_16s <- phy_tree(phylo_16s)

comm_its <- t(otu_table(phylo_its))
tree_its <- phy_tree(phylo_its)

comm_pro <- t(otu_table(phylo_protist))
tree_pro <- phy_tree(phylo_protist)

comm_anim <- t(otu_table(phylo_animal))
tree_anim <- phy_tree(phylo_animal)

#iCAMP
setwd(save.dir)
library(iCAMP)
set.seed(123)
icamp_16s <- icamp.big(comm = comm_16s, tree = tree_16s, pd.wd = paste0(save.dir,"/tables/null_model/16s_all_1"),
                      ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24,
                      rand = 1000, nworker = 8)
head(icamp_16s$CbMPDiCBraya)
# write.csv(icamp_16s$CbMPDiCBraya,
#           file.path(save.dir, './tables/null_model/16s_all/iCAMP.process.CbMPDiCBraya.csv'))


treatment <- data.frame(metadata$Group)
rownames(treatment) <- rownames(metadata)
colnames(treatment) <- "Group"
treatment$Group <- factor(treatment$Group, levels = c("Collapsed", "Uncollapsed"))
icampbt_16s <- icamp.boot(icamp.result = icamp_16s$CbMPDiCBraya, treat = treatment, 
                   rand.time = 50, compare = TRUE, between.group = F, ST.estimation = TRUE)


icampbt_16s$compare


data("icamp.out")
data("example.data")
treatment=example.data$treat
rand.time=20 # usually use 1000 for real data.
icampbt=icamp.boot(icamp.result = icamp.out$bNRIiRCa, treat = treatment,
                   rand.time = rand.time, compare = TRUE,
                   between.group = TRUE, ST.estimation = TRUE)




set.seed(123)
icamp_its <- icamp.big(comm = comm_its, tree = tree_its, pd.wd = paste0(save.dir,"/tables/null_model/its_all"),
                        ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24,
                        rand = 1000, nworker = 8)
head(icamp_its$CbMPDiCBraya)
# write.csv(icamp_its$CbMPDiCBraya,
#           file.path(save.dir, './tables/null_model/its_all/iCAMP.process.CbMPDiCBraya.csv'))

set.seed(123)
icamp_pro <- icamp.big(comm = comm_pro, tree = tree_pro, pd.wd = paste0(save.dir,"/tables/null_model/pro_all"),
                        ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24,
                       rand = 1000, nworker = 8)
head(icamp_pro$CbMPDiCBraya)
# write.csv(icamp_pro$CbMPDiCBraya,
#           file.path(save.dir, './tables/null_model/pro_all/iCAMP.process.CbMPDiCBraya.csv'))

set.seed(123)
icamp_anim <- icamp.big(comm = comm_anim, tree = tree_anim, pd.wd = paste0(save.dir,"/tables/null_model/anim_all"),
                       ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24,
                       rand = 1000, nworker = 8)
head(icamp_anim$CbMPDiCBraya)
# write.csv(icamp_anim$CbMPDiCBraya,
#           file.path(save.dir, './tables/null_model/anim_all/iCAMP.process.CbMPDiCBraya.csv'))