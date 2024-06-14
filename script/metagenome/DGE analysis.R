save.dir <- file.path(getwd(),"result")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}
#loading packages
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(edgeR)
library(DESeq2)

# 1. edgeR analysis. reference:https://www.reneshbedre.com/blog/edger-tutorial.html
#The field in the class definition file that defines the classes of the data.
data_classes <- "Group"

# Load Expression Data
RNASeq <- ko_count_table

# Load group information
RNASeq[1:5, 1:5]
classDefinitions_RNASeq <- metadata
classDefinitions_RNASeq[1:5, 1:4]

# Filter Data (RNA-seq read counts are converted to CPM values 
# and genes with CPM > 1 in at least 50 of the samples are 
# retained for further study, a gene mush have at least 50 measurements 
# with more than 1 CPM in one of the classes to be included in the analysis)
cpms <- cpm(RNASeq)
keep <- rowSums(cpms > 1) >= 1
counts <- RNASeq[keep,]

# Normalization and Dispersion
# create data structure to hold counts and group information for each sample.
d <- DGEList(counts = counts, group =  factor(classDefinitions_RNASeq$Group))
# Filter out the genes with low counts
# #filterByExpr function to remove the low count genes. 
# #This function keeps the genes with a worthwhile counts (at least 10 read counts) in a minimal number of samples.
# keep <- filterByExpr(y = d)
# d <- d[keep, , keep.lib.sizes=FALSE]
# Normalization and effective library size
d <- calcNormFactors(d)
# Model fitting and estimating dispersions
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

# Testing for differential gene expression
et <- exactTest(object = d)
# topTags() function is useful to extract the table with adjusted p values (FDR).
top_degs <- topTags(object = et, n = "Inf")
top_degs
#Get a summary DGE table (returns significant genes with absolute log fold change at least 1 and adjusted p value < 0.05)
summary(decideTests(object = et, lfc = 1))
# write.csv(as.data.frame(top_degs), file = file.path(save.dir, "tables/kegg_Collapsed_vs_control_edge.csv"))
# DGE Visualization
# Create a volcano plot
# set the plot theme
main_theme = theme_linedraw() + 
  theme(panel.grid=element_blank(), 
        strip.text = element_text(colour = 'black', size = 12),
        strip.background = element_rect(colour = 'grey', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        legend.position = c(0.15, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))
data <- top_degs %>% as.data.frame() %>% 
  mutate(Expression = case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))
ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  main_theme +
  theme(legend.background = element_blank(),
        legend.key = element_blank())


data <- top_degs %>% as.data.frame() %>%
  mutate(
    Significance = case_when(
      abs(logFC) >= log(2) & FDR <= 0.05 & FDR > 0.01 ~ "FDR 0.05", 
      abs(logFC) >= log(2) & FDR <= 0.01 & FDR > 0.001 ~ "FDR 0.01",
      abs(logFC) >= log(2) & FDR <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )

ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  main_theme +
  theme(legend.background = element_blank(),
        legend.key = element_blank())



# 2. DESeq2 analysis
group_df <- metadata
rownames(group_df) <- metadata$Sample_name # make sure the colnames of count table match the rownames of group_df.
group_df$Group <- factor(group_df$Group, levels = c("Control", "Collapsed"))
KEGG_dds <- DESeqDataSetFromMatrix(countData = round(ko_count_table + 1), 
                                   colData = group_df, design = ~ Group)
KEGG_deseq <- DESeq(KEGG_dds)
KEGG_res <- results(KEGG_deseq, contrast = c("Group", 'Collapsed', 'Control'))
KEGG_res$padj[is.na(KEGG_res$padj)] = 1
KEGG_significant = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]
length(KEGG_significant)
write.csv(as.data.frame(KEGG_res), file = file.path(save.dir, "tables/kegg_Collapsed_vs_control_deseq.csv"))
# Volcano plot
data <- KEGG_res %>% as.data.frame() %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
ggplot(data, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"padj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 


EnhancedVolcano(KEGG_res,
                lab = rownames(KEGG_res),
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey30", "grey30", "#23A671"),
                x = 'log2FoldChange',
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                y = 'padj')
ggsave('3_Functional_analysis/3_1_Enrichment/KEGG_enriched_layer.pdf', width = 7, height = 7)

#t treatment, c control
kegg_sig_fun <- function(t, c) {
  KEGG_res <- results(KEGG_deseq, contrast = c("Group", t, c))
  KEGG_res$padj[is.na(KEGG_res$padj)] = 1
  KEGG_sig_in_t = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]
  KEGG_sig_in_c = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange < -1)]
  KEGG_sig <- list(KEGG_sig_in_t, KEGG_sig_in_c)
  names(KEGG_sig) <- c(t, c)
  return(KEGG_sig)
}
kegg_sig <- kegg_sig_fun('Collapsed', 'Control')

# 3. KEGG decoder
kegg_decoder_fun <- function(layer, kegg_sig) {
  samples <- grep(layer, meta_dat$sample_id, value = T)
  kegg_decoder_df = data.frame(Sample = c(), KEGG = c())
  for (kegg_id in kegg_sig) {
    for (sample in samples) {
      if (ko_count_table[rownames(ko_count_table) == kegg_id, sample] > 0) {
        kegg_decoder_df = rbind(kegg_decoder_df, data.frame(Sample = sample, KEGG = kegg_id))
      }
    }
  }
  return(kegg_decoder_df)
}
kegg_decoder_sur_df <- kegg_decoder_fun('SUR', select_sur_genes)
kegg_decoder_sub_df <- kegg_decoder_fun('SUB', select_sub_genes)
kegg_decoder_pl_df <- kegg_decoder_fun('PL', select_pl_genes)

kegg_decoder_sur_df$Sample <- sapply(kegg_decoder_sur_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_sub_df$Sample <- sapply(kegg_decoder_sub_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_pl_df$Sample <- sapply(kegg_decoder_pl_df$Sample, FUN = function(x){sub('_', '', x)})

# write.table(kegg_decoder_sur_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_sur_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_sub_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_sub_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_pl_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_pl_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)

# in the Data folder using the CryoBiome-kegg-decoder env: 
# python KEGG-decoder.py -i KEGG_decoder_sur_df_in.tsv -o KEGG_decoder_sur_df_output.tsv
# python KEGG-decoder.py -i KEGG_decoder_sub_df_in.tsv -o KEGG_decoder_sub_df_output.tsv
# python KEGG-decoder.py -i KEGG_decoder_pl_df_in.tsv -o KEGG_decoder_pl_df_output.tsv
# python -m pip install xxxx
KEGG_decoder_sur <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_sur_df_output.tsv"), 
                               sep = "\t", header = TRUE)
KEGG_decoder_sub <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_sub_df_output.tsv"),
                               sep = "\t", header = TRUE)
KEGG_decoder_pl <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_pl_df_output.tsv"),
                              sep = "\t", header = TRUE)
KEGG_decoder <- rbind(KEGG_decoder_sur, KEGG_decoder_sub, KEGG_decoder_pl)

colsums <- colSums(KEGG_decoder[, !colnames(KEGG_decoder) %in% c("Function")])
KEGG_decoder <- KEGG_decoder[, c("Function", names(colsums[colsums > 0]))]
names(KEGG_decoder)[names(KEGG_decoder) == "Function"] <- "Sample"
KEGG_decoder <- reshape2::melt(KEGG_decoder, variable.name = "Function", value.name = "Completion")

# Renaming pathways for plotting
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "TCA.Cycle"] <- "TCA Cycle"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "NAD.P.H.quinone.oxidoreductase"] <- "NAD(P)H-quinone oxidoreductase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "X3.Hydroxypropionate.Bicycle"] <- "3-Hydroxypropionate Bicycle"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "basic.endochitinase.B"] <- "Basic endochitinase B"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "DMS.dehydrogenase"] <- "DMS dehydrogenase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "thiamin.biosynthesis"] <- "Thiamin biosynthesis"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "transporter..vitamin.B12"] <- "Transporter: vitamin B12"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Flagellum"] <- "Flagellum"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Photosystem.II"] <- "Photosystem II"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Photosystem.I"] <- "Photosystem I"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Cytochrome.b6.f.complex"] <- "Cytochrome b6/f complex"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Entner.Doudoroff.Pathway"] <- "Entner-Doudoroff Pathway"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Mixed.acid..Ethanol..Acetate.to.Acetylaldehyde"] <- "Mixed acid: Ethanol, Acetate to Acetylaldehyde"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Mixed.acid..Ethanol..Acetyl.CoA.to.Acetylaldehyde..reversible."] <- "Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Naphthalene.degradation.to.salicylate"] <- "Naphthalene degradation to salicylate"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Curli.fimbriae.biosynthesis"] <- "Curli fimbriae biosynthesis"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Competence.related.core.components"] <- "Competence-related core components"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "CP.lyase.operon"] <- "CP-lyase operon"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.III.Secretion"] <- "Type III Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.IV.Secretion"] <- "Type IV Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.Vabc.Secretion"] <- "Type Vabc Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "serine"] <- "Serine"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "tyrosine"] <- "Tyrosine"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "starch.glycogen.degradation"] <- "Starch/glycogen degradation"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "end.product.astaxanthin"] <- "End-product astaxanthin"

KEGG_decoder$mean_comp = vapply(1:nrow(KEGG_decoder), function(x) mean(KEGG_decoder$Completion[KEGG_decoder$Function == KEGG_decoder$Function[x]]), FUN.VALUE = numeric(1))
KEGG_decoder$Function <- factor(KEGG_decoder$Function, levels=unique((KEGG_decoder$Function)[order(KEGG_decoder$mean_comp)]))
meta_dat$Sample <- sapply(meta_dat$sample_id, FUN = function(x){sub('_', '', x)})
KEGG_decoder$layer = sapply(1:nrow(KEGG_decoder), function(x) meta_dat$layer[meta_dat$Sample == KEGG_decoder$Sample[x]])
KEGG_decoder$Sample <- factor(KEGG_decoder$Sample, levels = meta_dat$Sample, ordered = T)
KEGG_decoder$layer <- factor(KEGG_decoder$layer, levels = c('SUR', 'SUB', 'PL'), ordered = T)

# plot
library(ggplot2)
p <- ggplot(KEGG_decoder) + 
  geom_tile(aes(x = as.factor(Sample), y = as.factor(Function), fill = Completion)) +
  xlab("") + ylab("") + 
  theme_linedraw() + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10)) +
  scale_fill_gradientn(colours = c("white", "#00A087FF"), values = c(0, 1)) + 
  facet_grid(~layer, scales = "free_x", space = "free_x") + 
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), 
        axis.text.y = element_text(size = 6), 
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave(file.path(save.dir, '/figs/metagenome/kegg_enrichment_edgeR.pdf'), 
       p, width = 10, height = 5)

############################################################################################################
# 4. enrichment analysis with clusterProfiler package
library(clusterProfiler)
KO_pathway <- read.csv("e:/permafrost/data/metagenome/fun/keggKO_ko_metabolism.csv", header = T)
termgene <- KO_pathway[, c(2, 1)]
termname <- KO_pathway[, c(2, 3)]
enrichpathway_sur <- enricher(select_sur_genes, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_sub <- enricher(select_sub_genes, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_pl <- enricher(select_pl_genes, TERM2GENE = termgene, TERM2NAME = termname)

write.csv(enrichpathway_sur, file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_sur_edgeR.csv"))
write.csv(enrichpathway_sub,  file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_sub_edgeR.csv"))
write.csv(enrichpathway_pl,  file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_pl_edgeR.csv"))
