library(tidyverse)
wd_fun <- file.path(getwd(),"data/metagenome")
sel_ko <- read_csv(file.path(wd_fun, 'CN_ko_input.csv'), col_names = T)
ko_tax_TPM_table <- read_delim(file.path(wd_fun, './fun/parse_dat.txt'), col_names = T)

C_N_ko_count_tab <- ko_count_table[rownames(ko_count_table) %in% sel_ko$KO, ]
C_N_ko_count_tab <- cbind(KO = rownames(C_N_ko_count_tab), C_N_ko_count_tab)

C_N_ko_tpm_tab <- ko_tpm_table[rownames(ko_tpm_table) %in% sel_ko$KO, ]
C_N_ko_tpm_tab <- cbind(KO = rownames(C_N_ko_tpm_tab), C_N_ko_tpm_tab)


C_N_ko_tax_tab <- ko_tax_TPM_table[ko_tax_TPM_table$KO %in% sel_ko$KO, ]
C_N_ko_tax_path_tab <- merge(C_N_ko_tax_tab, sel_ko, by="KO", all = T)

nrow(C_N_ko_tax_tab)
nrow(C_N_ko_tax_path_tab)
colnames(C_N_ko_tax_path_tab)


# heatmap
C_N_path1 <- C_N_ko_count_tab %>% 
  inner_join(sel_ko, c("KO" = "KO")) %>% 
  select(-c(1)) %>%
  group_by(Category, Subcategory, Enzyme_protein_encoded) %>%
  summarise(across(everything(), sum)) %>% data.frame() %>%
  mutate(Subcategory = factor(Subcategory, levels = unique(sel_ko$Subcategory), ordered = T))
C_N_path1 %>% filter(Category == "Carbon") %>% select(c(2, 3))
C_N_path2 <- C_N_path1 %>% select(c(3:63)) %>% 
  tibble::column_to_rownames('Enzyme_protein_encoded')
  
# C_N_path2 <- C_N_ko_count_tab %>% 
#   inner_join(sel_ko, c("KO" = "KO")) %>% 
#   select(-c(Subcategory, KO)) %>%
#   group_by(Enzyme_protein_encoded) %>%
#   summarise(across(everything(), sum)) %>% data.frame() %>%
#   pivot_longer(cols = -Enzyme_protein_encoded, names_to = 'sample_id', values_to = 'counts') %>%
#   mutate(layer = gsub("_.+$", "", sample_id))

# path_result <- C_N_path2 %>% select(-sample_id) %>% 
#   group_by(Enzyme_protein_encoded, layer) %>%
#   summarise_all(list(mean = mean, sd = sd, se = ~sd(./sqrt(.))))
library(reshape2)
plot_dat <- apply(C_N_path2, MARGIN = 2, FUN = scale)
rownames(plot_dat) <- rownames(C_N_path2)
plot_dat <- t(plot_dat)
plot_dat <-  setNames(melt(plot_dat), c('samples', 'Enzyme_protein_encoded', 'values'))
plot_dat$Enzyme_protein_encoded <- factor(plot_dat$Enzyme_protein_encoded, ordered = T,
                                          levels = unique(sel_ko$Enzyme_protein_encoded))
ggplot(plot_dat) +
  geom_tile(aes(x = as.factor(samples), y = as.factor(Enzyme_protein_encoded), fill = values)) +
  xlab('') + ylab('') + theme_linedraw() + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10)) +
  scale_fill_gradientn(colours = c('white','#00A087FF'), values = c(0,1)) + 
  #facet_grid(~Ecosystem, scales = 'free_x', space = 'free_x') +
  theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
        panel.spacing = unit(0, "lines"), strip.background = element_blank())

##############################################################################################################
# edgeR analysis
save.dir <- file.path(getwd(),"result")
# 1. alternately edgeR analysi.
library("edgeR")

#The field in the class definition file that defines the classes of the data.
data_classes <- "Layer"

# Load Expression Data
RNASeq <- C_N_path2

# Load subtype information
RNASeq[1:5, 1:5]
classDefinitions_RNASeq <- meta_dat
classDefinitions_RNASeq[1:5, 1:3]

# Filter Data (RNA-seq read counts are converted to CPM values 
# and genes with CPM > 1 in at least 50 of the samples are 
# retained for further study, a gene mush have at least 50 measurements 
# with more than 1 CPM in one of the classes to be included in the analysis)
cpms <- cpm(RNASeq)
keep <- rowSums(cpms > 1) >= 1
counts <- RNASeq[keep,]

# Normalization and Dispersion
# create data structure to hold counts and subtype information for each sample.
d <- DGEList(counts = counts, group =  factor(classDefinitions_RNASeq$Group))
d <- calcNormFactors(d)
#calculate dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
# Testing for differential gene expression
et <- exactTest(object = d)
# topTags() function is useful to extract the table with adjusted p values (FDR).
top_degs = topTags(object = et, n = "Inf")
top_degs
#get the indices of scored dataset that have FDR < 0.05
select_genes = top_degs$table %>%  filter(FDR < 0.05)
#output how many genes there are in the set that have FDR < 0.05
nrow(select_genes)



#compare PL layer to the remaining two layers.
contrast_pl <- makeContrasts(
  plvsrest = "classesPL-(classesSUR + classesSUB)/2", levels = modelDesign)
fit_glm <- glmFit(d, modelDesign)
plvsrest <- glmLRT(fit_glm , contrast = contrast_pl)
tt_plvsrest <- topTags(plvsrest, n = nrow(d))
#get the indices of scored dataset that have FDR < 0.05
select_pl_genes = tt_plvsrest$table %>% filter(FDR < 0.05)
#output how many genes there are in the set that have FDR < 0.05
nrow(select_pl_genes)

# arrange the log fold change table
C_names <- c("Alpha-amylase", "Glucoamylase", "Pullulanase", "Isopullulanase",
             "Arabinofuranosidase", "Beta_mannanase", "Xylanase", "Xylose isomerase",
             "Beta-glucosidase",  "Endoglucanase", "Exoglucanase",
             "Acetylglucosaminidase", "Endochitinase", "Exochitinase",
             "Pectinase", 
             "Aryl-aldehyde oxidase", "Isocitrate lyase", "Limonene-1, 2-epoxide hydrolase", "Malate synthase", "Vanillate demethylase", 
             "Glyoxal oxidase", "Phenol oxidase (tyrosinase)",
             "Methyl coenzyme M reductase", "Particulate methane monooxygenase", "Soluable methane monooxygenase",
             "Carbon monoxide dehydrogenase", "Tetrahydrofolate formylase", "ATP citrate lyase", 
             "Propionyl-CoA carboxylase", "RuBisCo")
N_names <- c("amoA",	"amoB",	"amoC",	"hao",	"nxrA",	"nxrB", "narG",	"narH",	"narI", "napA",	
             "napB", "nirK", "nirS",	"norB",	"norC",	"nosZ", "nrfA",	"nrfH", "nirB",	"nirD",	
             "nasA",	"nasB", "narB",	"NR", "NIT-6", "nirA",	"nifD", "nifK", "nifH", "nrtA",	"nrtB", "nrtC",	
             "nrtD", "nmo", "gdh_K00261", "gdh_K00262", "gdh_K15371", "glsA", "ureA", "ureC", "glnA")
logFC_table <- data.frame(pathway = rownames(top_degs$table), logFC = top_degs$table$logFC)

#heatmap
# set the plot theme
main_theme = theme_linedraw() + 
  theme(panel.grid=element_blank(), 
        strip.text = element_text(colour = 'black', size = 12),
        strip.background = element_rect(colour = 'grey', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12, angle = 90, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))

p_C_enrich <- logFC_table %>% filter(pathway %in% C_names) %>%
  mutate(pathway = factor(pathway, levels = C_names, ordered = T)) %>%
  mutate(process = rep('process', nrow(.))) %>%
  ggplot(aes(x = pathway, y = process, fill = logFC)) +
  geom_tile() + 
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +  # low="#2C7BB6", mid="white", high="#D7191C" or low = "#009E73", mid = "white", high = "#E69F00"
  geom_text(aes(label = round(logFC, 2)), 
            color = "black", size = 2.5) +
  labs(x = 'Pathway', y = NULL, fill = "logFC")+
  main_theme

plot.name <- paste0(save.dir, "/figs/metagenome/edgeR_Carbon_heatmap.pdf")
print(plot.name)
cairo_pdf(filename = plot.name, width = 6.9, height = 2.2, onefile = TRUE)
p_C_enrich
dev.off()

p_N_enrich <- logFC_table %>% filter(pathway %in% N_names) %>%
  mutate(pathway = factor(pathway, levels = rev(N_names), ordered = T)) %>%
  mutate(process = rep('process', nrow(.))) %>%
  ggplot(aes(x = pathway, y = process, fill = logFC))+
  geom_tile() + 
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +  # low="#2C7BB6", mid="white", high="#D7191C" or low = "#009E73", mid = "white", high = "#E69F00"
  geom_text(aes(label = round(logFC, 2)), 
            color = "black", size = 2)+
  labs(x = 'Pathway', y = NULL, fill = "logFC")+
  main_theme





















# Draw barplot with grouping & stacking
library(ggplot2)
library(dplyr)
d1 <- C_N_ko_tax_path_tab %>% select(c('Subcategory', 'Genu', 
                                       grep('SUR', colnames(C_N_ko_tax_path_tab)),
                                       grep('SUB', colnames(C_N_ko_tax_path_tab)),
                                       grep('PL', colnames(C_N_ko_tax_path_tab)))) %>%
  pivot_longer(cols = -c(Subcategory, Genu),
               names_to = "sample", values_to = 'abundance') %>%
  mutate(Layers = gsub("_.+$", "", sample)) %>%
  filter(Subcategory == 'ANRA') %>% select('sample', 'Layers', 'Genu', 'abundance') %>%
  group_by(Genu, sample, Layers) %>%
  summarise(across(everything(), sum)) %>%
  group_by(Layers) %>%
  arrange(sample, desc(abundance)) %>%
  ggplot(.,aes(x = sample, y = abundance, fill = Genu)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  facet_grid(~ Layers, scales="free_x") +
  theme(legend.position = 'none')


# 2. DESeq2 analysis
C_N_ko_tab <- ko_count_table[rownames(ko_count_table) %in% sel_ko$KO, ]
C_N_ko_tab <- cbind(KO = rownames(C_N_ko_tab), C_N_ko_tab)
merged_dat <- dplyr::inner_join(C_N_ko_tab, sel_ko, 
                                c("KO" = "KO"))
C_N_path2 <- merged_dat %>% dplyr::select(-c(1)) %>% 
  group_by(Subcategory, Enzyme_protein_encoded) %>%
  dplyr::summarise(across(everything(), sum)) %>% data.frame() %>%
  dplyr::select(c(2:68)) %>%
  tibble::column_to_rownames('Enzyme_protein_encoded')

meta_dat$layer = as.factor(meta_dat$layer)

KEGG_dds <- DESeqDataSetFromMatrix(countData=round(C_N_path2), colData=meta_dat, design=~layer)
KEGG_deseq <- DESeq(KEGG_dds)
KEGG_res <- results(KEGG_deseq, contrast=c("layer", "SUR", "PL"))
KEGG_res$padj[is.na(KEGG_res$padj)] = 1
KEGG_significant = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]

# write.csv(as.data.frame(KEGG_res), file = 'Data/KEGG_deseq_results.csv')

############################################################################################################
# 2. volcano plot
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
