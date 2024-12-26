ko_tpm_table

ko_9_levels.file <- "E:/thermokarst_gully/data/metagenome/ko_9_levels.txt"
ko_9_levels <-  read.table(ko_9_levels.file, header = TRUE, sep = "\t", 
                           as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                           check.names = FALSE)

ko_9_levels_merged <- ko_tpm_table %>% rownames_to_column("KO") %>%
  inner_join(ko_9_levels, by = "KO") %>% 
  write.csv(file = "E:/thermokarst_gully/data/multifun/ko_9_levels_merged_ko_number.csv")

L2_cat <- c("Carbohydrate metabolism", "Lipid metabolism", "Amino acid metabolism", 
        "Metabolism of cofactors and vitamins", "Replication and repair",
        "Membrane transport", "Signal transduction", "Cell growth and death", "Cell motility", "Infectious disease: bacterial",
        "Infectious disease: parasitic", "Drug resistance: antimicrobial", "Drug resistance: antineoplastic")

L3_cat <- c("Methane metabolism", "Nitrogen metabolism", "Sulfur metabolism")

L2_df <- ko_9_levels_merged %>% 
  select(-c("KO", "Description", "L1_ID", "L1", "L2_ID", "L3_ID", "L3", "PathwayID")) %>% 
  distinct() %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% group_by(L2) %>%
  summarize(across(everything(), sum)) %>% filter(L2 != "") %>% column_to_rownames("L2") %>% t() %>%
  as.data.frame() %>% rownames_to_column("Sample_id") %>%
  mutate(Group = c(case_when(grepl("_C", Sample_id) ~ "Uncollapsed",
                             grepl("_T", Sample_id) ~ "Collapsed")),
         Gully_id = c(case_when(grepl("G1", Sample_id) ~ "G1",
                                grepl("G2", Sample_id) ~ "G2",
                                grepl("G3", Sample_id) ~ "G3",
                                grepl("G4", Sample_id) ~ "G4",
                                grepl("G5", Sample_id) ~ "G5",
                                grepl("G6", Sample_id) ~ "G6")))
L2_df 









# determine the effect size of the permafrost thawing for the L2 pathways
pathways_sel <- colnames(L2_df)[-c(1, 39, 40)]

L2_path_scale <- L2_df  %>% 
  select(c("Group", "Gully_id", all_of(pathways_sel))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
L2_path_S1 <- sapply(3:ncol(L2_path_scale), function(j) {
  message("Now j=", j, " in ", ncol(L2_path_scale), ". ", date())
  if (length(unique(L2_path_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(L2_path_scale[, j] ~ Group + (1 | Gully_id), data = L2_path_scale)
    
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
colnames(L2_path_S1)<-colnames(L2_path_scale)[-c(1:2)]
data.frame(L2_path_S1)





p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_L2_path_comparison <- L2_path_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(pathways_sel))) %>%
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
  # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#ffdb00") +
  # annotate("rect", xmin = 4.5, xmax = 9.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#f19837") +
  # annotate("rect", xmin = 9.5, xmax = 19.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#e56eee") +
  # annotate("rect", xmin = 19.5, xmax = 22.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#5fb236") +
  # annotate("rect", xmin = 22.5, xmax = 31.5, ymin = -2, ymax = 2, alpha = 0.1, fill = "#1ca9c9") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
ggsave(file.path(save.dir, "./figs/diversity/single_L2_path_comparison.pdf"),
       single_L2_path_comparison, width = 3.5, height = 5.25, units = "in")
single_L2_path_comparison







library(tidyverse)
wd_fun <- file.path(getwd(),"data/metagenome")
sel_ko <- read_csv(file.path(wd_fun, 'CN_ko_input.csv'), col_names = T)
ko_tax_TPM_table <- read_delim(file.path(wd_fun, 'fun/parse_dat.txt'), col_names = T)

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



