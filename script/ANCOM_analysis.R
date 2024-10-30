
## ANCOM test
library(microbiomeMarker)
ancom_fun <- function(ps) {
  ps_genus <- tax_glom(ps, "Genus", NArm = TRUE)
  ps_genus <- subset_taxa(ps_genus, Genus != "Unassigned")
  output <- run_ancom(
    ps_genus,
    group = "Group",
    confounders = character(0),
    taxa_rank = "Genus",
    transform = "log10",
    norm = "TSS",
    norm_para = list(),
    p_adjust = "BH",
    pvalue_cutoff = 0.01,
    W_cutoff = 0.75
  )
  tax_tab <- data.frame(tax_table(ps_genus)[, c(1, 2, 6)])
  colnames(tax_tab) <- c("Kingdom", "Phylum", "feature")
  marker_tab <- data.frame(marker_table(output))
  marker_tab <- inner_join(marker_tab, tax_tab, by = "feature")
  return(marker_tab)
}
marker_16s <- ancom_fun(phylo_16s)
marker_its <- ancom_fun(phylo_its)
marker_pro <- ancom_fun(phylo_protist)
marker_anim <- ancom_fun(phylo_animal)

marker_tab <- rbind(marker_16s, marker_its, marker_pro, marker_anim)

marker_tab <- marker_tab %>%
  mutate(dif = case_when(
    enrich_group == "Uncollapsed" ~ "(-)", 
    enrich_group == "Collapsed" ~ "(+)"))
marker_tab



## bar plot
# phylum: c('Acidobacteriota','Actinobacteriota','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Others,'Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota')
# colors: c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')
ANCOM_barplot <- ggplot(marker_tab, aes(x = ef_CLR_diff_mean, 
                                        y = reorder(as.factor(feature), ef_CLR_diff_mean), 
                                        fill = Phylum)) + 
  geom_bar(stat='identity', position = 'dodge', width = 0.7) + 
  geom_text(aes(label = dif, x = ef_CLR_diff_mean + 0.15), 
            position = position_dodge(0.7), vjust = 0,size = 1.5) +
  scale_fill_manual(values = c('#F39B7FFF','#DC0000FF','#00A087FF','grey','dimgrey','#3C5488FF','#8491B4FF', 
                               "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", 
                               "#FF9896", "#9467BD", "#C5B0D5")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2.4)) +
  ylab('') + xlab('CLR mean difference') +
  theme_bw()+
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.25), 
        axis.text.y = element_text(colour = 'black', size = 6, face = "italic"),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.position = "inside", 
        legend.position.inside =  c(0.7, 0.25),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_rect(colour = "white"))
ANCOM_barplot  