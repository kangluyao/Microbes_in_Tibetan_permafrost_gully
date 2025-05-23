## Diversity analysis

#Determine the alpha diversity including **Observed**, **Chao1**, **Shannon**, **Simpson**, **Faith index**, and beta diversity.

# estimate the alpha diversity using rarefied otu table
alpha_div_16s <- estimate_richness(phylo_16s_rare, measures = c("Observed", "Chao1", 'Shannon', 'Simpson'))
alpha_div_its <- estimate_richness(phylo_its_rare, measures = c("Observed", "Chao1", 'Shannon', 'Simpson'))

col_div_names <- c("Sample_name", 'Gully_id', 'Group', "Bacterial Richness", "Fungal Richness")
div_table <- cbind(metadata[, c("Sample_name", 'Gully_id', 'Group')], alpha_div_16s$Observed, alpha_div_its$Observed) %>%
  rename_with(~ col_div_names) %>%
  select(c("Sample_name", 'Gully_id', 'Group'), col_div_names) %>%
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')))


# Calculate alpha diversity measures for a specific taxon at a specified rank.
require("phyloseq")
require("dplyr")
# You can pass any parameters that you normally pass to `estimate_richness`
estimate_diversity_for_phyla <- function(ps, taxon_name, tax_rank = "Phylum", ...){
  # Subset to taxon of interest
  tax_tbl <- as.data.frame(tax_table(ps))
  keep <- tax_tbl[,tax_rank] == taxon_name
  keep[is.na(keep)] <- FALSE
  ps_phylum <- prune_taxa(keep, ps)
  
  # Calculate alpha diversity and generate a table
  alpha_diversity <- estimate_richness(ps_phylum, ...)
  alpha_diversity$taxon <- taxon_name
  alpha_diversity$sample_id <- row.names(alpha_diversity)
  return(alpha_diversity)
}

# Estimate alpha diversity for each phylum
estimate_alpha_fun <- function(ps) {
  alpha <- data.frame()
  phyla <- get_taxa_unique(ps, taxonomic.rank = 'Phylum')
  phyla <- phyla[!is.na(phyla)]
  for (phylum in phyla){
    a <- estimate_diversity_for_phyla(ps = ps, 
                                      taxon_name = phylum,
                                      measure = c("Shannon", "Observed"))
    alpha <- rbind(alpha, a)
  }
  alpha <- alpha %>% select(Observed, taxon, sample_id) %>%
    pivot_wider(names_from = taxon, values_from = Observed)
  return(alpha)
}

richness_16s_phylum <- estimate_alpha_fun(phylo_16s)
richness_its_phylum <- estimate_alpha_fun(phylo_its)

# estimate the richness at the class level for bacterial proteobacteria and soil animals
estimate_diversity_for_class <- function(ps, taxon_name, tax_rank = "Class", ...){
  # Subset to taxon of interest
  tax_tbl <- as.data.frame(tax_table(ps))
  keep <- tax_tbl[,tax_rank] == taxon_name
  keep[is.na(keep)] <- FALSE
  ps_phylum <- prune_taxa(keep, ps)
  
  # Calculate alpha diversity and generate a table
  alpha_diversity <- estimate_richness(ps_phylum, ...)
  alpha_diversity$taxon <- taxon_name
  alpha_diversity$sample_id <- row.names(alpha_diversity)
  return(alpha_diversity)
}

# Estimate alpha diversity for each phylum
estimate_fun_class <- function(ps) {
  alpha <- data.frame()
  phyla <- get_taxa_unique(ps, taxonomic.rank = 'Class')
  phyla <- phyla[!is.na(phyla)]
  for (class in phyla){
    a <- estimate_diversity_for_class(ps = ps, 
                                      taxon_name = class,
                                      measure = c("Shannon", "Observed"))
    alpha <- rbind(alpha, a)
  }
  alpha <- alpha %>% select(Observed, taxon, sample_id) %>%
    pivot_wider(names_from = taxon, values_from = Observed)
  return(alpha)
}


richness_16s_class <- estimate_fun_class(phylo_16s)

# Combine all diversity table
div_tab_all <- cbind(div_table, richness_16s_phylum[, colnames(richness_16s_phylum) %in% names_16s],
                     richness_16s_class[, colnames(richness_16s_class) %in% names_16s],
                     richness_its_phylum[, colnames(richness_its_phylum) %in% names_its])

div_tab_all <- div_tab_all[, c("Gully_id", "Group", "Bacterial Richness", names_16s, "Fungal Richness", names_its)]

# write.table(div_tab_all, file.path(save.dir, 'tables/div_table_all.csv'),
#             sep=',',  col.names = T, row.names = F, quote = FALSE)

# Linear mixed models test the effect of collapsed and gully_id
tax_index <- c("Bacterial Richness", names_16s, "Fungal Richness", names_its)
div_scale <- div_tab_all %>% 
  cbind(metadata[, c("Time", "Slope", "MAP")]) %>%
  select(all_of(c("Group", "Gully_id",  "Time", "Slope", "MAP", tax_index))) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Un-collapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
div_S1 <- sapply(6:ncol(div_scale), function(j) {
  if (length(unique(div_scale[, j])) < 3) {
    result <- rep(NA, 23)
  } else {
    fm1 <- lmer(div_scale[, j] ~ Group + Time + Slope + MAP + (1 | Gully_id), data = div_scale)
    
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
colnames(div_S1)<-colnames(div_scale)[-c(1:5)]
data.frame(div_S1)

p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
# Create a plot
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 7),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 7),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.title = element_text(colour = 'black', size = 7),
        legend.text = element_text(colour = 'black', size = 6),
        legend.key.size = unit(0.5, 'cm'))
all_tax_div_comparison <- div_S1 %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>% filter(variables %in% tax_index) %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(tax_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.375,  colour = "gray2")+
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.1),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("#79ceb8", "grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2.2, 2.2)) +
  coord_flip() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = -2.2, ymax = 2.2, alpha = 0.1, fill = "#e56eee") +
  annotate("rect", xmin = 9.5, xmax = 24.5, ymin = -2.2, ymax = 2.2, alpha = 0.1, fill = "#5fb236") +
  main_theme +
  theme(legend.position = "none",
        strip.background = element_rect(fill = c("#FFF6E1")),
        # axis.text.y = element_blank()
  )

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path("E:/thermokarst_gully/result2/all_tax_div_comparison.pdf"),
#        all_tax_div_comparison, width = 2.5, height = 5, units = "in")
all_tax_div_comparison
