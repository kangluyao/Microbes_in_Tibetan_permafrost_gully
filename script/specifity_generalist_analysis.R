#determine the generalist and specialist species
library(EcolUtils)
spec_gen_16s <- spec.gen(phylo_16s@otu_table %>% t(), niche.width.method = "levins", perm.method = "quasiswap",
                         n = 1000, probs = c(0.05, 0.95))
# write.table(spec_gen_16s, file.path(save.dir, "/tables/spe_gen/spec_gen_16s_95.csv"), sep = ",")
spec_gen_its <- spec.gen(phylo_its@otu_table %>% t(), niche.width.method = "levins", perm.method = "quasiswap",
                         n = 1000, probs = c(0.05, 0.95))
# write.table(spec_gen_its, file.path(save.dir, "/tables/spe_gen/spec_gen_its_95.csv"), sep = ",")
spec_gen_pro <- spec.gen(phylo_protist@otu_table %>% t(), niche.width.method = "levins", perm.method = "quasiswap",
                         n = 1000, probs = c(0.05, 0.95))
# write.table(spec_gen_pro, file.path(save.dir, "/tables/spe_gen/spec_pro_95.csv"), sep = ",")
spec_gen_animal <- spec.gen(phylo_animal@otu_table %>% t(), niche.width.method = "levins", perm.method = "quasiswap",
                            n = 1000, probs = c(0.05, 0.95))
# write.table(spec_gen_animal, file.path(save.dir, "/tables/spe_gen/spec_animal_95.csv"), sep = ",")

#explore the biotic homogenelization and heteronazition
spec_gen_16s <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_16s_95.csv"), sep = ",", header = T)
spec_id_16s <- spec_gen_16s %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_16s <- spec_gen_16s %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_16s <- spec_gen_16s %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% spec_id_16s, ]
gen_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% gen_id_16s, ]
neu_16s_comm <- phylo_16s@otu_table[rownames(phylo_16s@otu_table) %in% neu_id_16s, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_16s_dist <-vegdist(t(spec_16s_comm), "bray")
tax_gen_16s_dist <-vegdist(t(gen_16s_comm), "bray")
tax_neu_16s_dist <-vegdist(t(neu_16s_comm), "bray")

# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_spec_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_gen_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_neu_16s_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_16s_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_16s_data <- data.frame(Group = similar_spec_16s_data$Group,
                           Gully_id = similar_spec_16s_data$Gully_id, 
                           distance_16s_spec = similar_spec_16s_data$distance,
                           distance_16s_gen = similar_gen_16s_data$distance,
                           distance_16s_neu = similar_neu_16s_data$distance)

#Descriptive statistics for all diversity indexes
distance_16s_index <- c("distance_16s_spec", "distance_16s_gen", "distance_16s_neu")
distance_16s_stats <- similar_16s_data %>% group_by(Group) %>%
  get_summary_stats(distance_16s_index, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  arrange(variable, Group)
distance_16s_stats

library(lme4)
library(lmerTest)
# determine the effect size of the permafrost thawing for the diversity indexes
dist_16s_scale <- similar_16s_data %>% 
  select(c("Group", "Gully_id", distance_16s_index)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
dist_16s_S1 <- sapply(3:ncol(dist_16s_scale), function(j) {
  message("Now j=", j, " in ", ncol(dist_16s_scale), ". ", date())
  if (length(unique(dist_16s_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(dist_16s_scale[, j] ~ Group + (1 | Gully_id), data = dist_16s_scale)
    
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
colnames(dist_16s_S1)<-colnames(dist_16s_scale)[-c(1:2)]
data.frame(dist_16s_S1)

#plot
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
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_dist_16s_comparison <- dist_16s_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(distance_16s_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.7,  colour = "gray2")+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom") +
  theme(legend.position = "none", 
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.1), 
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_rect(colour = "white"))

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_div_comparison, width = 2.7, height = 5, units = "in")
single_dist_16s_comparison





################################################################

#explore the biotic homogenelization and heteronazition
spec_gen_its <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_its_95.csv"), sep = ",", header = T)
spec_id_its <- spec_gen_its %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_its <- spec_gen_its %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_its <- spec_gen_its %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% spec_id_its, ]
gen_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% gen_id_its, ]
neu_its_comm <- phylo_its@otu_table[rownames(phylo_its@otu_table) %in% neu_id_its, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_its_dist <-vegdist(t(spec_its_comm), "bray")
tax_gen_its_dist <-vegdist(t(gen_its_comm), "bray")
tax_neu_its_dist <-vegdist(t(neu_its_comm), "bray")

# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_spec_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_gen_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_neu_its_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_its_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_its_data <- data.frame(Group = similar_spec_its_data$Group,
                               Gully_id = similar_spec_its_data$Gully_id, 
                               distance_its_spec = similar_spec_its_data$distance,
                               distance_its_gen = similar_gen_its_data$distance,
                               distance_its_neu = similar_neu_its_data$distance)

#Descriptive statistics for all diversity indexes
distance_its_index <- c("distance_its_spec", "distance_its_gen", "distance_its_neu")
distance_its_stats <- similar_its_data %>% group_by(Group) %>%
  get_summary_stats(distance_its_index, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  arrange(variable, Group)
distance_its_stats

library(lme4)
library(lmerTest)
# determine the effect size of the permafrost thawing for the diversity indexes
dist_its_scale <- similar_its_data %>% 
  select(c("Group", "Gully_id", distance_its_index)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
dist_its_S1 <- sapply(3:ncol(dist_its_scale), function(j) {
  message("Now j=", j, " in ", ncol(dist_its_scale), ". ", date())
  if (length(unique(dist_its_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(dist_its_scale[, j] ~ Group + (1 | Gully_id), data = dist_its_scale)
    
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
colnames(dist_its_S1)<-colnames(dist_its_scale)[-c(1:2)]
data.frame(dist_its_S1)

#plot
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
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_dist_its_comparison <- dist_its_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(distance_its_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.7,  colour = "gray2")+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom") +
  theme(legend.position = "none", 
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.1), 
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_rect(colour = "white"))

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_div_comparison, width = 2.7, height = 5, units = "in")
single_dist_its_comparison


########################################################################
#explore the biotic homogenelization and heteronazition
spec_gen_protist <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_pro_95.csv"), sep = ",", header = T)
spec_id_protist <- spec_gen_protist %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_protist <- spec_gen_protist %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_protist <- spec_gen_protist %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_protist_comm <- phylo_protist@otu_table[rownames(phylo_protist@otu_table) %in% spec_id_protist, ]
gen_protist_comm <- phylo_protist@otu_table[rownames(phylo_protist@otu_table) %in% gen_id_protist, ]
neu_protist_comm <- phylo_protist@otu_table[rownames(phylo_protist@otu_table) %in% neu_id_protist, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_protist_dist <-vegdist(t(spec_protist_comm), "bray")
tax_gen_protist_dist <-vegdist(t(gen_protist_comm), "bray")
tax_neu_protist_dist <-vegdist(t(neu_protist_comm), "bray")

# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_spec_protist_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_protist_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_gen_protist_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_protist_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_neu_protist_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_protist_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_protist_data <- data.frame(Group = similar_spec_protist_data$Group,
                               Gully_id = similar_spec_protist_data$Gully_id, 
                               distance_protist_spec = similar_spec_protist_data$distance,
                               distance_protist_gen = similar_gen_protist_data$distance,
                               distance_protist_neu = similar_neu_protist_data$distance)

#Descriptive statistics for all diversity indexes
distance_protist_index <- c("distance_protist_spec", "distance_protist_gen", "distance_protist_neu")
distance_protist_stats <- similar_protist_data %>% group_by(Group) %>%
  get_summary_stats(distance_protist_index, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  arrange(variable, Group)
distance_protist_stats

library(lme4)
library(lmerTest)
# determine the effect size of the permafrost thawing for the diversity indexes
dist_protist_scale <- similar_protist_data %>% 
  select(c("Group", "Gully_id", distance_protist_index)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
dist_protist_S1 <- sapply(3:ncol(dist_protist_scale), function(j) {
  message("Now j=", j, " in ", ncol(dist_protist_scale), ". ", date())
  if (length(unique(dist_protist_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(dist_protist_scale[, j] ~ Group + (1 | Gully_id), data = dist_protist_scale)
    
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
colnames(dist_protist_S1)<-colnames(dist_protist_scale)[-c(1:2)]
data.frame(dist_protist_S1)

#plot
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
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_dist_protist_comparison <- dist_protist_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(distance_protist_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.7,  colour = "gray2")+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom") +
  theme(legend.position = "none", 
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.1), 
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_rect(colour = "white"))

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_div_comparison, width = 2.7, height = 5, units = "in")
single_dist_protist_comparison


#############################################################################
#explore the biotic homogenelization and heteronazition
spec_gen_anim <- read.table(file.path(save.dir, "/tables/spe_gen/spec_gen_animal_95.csv"), sep = ",", header = T)
spec_id_anim <- spec_gen_anim %>% filter(sign == "SPECIALIST") %>% pull(OTU)
gen_id_anim <- spec_gen_anim %>% filter(sign == "GENERALIST") %>% pull(OTU)
neu_id_anim <- spec_gen_anim %>% filter(sign == "NON SIGNIFICANT") %>% pull(OTU)

spec_anim_comm <- phylo_animal@otu_table[rownames(phylo_animal@otu_table) %in% spec_id_anim, ]
gen_anim_comm <- phylo_animal@otu_table[rownames(phylo_animal@otu_table) %in% gen_id_anim, ]
neu_anim_comm <- phylo_animal@otu_table[rownames(phylo_animal@otu_table) %in% neu_id_anim, ]


library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
tax_spec_anim_dist <-vegdist(t(spec_anim_comm), "bray")
tax_gen_anim_dist <-vegdist(t(gen_anim_comm), "bray")
tax_neu_anim_dist <-vegdist(t(neu_anim_comm), "bray")

# difference in taxonomic variance among Group
vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
# Assuming vars is defined somewhere earlier in your code
similar_spec_anim_data <- lapply(vars, function(x) usedist::dist_subset(tax_spec_anim_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_gen_anim_data <- lapply(vars, function(x) usedist::dist_subset(tax_gen_anim_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_neu_anim_data <- lapply(vars, function(x) usedist::dist_subset(tax_neu_anim_dist, grep(x, metadata$Sample_name, value = TRUE))) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  gather("tem_group", "distance") %>%
  cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
        Group = rep(c('Uncollapsed', 'Collapsed'), each = 10, times = 6)) %>%
  select(-tem_group) %>%
  mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
         Group = factor(Group, levels = c('Uncollapsed', 'Collapsed')))

similar_anim_data <- data.frame(Group = similar_spec_anim_data$Group,
                                   Gully_id = similar_spec_anim_data$Gully_id, 
                                   distance_anim_spec = similar_spec_anim_data$distance,
                                   distance_anim_gen = similar_gen_anim_data$distance,
                                   distance_anim_neu = similar_neu_anim_data$distance)

#Descriptive statistics for all diversity indexes
distance_anim_index <- c("distance_anim_spec", "distance_anim_gen", "distance_anim_neu")
distance_anim_stats <- similar_anim_data %>% group_by(Group) %>%
  get_summary_stats(distance_anim_index, type = "common") %>% #or using type = "mean_sd"
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  arrange(variable, Group)
distance_anim_stats

library(lme4)
library(lmerTest)
# determine the effect size of the permafrost thawing for the diversity indexes
dist_anim_scale <- similar_anim_data %>% 
  select(c("Group", "Gully_id", distance_anim_index)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(Group = factor(Group, levels = c("Uncollapsed", "Collapsed"))) %>%
  select(where(~ !any(is.na(.))))

# codes for calculating the effect size refer to wu et al. 2022:https://github.com/Linwei-Wu/warming_soil_biodiversity.
dist_anim_S1 <- sapply(3:ncol(dist_anim_scale), function(j) {
  message("Now j=", j, " in ", ncol(dist_anim_scale), ". ", date())
  if (length(unique(dist_anim_scale[, j])) < 3) {
    result <- rep(NA, 10)
  } else {
    fm1 <- lmer(dist_anim_scale[, j] ~ Group + (1 | Gully_id), data = dist_anim_scale)
    
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
colnames(dist_anim_S1)<-colnames(dist_anim_scale)[-c(1:2)]
data.frame(dist_anim_S1)

#plot
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
p.stars <- function(p.values) {
  unclass(symnum(p.values, corr = FALSE, 
                 na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))}
single_dist_anim_comparison <- dist_anim_S1 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "variables") %>%
  mutate(sig = as.vector(unlist(lapply(Group.P, p.stars)))) %>%
  mutate(variables = factor(variables, levels = rev(distance_anim_index))) %>%
  mutate(colour = case_when(GroupCollapsed.mean <= 0 & Group.P <= 0.05 ~ "Negative",
                            GroupCollapsed.mean > 0 & Group.P <= 0.05 ~ "Positvie",
                            Group.P > 0.05 ~ "Neutral")) %>%
  ggplot(aes(x = variables, y = GroupCollapsed.mean, colour = colour)) +
  geom_hline(aes(yintercept = 0), size = 0.7,  colour = "gray2")+
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = GroupCollapsed.mean - GroupCollapsed.se, 
                    ymax = GroupCollapsed.mean + GroupCollapsed.se), 
                width = 0, position = position_dodge(width = 0.7), cex = 0.9) +
  geom_text(aes(label = sig, x = variables, y = (GroupCollapsed.mean/abs(GroupCollapsed.mean))*(abs(GroupCollapsed.mean) + GroupCollapsed.se)*1.2),
            position = position_dodge(0.1), vjust = 0.55) +
  labs(x = NULL, y = "Effect size") +
  scale_color_manual(values=c("grey", "#e95f5c")) +
  scale_y_continuous(expand = c(0, 0), limit = c(-2, 2)) +
  theme_bw() + coord_flip() + scale_x_discrete(position = "bottom") +
  theme(legend.position = "none", 
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.1, "lines"), axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.1), 
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_rect(colour = "white"))

# if (!dir.exists(file.path(save.dir, "figs/env/"))) {
#   dir.create(file.path(save.dir, "figs/env/"))
# }
# ggsave(file.path(save.dir.multifunc, "./single_div_comparison.pdf"),
#        single_div_comparison, width = 2.7, height = 5, units = "in")
single_dist_anim_comparison










comm.tab.bin<-ceiling(t(otu_16s)/max(t(otu_16s)))
plot(colMeans(t(otu_16s)),colSums(comm.tab.bin)/dim(comm.tab.bin)[1],
     col=spec_gen$sign,pch=19,log="x",xlab="Abundance",ylab="Occurrence")
legend("bottomright",levels(spec_gen$sign),col=1:3,pch=19,inset=0.01,cex=0.7)

spe_tax <- c("GENERALIST", "NON SIGNIFICANT", "SPECIALIST")
spe_gen_data <- merge(spec_gen_animal, phylo_animal@otu_table, by = 0, all = TRUE) %>%
  select(-c(1:5)) %>% group_by(sign) %>%
  summarise(across(everything(), sum)) %>% 
  column_to_rownames("sign") %>%
  t() %>% as.data.frame() %>%
  mutate(Gully_id = metadata$Gully_id) %>%
  mutate(Group = factor(metadata$Group, ordered = T, levels = c("Uncollapsed", "Collapsed")))

spe_gen_summary <- spe_gen_data %>%
  group_by(Group) %>%
  get_summary_stats(spe_tax, type = "common") %>% #or using type = "mean_sd"
  arrange(variable, Group)
spe_gen_summary
