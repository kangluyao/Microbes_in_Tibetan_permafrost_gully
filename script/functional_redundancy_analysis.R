##' Functional redundancy via Contribution Evenness
##'
##' @export
##' @description Calculates functional redundancy via the Royalty method, as described in INSERT REFERENCE HERE.
##' @param abundance.matrix A traditional taxa abundance matrix. Rows correspond to different samples and columns correspond to different taxa.
##' @param trait.matrix #A matrix corresponding to traits within taxa. Columns correspond to traits while rows correspond to taxa 
##' @param q.range The diversity order to evaluate. The default setting evaluates diversity orders 0 through 2.
##' @details ANYTHING THAT NEEDS TO GO INTO THE DETAILS SECTION
##' @return Returns a list, where individual elements correspond to samples in the abundance matrix. The list includes an estimate of functional redundancy and functional diversity for each diversity order, q, as well as the species richness.
##' @examples # 
##' abundance.matrix <- read.csv('data/MAG_abundance_table.csv', row.names = 1 ) 
##' abundance.matrix <- round(abundance.matrix/min(abundance.matrix[abundance.matrix>0])) 
##' sample.effort <- min(rowSums(abundance.matrix))
##' abundance.matrix <- vegan::rrarefy(abundance.matrix,sample.effort)
##' abundance.matrix <- sweep(abundance.matrix,1,rowSums(abundance.matrix),'/')
##' trait.matrix <- read.csv('data/MAG_enzyme_gene_copies.csv', row.names = 1 )
##' fr <- royalty_fr(abundance.matrix, trait.matrix, q = 0.5)
##' fr <- tidyr::separate(fr,sample,into = c("site","size_fraction","depth"), sep = '_')
##' ggplot2::ggplot(fr,aes(x=trait,y=fr,color=depth))+geom_boxplot()+ylim(0,0.7)  
##' @references ADD MANUSCRIPT REFERENCE HERE

royalty_fr <- function(abundance.matrix, trait.matrix, q = 0.5) {
  
  # ADD SOME CODE TO MAKE SURE THAT PARAMS ARE VALID
  if(!all(apply(trait.matrix,1,FUN=is.numeric) == TRUE)) {
    stop("The trait.matrix argument must be a numeric matrix")
  }
  if(!all(apply(abundance.matrix,1,FUN=is.numeric) == TRUE)) {
    stop("The abundance.matrix argument must be a numeric matrix")
  }
  
  
  fr_single_sample <- function(abundance.vector,trait.level, q) {
    
    indx <- abundance.vector > 0
    abundance.vector <- abundance.vector[indx]
    trait.level <- trait.level[indx]
    
    trait.contribution <- abundance.vector * trait.level
    trait.contribution <- trait.contribution / sum(trait.contribution)
    
    S <- length(abundance.vector)
    
    fr_sample <- data.frame(fr = NaN,
                            D = NaN,
                            q = NaN,
                            S = NaN)
    
    if (!is.na(sum(trait.contribution))) {
      if (q == 1) {
        D.q <- exp(-1 * sum(trait.contribution[trait.contribution > 0] * log(trait.contribution[trait.contribution >0])))
      } else {
        D.q <- sum(trait.contribution[trait.contribution > 0] ^ q) ^ (1 / (1 - q))
      }
      
      fr <- (D.q - 1) / (S - 1)
      fr_sample$D <- D.q
      fr_sample$fr <- fr
    } else {
      fr_sample$D <- 0
    }
    fr_sample$q <- q
    fr_sample$S <- S
    
    return(fr_sample)
  }
  
  fr_multiple_sample <- function(abundance.matrix, trait.level, q) {
    
    fr <- apply(X = abundance.matrix,
                MARGIN = 1,
                FUN = fr_single_sample,
                trait.level = trait.level,
                q = q)
    
    fr <- cbind(dplyr::bind_rows(fr),
                sample=names(fr))
    return(fr)
  }
  
  trait.taxa.names<-rownames(trait.matrix)
  trait.names<-colnames(trait.matrix)
  abundance.names <- colnames(abundance.matrix)
  sample.names <- rownames(abundance.matrix)
  overlap.names<-intersect(trait.taxa.names,abundance.names)
  
  n.taxa <- length(overlap.names)
  
  if(n.taxa != length(trait.taxa.names) | n.taxa != length(abundance.names)) {
    warning('The taxa in the trait matrix and abundance matrix do not match.\nUsing only taxa which appear in both.')
  }
  
  trait.matrix <- trait.matrix[overlap.names,] 
  trait.matrix <- t(trait.matrix) 
  colnames(trait.matrix) <- overlap.names
  rownames(trait.matrix) <- trait.names
  abundance.matrix <- abundance.matrix[,overlap.names]
  
  fr_out <- apply(X = trait.matrix, 
                  MARGIN = 1,
                  FUN = fr_multiple_sample,
                  abundance.matrix = abundance.matrix,
                  q = q)
  
  fr_out <- lapply(1:length(fr_out),function(x,fr_out,trait.list) {cbind(fr_out[[x]],trait=trait.list[x])},fr_out,names(fr_out))
  
  fr_out <- dplyr::bind_rows(fr_out)
  
  
  if(sum(is.na(fr_out$fr))>0) {
    warning('Some communities do not have any members contributing to a community-aggregated parameter.\nFunctional redundancy is not defined.')
  }
  
  
  return(fr_out)
}

# read tada
fr_dat <- read.table(file = "E:/thermokarst_gully/data/metagenome/kraken2/parse_dat_count.txt", 
                     header = T, row.names = NULL, sep = "\t")
nrow(fr_dat)
fr_dat[1:5, 1:9]

# function mapping
KO_9leves <- read.table(file = "E:/thermokarst_gully/data/metagenome/ko_9_levels.txt", 
                        header = T, row.names = NULL, sep = "\t")
nrow(KO_9leves)
KO_9leves[1:5, 1:9]

# mapping the KOs into L3 level
fr_dat <- fr_dat %>% left_join(KO_9leves, by = "KO")
nrow(fr_dat)
fr_dat[1:25, c(1:9, 69:75)]

# prepare the abundance table
species.abund <- fr_dat %>% select(c(8:68)) %>%
  group_by(Species) %>%
  summarise(across(everything(), sum)) %>%
  filter(!Species %in% c("", "Unassigned")) %>%
  column_to_rownames(var = "Species") %>% 
  apply(., 2, function(x) x/sum(x)) %>%
  t()

species.abund[1:5, 1:5]
nrow(species.abund)
ncol(species.abund)


# prepare the traits table
trait.levels <- fr_dat %>% select(c(8:68, 75)) %>%
  pivot_longer(cols = -c("L3", "Species"), names_to = "Sample", values_to = "Value") %>%
  select(-c("Sample")) %>%
  group_by(Species, L3) %>%
  summarise(across(everything(), sum)) %>%
  filter(!Species %in% c("", "Unassigned")) %>%
  pivot_wider(names_from = L3, values_from = Value) %>%
  column_to_rownames(var = "Species") %>%
  mutate_all(~replace(., is.na(.), 0))

trait.levels[1:5, 1:5]
nrow(trait.levels)
ncol(trait.levels)

# Caculate the functional reduntancy
fr_table <- royalty_fr(species.abund, trait.levels)
fr_table[1:100, 1:6]

fr <- fr_table %>%
  mutate(Gully_id = sapply(strsplit(sample, "_"), '[', 1)) %>%
  mutate(Group = case_when(
    grepl("_C", sample) ~ "Control",
    grepl("_T", sample) ~ "Collapsed"))

library(ggpubr)
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 6),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.position = "none")

# box plot
my_comparisons <- list(c('Control', 'Collapsed'))
ggplot2::ggplot(fr, aes(x = Group, y = fr)) + 
  geom_boxplot(width = 0.5, aes(fill = Group), outlier.shape = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  stat_compare_means(comparisons = my_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", bracket.size = 0.5,
                     size = 3.5, tip.length = 0.00, method = "wilcox.test") +
  labs(x = 'Group', y = "Functional redundancy", fill= 'Group') +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1))


# test the difference in alpha diversity
library(vegan)
richness_div <- vegan::estimateR(fr_dat %>% select(c(8:68)) %>%
                                   group_by(Species) %>%
                                   summarise(across(everything(), sum)) %>%
                                   filter(!Species %in% c("", "Unassigned")) %>%
                                   column_to_rownames(var = "Species") %>% 
                                   round(., digits = 0) %>%
                                   t()) %>% t()


richness_table <- richness_div %>% as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  mutate(Gully_id = sapply(strsplit(Sample, "_"), '[', 1)) %>%
  mutate(Group = case_when(
    grepl("_C", Sample) ~ "Control",
    grepl("_T", Sample) ~ "Collapsed"))

ggplot2::ggplot(richness_table, aes(x = Group, y = S.obs)) + 
  geom_boxplot(width = 0.5, aes(fill = Group), outlier.shape = NA) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  stat_compare_means(comparisons = my_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", bracket.size = 0.5,
                     size = 3.5, tip.length = 0.00, method = "wilcox.test") +
  labs(x = 'Group', y = "Richness", fill= 'Group') +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  main_theme +
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1))

library(rstatix)
data.frame(Richness = richness_table$S.obs, meanFunction = multi_ave_func_df$meanFunction) %>%
  ggplot(aes(x = Richness, y = meanFunction)) +
  geom_point(size = 2, alpha = 0.8, aes(colour = "#5cc3e8")) +
  geom_smooth(method = "lm", size = 1, se = T, colour = 'black') +
  # scale_colour_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  # scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  stat_poly_line(colour = 'black') +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), 
                                 after_stat(p.value.label), sep = "*\", \"*")),
               size = 2.4) +
  xlab("\nDiversity") + 
  ylab("Average Value of Standardized Functions\n") +
  main_theme + theme(legend.position = "none")


### PERMANOVA test
library(usedist)
item_groups <- rep(rep(c("Control", "Collapsed"), each = 5), 6)
tax_16s_dist <- as.matrix(vegdist(t(otu_16s), "bray")) %>%
  usedist::dist_groups(., item_groups)
phy_dist <- as.matrix(UniFrac(phylo_16s, weighted = TRUE, normalized = TRUE, parallel = T, fast = TRUE)) %>%
  usedist::dist_groups(., item_groups)
fun_dist <- as.matrix(vegdist(t(ko_tpm_table), "bray")) %>%
  usedist::dist_groups(., item_groups)


dist_dat <- data.frame(tax_16s_dist, phy_dist[ , 6], fun_dist[ , 6])
colnames(dist_dat) <- c("Item1", "Item2", "Group1", "Group2", "Label", "Taxnomic_distance", "Phylogenetic_distance", "Functional_distance")


p_linear_tax_fun <- dist_dat %>%
  filter(Label %in% c("Within Control", "Within Collapsed")) %>%
  mutate(Group1 = factor(Group1, levels = c('Control', 'Collapsed'))) %>%
  ggplot(aes(x = Taxnomic_distance, y = Functional_distance, fill = Group1)) +
  geom_point(size = 3.5, alpha = 0.8, aes(colour = Group1)) +
  geom_smooth(method = "lm", size=1, se = T, colour = 'black') +
  scale_colour_manual(values = c("#f8766d", "#a3a500", "#00b0f6")) +
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00b0f6")) +
  stat_poly_line(colour = 'black') +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"), colour = Group1)) +
  ylab("Functional distance")+xlab("Taxnomic distance") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(color = 'black', size = 14),
        axis.ticks.length = unit(0.2, "lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        strip.text = element_text(size = 14),
        legend.position = c(0.15, 0.65))
p_linear_tax_fun




p_linear_phy_fun <- dist_dat %>%
  filter(Label %in% c("Within Control", "Within Collapsed")) %>%
  mutate(Group1 = factor(Group1, levels = c('Control', 'Collapsed'))) %>%
  ggplot(aes(x = Phylogenetic_distance, y = Functional_distance, fill = Group1)) +
  geom_point(size = 3.5, alpha = 0.8, aes(colour = Group1)) +
  geom_smooth(method = "lm", size=1, se = T, colour = 'black') +
  scale_colour_manual(values = c("#f8766d", "#a3a500", "#00b0f6")) +
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00b0f6")) +
  stat_poly_line(colour = 'black') +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"), colour = Group1)) +
  ylab("Functional distance")+xlab("Phylogenetic distance") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(color = 'black', size = 14),
        axis.ticks.length = unit(0.2, "lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        strip.text = element_text(size = 14),
        legend.position = c(0.15, 0.65))
p_linear_phy_fun







uniqueness <- function(comm, dis, tol = 1e-8, abundance = TRUE){
  
  if(!is.null(colnames(comm)) & !is.null(attributes(dis)$Labels)){
    if(any(!colnames(comm)%in%attributes(dis)$Labels)) stop("One or several species in comm are not in dis; check species names in comm and in dis")
    else dis <- as.dist(as.matrix(dis)[colnames(comm), colnames(comm)])
  }
  else if(ncol(comm)!=attributes(dis)$Size) stop("the number of species in comm must be equal to that in dis") 		
  
  D <- as.matrix(dis)
  if(!abundance){
    comm[comm>0] <- 1
  }
  commt <- as.data.frame(t(comm))
  
  funK <- function(v){
    p <- v/sum(v)
    K <- apply(D, 1, function(x) sum(x*p))
    K[p<tol] <- NA
    return(K)
  }
  V <- cbind.data.frame(sapply(commt, funK))
  rownames(V) <- colnames(comm)
  colnames(V) <- rownames(comm)
  funKbar <- function(v){
    p <- v/sum(v)
    Kbar <- sapply(1:nrow(D), function(i) sum(D[i,]*v/sum(v[-i])))
    Kbar[p<tol] <- NA
    return(Kbar)
  }
  Kbar <- cbind.data.frame(sapply(commt, funKbar))
  rownames(Kbar) <- colnames(comm)
  colnames(Kbar) <- rownames(comm)
  funQ <- function(v){
    p <- v/sum(v)
    Q <- t(p)%*%D%*%p
    return(Q)
  }
  Q <- unlist(sapply(commt, funQ))
  
  funSim <- function(v){
    p <- v/sum(v)
    S <- 1-sum(p^2)
    return(S)
  }
  Sim <- unlist(sapply(commt, funSim))
  
  funN <- function(v){
    p <- v/sum(v)
    N <- length(p[p>tol])
    return(N)
  }
  N <- unlist(sapply(commt, funN))
  U <- Q/Sim
  
  
  red <- cbind.data.frame(N=N, D=Sim, Q=Q, U=U)
  rownames(red) <- rownames(comm)
  
  res <- list()
  res$Kbar <- Kbar
  res$V <- V
  res$red <- red
  return(res)
  
}


# Load in R the first data set analyzed in the main text. Name Com the data frame with species' abundances (species as rows and plots as columns as in Appendix S2). Name T the data frame with species' traits (species as rows and C, S, R as columns as in Appendix S2). 
Com <- species.abund # Species are columns and plots are rows
traits <- trait.levels
Com <- Com[, rownames(traits)]

# Here we define a function to calculate the Marczewski–Steinhaus coefficient among species' trait profiles.

dist.MS <- function (comm, diag = FALSE, upper = FALSE, tol = 1e-07) 
{
  
  df <- data.frame(comm)
  if (!inherits(df, "data.frame")) 
    stop("df is not a data.frame")
  nlig <- nrow(df)
  d <- matrix(0, nlig, nlig)
  d.names <- row.names(df)
  fun1 <- function(x) {
    sum(abs(df[x[1], ] - df[x[2], ]))/sum(apply(cbind.data.frame(df[x[1], ], df[x[2], ]), 1, max))
  }
  df <- as.matrix(df)
  index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
  
  d <- unlist(apply(index, 1, fun1))
  
  attr(d, "Size") <- nlig
  attr(d, "Labels") <- d.names
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- "Marczewski–Steinhaus"
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
}
# This function is a modification of function dist.quant from library ade4 (Dray & Dufour 2007) where other dissimilarity coefficients can be found.

Dis <- dist.MS(traits)
Uni <- uniqueness(Com, Dis, abundance = TRUE)

fac <- factor(rep(rep(c("Control", "Collapsed"), each = 5), 6), levels = c("Control", "Collapsed"))

#Mean for species richness N, the Rao index Q, the Simpson index D, and community-level functional uniqueness   for the three successional stages identified along the primary succession on the foreland of the Rutor glacier.
sapply(Uni$red, function(x) tapply(x, fac, mean))

#SD for species richness N, the Rao index Q, the Simpson index D, and community-level functional uniqueness   for the three successional stages identified along the primary succession on the foreland of the Rutor glacier.
sapply(Uni$red, function(x) tapply(x, fac, sd))



