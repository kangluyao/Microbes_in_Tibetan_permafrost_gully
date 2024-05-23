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

spe_abun_test <- load("species.abund.rda")
species.abund[1:5, 1:5]
trait_test <- load("trait.levels.rda")
trait.levels[1:5, 1:5]

fr_raw <- read.table(file = "E:/thermokarst_gully/data/metagenome/humann2/ko.tsv", header = T, row.names = 1, sep = "\t")
fr_dat <- fr_raw %>%
  mutate(KO = sapply(strsplit(rownames(.),"[|]"),'[',1),
         Species = sapply(strsplit(rownames(.),"[.]"),'[',2)) %>%
  filter(!KO %in% c("UNMAPPED", "UNGROUPED") & !Species %in% c(NA, "unclassified"))
  
fr_dat[1:5, 60:62]

nrow(fr_dat)


species.abund <- fr_dat %>% select(-c("KO")) %>%
  group_by(Species) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var = "Species")

species.abund[1:5, 1:5]

nrow(species.abund)


trait.levels <- fr_dat %>% pivot_longer(cols = -c("KO", "Species"), names_to = "Sample", values_to = "Value") %>%
  select(-c("Sample")) %>%
  group_by(Species, KO) %>%
  summarise(across(everything(), sum)) %>%
  pivot_wider(names_from = KO, values_from = Value) %>%
  column_to_rownames(var = "Species")

trait.levels[1:5, 1:5]
nrow(trait.levels)
