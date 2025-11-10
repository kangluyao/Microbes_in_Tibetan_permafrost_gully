# ================
# Packages
# ================
if(!requireNamespace("tidyverse", quietly=TRUE)) install.packages("tidyverse")
library(tidyverse)

# ================
# Inputs (replace with your data)
# ================
# d1 : uncollapsed (30 samples x S species)  (row names = sample IDs)
# d2 : collapsed   (30 samples x S species)
# meta: data.frame with columns: sample, block, treatment
#       treatment values should indicate "Uncollapsed" and "Collapsed" (or similar)
# Ensure colnames(d1) == colnames(d2)
comm <- t(singlem_genu)
meta <- data.frame(sample = metadata$Sample_name,
                   treatment = metadata$Group,
                   block = metadata$Gully_id) %>%
  mutate(treatment = gsub("-", "",  treatment))

d1 <- comm[grepl("_C", rownames(comm)), ]
d1[1:30, 1:5]
d2 <- comm[grepl("_T", rownames(comm)), ]
d2[1:30, 1:5]
stopifnot(identical(colnames(d1), colnames(d2)))
S <- ncol(d1)
species_names <- colnames(d1)

# ================
# 1) Extract Tatsumi's baselga code as a function that returns DBeta (4 x S) for two matrices d1,d2
#    This is essentially the "baselga" branch of ecopart.multi with components="sp"
# ================
ecopart_multi_baselga_sp_single <- function(d1_block, d2_block) {
  # d1_block and d2_block: matrices with same number of rows N and same columns S
  if(ncol(d1_block) != ncol(d2_block)) stop("d1 and d2 must have same species columns")
  if(nrow(d1_block) != nrow(d2_block)) stop("d1 and d2 must have same number of rows (N) per block")
  
  N <- nrow(d1_block)
  S <- ncol(d1_block)
  DBeta <- matrix(0, nrow = 4, ncol = S)
  
  not.zero <- function(x) sapply(1:length(x), function(i) all.equal(as.numeric(x[i]), 0)) != TRUE
  A <- sum(d1_block)
  M <- sum(apply(d1_block, 2, max))
  Y <- N * M / (N - 1) / sum(d2_block)
  
  for(s in 1:S) {
    Abn1 <- as.numeric(d1_block[, s])
    Abn2 <- as.numeric(d2_block[, s])
    Max1 <- max(Abn1)
    Max2 <- max(Abn2)
    # handle edge cases where maxima matches multiple positions
    A2_M1 <- ifelse(any(Abn1 == Max1), max(Abn2[Abn1 == Max1]), 0)
    A1_M2 <- ifelse(any(Abn2 == Max2), max(Abn1[Abn2 == Max2]), 0)
    AX <- (Abn1 > A2_M1 & Abn1 > A1_M2 & Abn2 > A2_M1 & Abn2 > A1_M2)
    AX <- ifelse(any(AX), max(pmin(Abn1[AX], Abn2[AX])), 0)
    LossL <- ifelse(Abn1 > Max2, Abn1 - Max2, 0)
    LossD <- ifelse(Max1 > A2_M1 & Max2 > A1_M2 & Abn1 > A2_M1 & Abn1 > A1_M2 & Abn1 > AX & Abn1 > Abn2,
                    pmin(Max1, Max2, Abn1) - max(A2_M1, A1_M2, AX), 0)
    LossS <- ifelse(Abn1 > Abn2, Abn1 - Abn2 - LossL - LossD, 0)
    GainL <- ifelse(Abn2 > Max1, Abn2 - Max1, 0)
    GainD <- ifelse(Max1 > A2_M1 & Max2 > A1_M2 & Abn2 > A2_M1 & Abn2 > A1_M2 & Abn2 > AX & Abn2 > Abn1,
                    pmin(Max1, Max2, Abn2) - max(A2_M1, A1_M2, AX), 0)
    GainS <- ifelse(Abn2 > Abn1, Abn2 - Abn1 - GainL - GainD, 0)
    
    # create level slices
    Lvl <- sort(c(0, Abn1, Abn2))
    if(length(Lvl) < 2) {
      # degenerate case: all zeros
      DBeta[, s] <- 0
      next
    }
    Lvl.mean <- sapply(1:(2*N), function(i) sum(Lvl[i:(i+1)])/2)
    Lvl.diff <- diff(Lvl)
    
    # counts per level
    LossL.up  <- ifelse(not.zero(LossL), Abn1, 0)
    LossL.low <- ifelse(not.zero(LossL), Abn1 - LossL, 0)
    LossD.up  <- ifelse(not.zero(LossD), Abn1 - LossL, 0)
    LossD.low <- ifelse(not.zero(LossD), Abn1 - LossL - LossD, 0)
    LossS.up  <- ifelse(not.zero(LossS), Abn1 - LossL - LossD, 0)
    LossS.low <- ifelse(not.zero(LossS), Abn1 - LossL - LossD - LossS, 0)
    GainL.up  <- ifelse(not.zero(GainL), Abn2, 0)
    GainL.low <- ifelse(not.zero(GainL), Abn2 - GainL, 0)
    GainD.up  <- ifelse(not.zero(GainD), Abn2 - GainL, 0)
    GainD.low <- ifelse(not.zero(GainD), Abn2 - GainL - GainD, 0)
    GainS.up  <- ifelse(not.zero(GainS), Abn2 - GainL - GainD, 0)
    GainS.low <- ifelse(not.zero(GainS), Abn2 - GainL - GainD - GainS, 0)
    
    LossL.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < LossL.up & Lvl.mean[i] > LossL.low))
    LossD.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < LossD.up & Lvl.mean[i] > LossD.low))
    LossS.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < LossS.up & Lvl.mean[i] > LossS.low))
    GainL.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < GainL.up & Lvl.mean[i] > GainL.low))
    GainD.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < GainD.up & Lvl.mean[i] > GainD.low))
    GainS.n <- sapply(1:(2*N), function(i) sum(Lvl.mean[i] < GainS.up & Lvl.mean[i] > GainS.low))
    
    Term.1 <- ifelse(LossL.n != 0, Y * (-1/M + LossL.n/A) * Lvl.diff, 0)
    Term.2 <- ifelse(LossD.n != 0, Y * (-1/M + LossD.n/A) * Lvl.diff, 0)
    Term.3 <- ifelse(LossS.n != 0, Y * (LossS.n/A) * Lvl.diff, 0)
    Term.4 <- ifelse(GainL.n != 0, Y * (1/M - GainL.n/A) * Lvl.diff, 0)
    Term.5 <- ifelse(GainD.n != 0, Y * (1/M - GainD.n/A) * Lvl.diff, 0)
    Term.6 <- ifelse(GainS.n != 0, -Y * (GainS.n/A) * Lvl.diff, 0)
    
    DBeta[1, s] <- sum(Term.1[Term.1 < 0], Term.2[Term.2 < 0])
    DBeta[2, s] <- sum(Term.1[Term.1 > 0], Term.2[Term.2 > 0], Term.3)
    DBeta[3, s] <- sum(Term.4[Term.4 < 0], Term.5[Term.5 < 0], Term.6)
    DBeta[4, s] <- sum(Term.4[Term.4 > 0], Term.5[Term.5 > 0])
  } # end species loop
  
  rownames(DBeta) <- c("Subtractive homogenization", "Subtractive differentiation",
                       "Additive homogenization", "Additive differentiation")
  colnames(DBeta) <- colnames(d1_block)
  return(DBeta)
}

source("https://raw.githubusercontent.com/communityecologist/ecopart/refs/heads/main/R/ecopart-multi.R")
# ================
# 2) Loop over blocks to compute DBeta per block, then average across blocks
# ================
# meta must map sample -> block and which treatment (Uncollapsed/Collapsed)
# Here assume you have meta with columns: sample, block, treatment (values exactly "Uncollapsed"/"Collapsed")
blocks <- unique(meta$block)
DBeta_blocks <- list()
for(b in blocks) {
  # indices of samples for this block and treatment
  idx1 <- meta %>% filter(block == b & treatment == "Uncollapsed") %>% pull(sample)
  idx2 <- meta %>% filter(block == b & treatment == "Collapsed") %>% pull(sample)
  # check counts
  if(length(idx1) < 1 || length(idx2) < 1) {
    warning("Block ", b, " missing samples in one treatment; skipping")
    next
  }
  # extract matrices
  d1_block <- d1[idx1, , drop = FALSE]
  d2_block <- d2[idx2, , drop = FALSE]
  # require same N (here both should be 5)
  if(nrow(d1_block) != nrow(d2_block)) {
    stop("Block ", b, ": unequal sample counts per treatment (need equal N per block for this implementation).")
  }
  DB <- ecopart.multi(d1_block, d2_block, index = "baselga", components="four")  # 4 x S
  DBeta_blocks[[as.character(b)]] <- DB
}

# Now aggregate across blocks: element-wise mean across the 6 block matrices
# Stack into array
arr <- array(NA, dim = c(4, S, length(DBeta_blocks)))
i <- 1
for(bn in names(DBeta_blocks)) {
  arr[,,i] <- DBeta_blocks[[bn]]
  i <- i + 1
}
# compute mean across blocks (third dimension)
res_sp <- apply(arr, c(1,2), mean, na.rm = TRUE)  # 4 x S matrix
rownames(res_sp) <- rownames(DBeta_blocks[[1]])
colnames(res_sp) <- colnames(DBeta_blocks[[1]])

# ================
# 3) Output and tidy result for plotting
# ================
# res_sp is 4 x S species matrix you requested
res_sp  # inspect

# Tidy for plotting
res_long <- as.data.frame(res_sp) %>%
  rownames_to_column("Component") %>%
  pivot_longer(-Component, names_to = "Species", values_to = "Contribution")

# compute percent contribution within each component (absolute)
res_long <- res_long %>%
  group_by(Component) %>%
  mutate(Percent = 100 * Contribution / sum(abs(Contribution), na.rm=TRUE)) %>%
  ungroup()

# ================
# 4) Visualize (top species per component)
# ================
library(ggplot2)
# Top contributors per component (absolute contribution)
rep_sixblocks_long <- rbind(DBeta_blocks[["EB"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "EB"),
                            DBeta_blocks[["ML"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "ML"),
                            DBeta_blocks[["RS"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "RS"),
                            DBeta_blocks[["SLH"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "SLH"),
                            DBeta_blocks[["HSX"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "HSX"),
                            DBeta_blocks[["HH"]] %>% as.data.frame() %>%
                              rownames_to_column("component") %>%
                              pivot_longer(-component, names_to = "genera",
                                           values_to = "contribution") %>%
                              mutate(blocks = "HH"))

library(tidyverse)
# 步骤1：计算每个component中每个genera的平均contribution和标准误
summary_data <- rep_sixblocks_long %>%
  group_by(component, genera) %>%
  summarise(
    mean_contribution = mean(contribution, na.rm = TRUE),
    se_contribution = sd(contribution, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# 步骤2：在每个component中选出平均contribution排名前20的genera
top20_data <- summary_data %>%
  group_by(component) %>%
  slice_max(order_by = abs(mean_contribution), n = 20) %>%
  ungroup()

# 步骤3：为了更好的可视化，按平均contribution对genera排序
top20_data <- top20_data %>%
  group_by(component) %>%
  mutate(genera = fct_reorder(genera, mean_contribution)) %>%
  ungroup()

# 步骤4：绘图
top20_4process_plot <- ggplot(top20_data, aes(x = interaction(genera, rank(mean_contribution)), y = mean_contribution)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_contribution - se_contribution, 
                    ymax = mean_contribution + se_contribution),
                width = 0.3, linewidth = 0.25) +
  geom_point(size = 1, color = "#2C7BB6", shape = 16) +
  facet_wrap(~ component, scales = "free") +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
  coord_flip() +
  labs(x = "Genera", 
       y = "Contribution (mean ± SE)") +
  main_theme +
  theme(strip.text = element_text(face = "bold", size = 6),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

top20_4process_plot


# ggsave(filename = "E:/thermokarst_gully/revision/result/top20_4process_plot.pdf",
#        top20_4process_plot, width = 89, height = 89, units = "mm", dpi = 600)

# Compare the difference in genera traits
bac_gtdb <- fread("E:/thermokarst_gully/data/metagenome/singlem/bac120_metadata_r226.tsv") %>%
  select(gtdb_taxonomy, ncbi_species_taxid, gc_percentage, 
         genome_size, ncbi_rrna_count, checkm2_completeness)

arc_gtdb <- fread("E:/thermokarst_gully/data/metagenome/singlem/ar53_metadata_r226.tsv") %>%
  select(gtdb_taxonomy, ncbi_species_taxid, gc_percentage, 
         genome_size, ncbi_rrna_count, checkm2_completeness)

gtdb_trait <- rbind(bac_gtdb, arc_gtdb) %>%
  mutate(
    gtdb_taxonomy = str_extract(gtdb_taxonomy, "g__[^;]+"),
  ) %>%
  mutate(across(c(gtdb_taxonomy),  ~ str_replace_all(., ".*__", ""))) %>%
  mutate(
    # Replace any non-numeric entries (e.g., text like 'NA' or 'none') with R's NA
    ncbi_rrna_count = na_if(ncbi_rrna_count, "none"), # Replace specific text
    ncbi_rrna_count = as.numeric(ncbi_rrna_count)
  ) %>% # correct the genome size and rrn copies using completeness
  mutate(genome_size = genome_size/(checkm2_completeness/100),
         ncbi_rrna_count = ncbi_rrna_count/(checkm2_completeness/100))

MAGs_genera <- unique(tax_bin$Genus) %>% gsub("g__", "", .)
singlem_genera <- unique(top20_data$genera)
intersect(MAGs_genera, singlem_genera)
setdiff(MAGs_genera, singlem_genera)
setdiff(singlem_genera, MAGs_genera)



gtdb_trait_top20_genera_ave <- gtdb_trait %>%
  rename("genera" = "gtdb_taxonomy") %>%
  group_by(genera) %>%
  summarise(across(c(gc_percentage, genome_size, 
                     ncbi_rrna_count), mean, na.rm = TRUE)) %>%
  inner_join(top20_data, by = "genera") %>%
  mutate(genome_size = genome_size/1000000)


head(gtdb_trait_top20_genera_ave)


traits_id <- c("gc_percentage", "genome_size", "ncbi_rrna_count")
compo_comparisons <- list(c("Subtractive differentiation", "Subtractive homogenization"),
                          c("Subtractive homogenization", "Additive differentiation"),
                          c("Additive differentiation", "Additive homogenization"),
                          c("Subtractive differentiation", "Additive differentiation"),
                          c("Subtractive homogenization", "Additive homogenization"),
                          c("Subtractive differentiation", "Additive homogenization"))
singlem_traits_compo_comparison <- gtdb_trait_top20_genera_ave %>% 
  select(c("component", all_of(traits_id))) %>%
  pivot_longer(-component, names_to = "traits", values_to = "value") %>%
  mutate(component = factor(component, 
                            levels = c("Subtractive differentiation",
                                       "Subtractive homogenization",
                                       "Additive differentiation", 
                                       "Additive homogenization"))) %>%
  mutate(traits = factor(traits, levels = traits_id)) %>%
  ggplot(aes(component, value)) +
  geom_half_violin(aes(fill = component), color = NA, 
                   position = position_nudge(x = 0.25), 
                   side = "r", width = 0.5, alpha = 0.65) +
  # geom_boxplot(aes(fill = component), width = 0.35, 
  #              size = 0.3, outlier.color = NA, alpha = 0.65,) +
  geom_jitter(aes(colour = component), shape = 16, size = 0.5,
              width = 0.15, alpha = 0.65) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.35, linewidth = 0.15) +
  stat_compare_means(comparisons = compo_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", 
                     bracket.size = 0.25, bracket.width = 0.1,  size = 2,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~traits, scales = "free_y", ncol = 6) + #
  main_theme +
  theme(legend.position = "none",
        axis.ticks.length=unit(.1, "cm"),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
singlem_traits_compo_comparison

library(patchwork)
combined_beita_par_plot <- disp_plot + 
  top20_4process_plot +
  plot_layout(ncol = 1, heights = c(1, 2.5))
ggsave(filename = "E:/thermokarst_gully/revision/result/combined_beita_par_plot.pdf",
       combined_beita_par_plot, width = 165, height = 150, units = "mm", dpi = 600)


top20_4process_traits_plot <- plot_grid(top20_4process_plot, 
                                        singlem_traits_compo_comparison,
                                        ncol = 1, rel_heights = c(2,1))
# ggsave(filename = "E:/thermokarst_gully/revision/result/top20_4process_traits_plot.pdf",
#        top20_4process_traits_plot, width = 150, height = 150, units = "mm", dpi = 600)


genera_traits_bar_plot <- gtdb_trait_top20_genera_ave %>% 
  select(c("genera", all_of(traits_id))) %>%
  pivot_longer(-genera, names_to = "traits", values_to = "value") %>%
  group_by(genera, traits) %>%
  summarise(value = mean(value, na.rm =T), .groups = "drop") %>%
  mutate(traits = factor(traits, levels = traits_id)) %>%
  ggplot(aes(genera, value)) +
  geom_bar(stat = "identity", aes(fill = "#5cc3e8"), color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_wrap(~traits, scales = "free_y", ncol = 6) + #
  main_theme +
  theme(legend.position = "none",
        axis.ticks.length=unit(.1, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))
genera_traits_bar_plot

# ggsave(filename = "E:/thermokarst_gully/revision/result/genera_traits_bar_plot.pdf",
#        genera_traits_bar_plot, width = 89, height = 89, units = "mm", dpi = 600)




plot_correlation(
  data = gtdb_trait_top20_genera_ave %>%
    filter(component == "Subtractive homogenization"),
  x_var = "genome_size",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = gtdb_trait_top20_genera_ave %>%
    filter(component == "Subtractive differentiation"),
  x_var = "gc_percentage",
  y_var = "mean_contribution",
  method = "spearman")

plot_correlation(
  data = gtdb_trait_top20_genera_ave %>%
    filter(component == "Subtractive differentiation"),
  x_var = "ncbi_rrna_count",
  y_var = "mean_contribution",
  method = "spearman")
