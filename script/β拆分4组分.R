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
  DB <- ecopart_multi_baselga_sp_single(d1_block, d2_block)  # 4 x S
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
# 4) Visualize (heatmap + top species per component)
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
ggplot(top20_data, aes(x = interaction(genera, rank(mean_contribution)), y = mean_contribution)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_contribution - se_contribution, 
                    ymax = mean_contribution + se_contribution),
                width = 0.3, linewidth = 0.5) +
  geom_point(size = 2.5, color = "#2C7BB6") +
  facet_wrap(~ component, scales = "free") +
  scale_x_discrete(labels = function(x) str_remove(x, "\\.[^.]*$")) +
  coord_flip() +
  labs(x = "Genera", 
       y = "Contribution (mean ± SE)",
       title = "Top 20 Genera by Average Contribution in Each Component") +
  main_theme +
  theme(strip.text = element_text(face = "bold", size = 6),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p)
