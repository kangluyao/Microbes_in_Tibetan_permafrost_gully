# Required packages
install.packages(c("dplyr","tidyr","ggplot2","purrr","coin"), dependencies = TRUE)
library(dplyr); library(tidyr); library(ggplot2); library(purrr); library(coin)

# ---------- Input ----------
# comm: data.frame/matrix, rows = samples, cols = species (abundance or counts)
# meta: data.frame with columns: sample, treatment (Uncollapsed/Collapsed), block (1..6 or B1..)
# Ensure sample names of comm match meta$sample, in same order or map by name.

# Example placeholders (remove if you have real data)
comm <- t(singlem_genu)
meta <- data.frame(sample = metadata$Sample_name,
                   treatment = metadata$Group,
                   block = metadata$Gully_id) %>%
  mutate(treatment = gsub("-", "",  treatment))
# comm <- read.table("comm.tsv", header=TRUE, row.names=1, sep="\t")
# meta <- read.table("meta.tsv", header=TRUE, sep="\t")

# For this script we will:
# 1) binarize abundances to presence/absence (change if you want abundance-based)
comm_pa <- (comm > 0) * 1
samples <- rownames(comm_pa)
n <- nrow(comm_pa)

# map sample -> index
sample_idx <- setNames(1:n, samples)

# ---------- Helper: compute a,b,c for unordered pair i<j ----------
pair_abc <- function(i_idx, j_idx, mat) {
  xi <- mat[i_idx, ]
  xj <- mat[j_idx, ]
  a <- sum(xi & xj)
  b <- sum(xi & (!xj))  # species in i only
  c <- sum((!xi) & xj)  # species in j only
  total <- a + b + c
  if (total == 0) {
    jacc <- NA
    lossc <- NA
    gainc <- NA
  } else {
    jacc <- (b + c) / total
    lossc <- b / total
    gainc <- c / total
  }
  data.frame(
    s1 = rownames(mat)[i_idx],
    s2 = rownames(mat)[j_idx],
    a = a, b = b, c = c, tot = total,
    jacc = jacc,
    loss_contrib = lossc,
    gain_contrib = gainc,
    stringsAsFactors = FALSE
  )
}

# ---------- compute all unique unordered pairs (upper triangle) ----------
pairs_list <- list()
k <- 1
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    pairs_list[[k]] <- pair_abc(i, j, comm_pa)
    k <- k + 1
  }
}
pairs_df <- bind_rows(pairs_list)

# attach meta info for each sample in pair
pairs_df <- pairs_df %>%
  left_join(meta %>% rename(s1 = sample, treatment1 = treatment, block1 = block), by = "s1") %>%
  left_join(meta %>% rename(s2 = sample, treatment2 = treatment, block2 = block), by = "s2")

# ---------- focus on within-block pairs only ----------
within_block_pairs <- pairs_df %>%
  filter(block1 == block2) %>%
  mutate(block = block1)

# ---------- within each block: split into within-treatment pairs ----------
# That is: pairs where both are Uncollapsed OR both are Collapsed
within_block_pairs <- within_block_pairs %>%
  mutate(within_group = case_when(
    treatment1 == treatment2 ~ treatment1,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(within_group))

# compute per-block & per-treatment means
block_treatment_summary <- within_block_pairs %>%
  group_by(block, within_group) %>%
  summarise(
    n_pairs = n(),
    mean_jacc = mean(jacc, na.rm=TRUE),
    mean_loss = mean(loss_contrib, na.rm=TRUE),
    mean_gain = mean(gain_contrib, na.rm=TRUE),
    .groups = "drop"
  )

# reshape so each block has two rows (Uncollapsed, Collapsed)
block_wide <- block_treatment_summary %>%
  pivot_wider(names_from = within_group,
              values_from = c(n_pairs, mean_jacc, mean_loss, mean_gain),
              names_sep = "_")

# ensure block ordering
block_wide <- block_wide %>% arrange(block)

# ---------- compute per-block deltas: Collapsed minus Uncollapsed ----------
block_wide <- block_wide %>%
  mutate(
    delta_jacc = mean_jacc_Collapsed - mean_jacc_Uncollapsed,
    delta_loss = mean_loss_Collapsed - mean_loss_Uncollapsed,
    delta_gain = mean_gain_Collapsed - mean_gain_Uncollapsed
  )

# ---------- classify each block into the four categories ----------
block_wide <- block_wide %>%
  mutate(
    category = case_when(
      delta_jacc > 0 & delta_loss > 0 & delta_gain <= 0 ~ "loss-driven heterogenization",
      delta_jacc > 0 & delta_gain > 0 & delta_loss <= 0 ~ "gain-driven heterogenization",
      delta_jacc > 0 & delta_loss > 0 & delta_gain > 0 ~ "mixed heterogenization",
      delta_jacc < 0 & delta_loss < 0 & delta_gain >= 0 ~ "loss-driven homogenization",
      delta_jacc < 0 & delta_gain < 0 & delta_loss >= 0 ~ "gain-driven homogenization",
      delta_jacc < 0 & delta_loss < 0 & delta_gain < 0 ~ "mixed homogenization",
      TRUE ~ "ambiguous/mixed"
    )
  )

# ---------- paired statistical tests across blocks ----------
# We perform paired Wilcoxon tests across blocks: collapse vs uncollapsed
# Create block-level vectors
delta_jacc_vec <- block_wide$delta_jacc
delta_loss_vec <- block_wide$delta_loss
delta_gain_vec <- block_wide$delta_gain

# Paired tests (Collapsed - Uncollapsed)
# Wilcoxon signed-rank (nonparametric) across the 6 blocks
test_jacc <- wilcox.test(block_wide$mean_jacc_Collapsed, block_wide$mean_jacc_Uncollapsed, paired = TRUE, exact=FALSE)
test_loss <- wilcox.test(block_wide$mean_loss_Collapsed, block_wide$mean_loss_Uncollapsed, paired = TRUE, exact=FALSE)
test_gain <- wilcox.test(block_wide$mean_gain_Collapsed, block_wide$mean_gain_Uncollapsed, paired = TRUE, exact=FALSE)

# ---------- permutation test (constrained within-block) ----------
# Permutation strategy: within each block, shuffle treatment labels among the samples,
# recompute block summaries and deltas, repeat many times to build null distribution.
permute_block_deltas <- function(comm_pa, meta, nperm = 1000, seed = 42) {
  set.seed(seed)
  blocks <- unique(meta$block)
  perm_results <- matrix(NA, nrow = nperm, ncol = 3) # delta_jacc, delta_loss, delta_gain averaged across blocks
  for (p in 1:nperm) {
    # permute treatment within each block
    meta_perm <- meta
    for (b in blocks) {
      idx <- which(meta$block == b)
      meta_perm$treatment[idx] <- sample(meta$treatment[idx])
    }
    # recompute pairs and block summaries (simpler approach: reuse pairwise results but re-evaluate within_group)
    pairs_perm <- within_block_pairs %>%
      select(s1,s2,a,b,c,tot,jacc,loss_contrib,gain_contrib,block) %>%
      left_join(meta_perm %>% rename(s1 = sample, treatment1 = treatment), by="s1") %>%
      left_join(meta_perm %>% rename(s2 = sample, treatment2 = treatment), by="s2") %>%
      filter(treatment1 == treatment2) %>%
      group_by(block, treatment1) %>%
      summarise(
        mean_jacc = mean(jacc, na.rm=TRUE),
        mean_loss = mean(loss_contrib, na.rm=TRUE),
        mean_gain = mean(gain_contrib, na.rm=TRUE),
        .groups="drop"
      ) %>%
      pivot_wider(names_from = treatment1, values_from = c(mean_jacc, mean_loss, mean_gain))
    # compute deltas per block and average
    if (nrow(pairs_perm) == 0) {
      perm_results[p, ] <- c(NA, NA, NA)
    } else {
      # ensure same ordering of blocks as block_wide
      perm_wide <- pairs_perm %>% arrange(block)
      dj <- mean(perm_wide[[paste0("mean_jacc_Collapsed")]] - perm_wide[[paste0("mean_jacc_Uncollapsed")]], na.rm=TRUE)
      dl <- mean(perm_wide[[paste0("mean_loss_Collapsed")]] - perm_wide[[paste0("mean_loss_Uncollapsed")]], na.rm=TRUE)
      dg <- mean(perm_wide[[paste0("mean_gain_Collapsed")]] - perm_wide[[paste0("mean_gain_Uncollapsed")]], na.rm=TRUE)
      perm_results[p, ] <- c(dj, dl, dg)
    }
  }
  colnames(perm_results) <- c("delta_jacc", "delta_loss", "delta_gain")
  as.data.frame(perm_results)
}

# Run permutation (nperm can be increased, here set small for speed)
perm_out <- permute_block_deltas(comm_pa, meta, nperm = 999)

# empirical mean deltas
emp_mean_delta <- c(mean(block_wide$delta_jacc, na.rm=TRUE), mean(block_wide$delta_loss, na.rm=TRUE), mean(block_wide$delta_gain, na.rm=TRUE))
names(emp_mean_delta) <- c("delta_jacc","delta_loss","delta_gain")

# p-values: proportion of permuted deltas more extreme than empirical (two-sided)
perm_pvals <- sapply(1:3, function(i) {
  col <- perm_out[[i]]
  mean(abs(col) >= abs(emp_mean_delta[i]), na.rm=TRUE)
})
names(perm_pvals) <- names(emp_mean_delta)

# ---------- outputs ----------
list(
  block_summary = block_wide,
  paired_tests = list(jacc = test_jacc, loss = test_loss, gain = test_gain),
  permutation = list(emp_mean_delta = emp_mean_delta, perm_pvals = perm_pvals)
)

# ---------- simple plots ----------
# barplot of per-block deltas
ggplot(block_wide, aes(x = factor(block), y = delta_jacc)) +
  geom_bar(stat="identity") + geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Block (Gully)", y="Delta Jaccard (Collapsed - Uncollapsed)", title="Per-block Δβ") +
  theme_minimal()

# heatmap-like visualization for delta_loss and delta_gain
block_wide_long <- block_wide %>%
  select(block, delta_jacc, delta_loss, delta_gain) %>%
  pivot_longer(cols = starts_with("delta"), names_to = "metric", values_to = "value")

ggplot(block_wide_long, aes(x = factor(block), y = metric, fill = value)) +
  geom_tile() + scale_fill_gradient2(low="blue", mid="white", high="red") +
  labs(x="Block", y="", fill="Delta value")





# Required packages
install.packages(c("vegan","betapart","dplyr","tidyr","ggplot2"))
library(vegan)
library(betapart)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------- Input ----------
# comm: abundance data (samples × species)
# meta: metadata with columns sample, treatment (Uncollapsed/Collapsed), block

# Example:
# comm <- read.table("comm_abund.tsv", header=TRUE, row.names=1, sep="\t")
# meta <- read.table("meta.tsv", header=TRUE, sep="\t")

# Ensure sample order matches metadata
comm <- comm[match(meta$sample, rownames(comm)), ]

# ---------- Step 1. Compute Bray–Curtis partition ----------
beta_abund <- beta.pair.abund(comm, index.family="bray")

# beta_abund is a list with:
# $beta.bray.bal  (balanced changes)
# $beta.bray.gra  (abundance gradient)
# $beta.bray      (total Bray-Curtis)

# Convert to matrix → long format
to_long <- function(mat, metric) {
  as.data.frame(as.table(as.matrix(mat))) %>%
    rename(s1 = Var1, s2 = Var2, value = Freq) %>%
    mutate(metric = metric)
}

beta_long <- bind_rows(
  to_long(beta_abund$beta.bray.bal, "bal"),
  to_long(beta_abund$beta.bray.gra, "gra"),
  to_long(beta_abund$beta.bray, "total")
) %>% filter(s1 != s2)

# Attach meta info
beta_long <- beta_long %>%
  left_join(meta %>% rename(s1 = sample, trt1 = treatment, blk1 = block), by="s1") %>%
  left_join(meta %>% rename(s2 = sample, trt2 = treatment, blk2 = block), by="s2") %>%
  filter(blk1 == blk2) %>% mutate(block = blk1)

# ---------- Step 2. Within-block pair summaries ----------
within_block <- beta_long %>%
  filter(trt1 == trt2) %>%
  group_by(block, trt1, metric) %>%
  summarise(mean_val = mean(value, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = trt1, values_from = mean_val) %>%
  mutate(delta = Collapsed - Uncollapsed)

# ---------- Step 3. Spread out ----------
beta_summary <- within_block %>%
  pivot_wider(names_from = metric, values_from = delta,
              names_prefix="delta_") %>%
  mutate(
    category = case_when(
      delta_total > 0 & delta_gra > 0 ~ "Gain-driven heterogenization",
      delta_total > 0 & delta_gra < 0 ~ "Loss-driven heterogenization",
      delta_total < 0 & delta_gra > 0 ~ "Gain-driven homogenization",
      delta_total < 0 & delta_gra < 0 ~ "Loss-driven homogenization",
      TRUE ~ "Ambiguous"
    )
  )

# ---------- Step 4. Visualize ----------
# (a) stacked Δ components
beta_long_plot <- beta_summary %>%
  pivot_longer(cols = starts_with("delta_"), names_to="component", values_to="value")

ggplot(beta_long_plot, aes(x=factor(block), y=value, fill=component)) +
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Block (Gully)", y="Δ Value (Collapsed - Uncollapsed)",
       title="Decomposition of β-diversity change (abundance-based)") +
  scale_fill_manual(values=c("grey40","skyblue","tomato"),
                    labels=c("Balanced (replacement)","Gradient (gain/loss)","Total Bray-Curtis")) +
  theme_minimal()

# (b) classification summary
ggplot(beta_summary, aes(x=factor(block), y=delta_total, fill=category)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Block", y="Δβ (Collapsed - Uncollapsed)", title="β-diversity change classification") +
  scale_fill_brewer(palette="Set2") +
  theme_minimal()

# ---------- Output ----------
beta_summary
