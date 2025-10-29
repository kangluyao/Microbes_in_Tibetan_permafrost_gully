# --- Load required libraries ---
install.packages(c("betapart", "vegan", "permute", "lme4", 
                   "reshape2", "ggplot2", "indicspecies"), 
                 dependencies = TRUE)
library(betapart)
library(vegan)
library(permute)
library(lme4)
library(reshape2)
library(ggplot2)
library(indicspecies)

# --- Input data assumptions ---
# comm: community matrix (samples as rows, taxa as columns), containing abundance or relative abundance
# treatment: factor vector, length = nrow(comm), levels = "Control", "Treat"
# block: factor vector, length = nrow(comm), levels = "B1"..."B6"
# Example:
comm <- t(singlem_genu)
treatment <- metadata$Group
block <- metadata$Gully_id

# --- 1) Beta-diversity decomposition ---
# For abundance data (recommended)
bpair_abund <- beta.pair.abund(comm, index.family = "bray")
# Returned object includes:
# bpair_abund$beta.bray       : total Bray–Curtis dissimilarity
# bpair_abund$beta.bray.bal   : balanced variation in abundance (turnover-like)
# bpair_abund$beta.bray.gra   : abundance gradient (nestedness-like)

# Convert to dist objects for betadisper()
dist_total <- as.dist(bpair_abund$beta.bray)
dist_turn  <- as.dist(bpair_abund$beta.bray.bal)
dist_nest  <- as.dist(bpair_abund$beta.bray.gra)

# --- 2) Calculate dispersion (distance to group centroid) for each treatment ---
bd_total <- betadisper(dist_total, treatment)
bd_turn  <- betadisper(dist_turn, treatment)
bd_nest  <- betadisper(dist_nest, treatment)

# --- 3) Permutation test with block-constrained randomization ---
# Create permutation control object
# Permutations are restricted within each block (so that treatment labels are shuffled only within blocks)
ctrl_how <- how(nperm = 999, blocks = block, within = Within(type = "free"))

# Run permutational tests
pt_total <- permutest(bd_total, permutations = ctrl_how)
pt_turn  <- permutest(bd_turn, permutations = ctrl_how)
pt_nest  <- permutest(bd_nest, permutations = ctrl_how)

# Inspect results
pt_total
pt_turn
pt_nest

# --- 4) Visualization: boxplots of sample dispersion values ---
df_disp <- data.frame(
  sample = names(bd_total$distances),
  disp_total = bd_total$distances,
  disp_turn  = bd_turn$distances,
  disp_nest  = bd_nest$distances,
  treatment = treatment,
  block = block
)

# Plot dispersion (divergence) distributions
p1 <- ggplot(df_disp, aes(x = treatment, y = disp_total)) + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("Total dispersion")
p2 <- ggplot(df_disp, aes(x = treatment, y = disp_turn))  + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("Turnover dispersion")
p3 <- ggplot(df_disp, aes(x = treatment, y = disp_nest))  + geom_boxplot() + geom_jitter(width = 0.1) + ggtitle("Nestedness dispersion")
print(p1); print(p2); print(p3)

# --- 5) Alternative test: Linear Mixed Model (LMM) on within-block mean pairwise distances ---
# Compute pairwise distances among samples within the same block and treatment
long_turn <- melt(as.matrix(dist_turn), varnames = c("s1", "s2"), value.name = "turn")
long_total <- melt(as.matrix(dist_total), varnames = c("s1", "s2"), value.name = "total")
# Keep only upper triangle (no duplicates)
long_turn <- subset(long_turn, as.character(s1) < as.character(s2))
long_total <- subset(long_total, as.character(s1) < as.character(s2))

# Map metadata (treatment and block) to sample IDs
meta <- data.frame(sample = rownames(comm), treatment = treatment, block = block, stringsAsFactors = FALSE)
long_turn <- merge(long_turn, meta, by.x = "s1", by.y = "sample", all.x = TRUE)
long_turn <- merge(long_turn, meta, by.x = "s2", by.y = "sample", suffixes = c(".1", ".2"))

# Select only pairs within the same block and treatment
long_turn_within <- subset(long_turn, block.1 == block.2 & treatment.1 == treatment.2)

# Calculate mean pairwise distance per block × treatment
mean_within <- aggregate(turn ~ block.1 + treatment.1, data = long_turn_within, FUN = mean)
colnames(mean_within) <- c("block", "treatment", "mean_turn")

# Fit linear mixed model: treatment as fixed effect, block as random effect
lmm <- lmer(mean_turn ~ treatment + (1 | block), data = mean_within)
summary(lmm)
anova(lmm)

# --- 6) Identify taxa contributing to divergence ---
# SIMPER analysis: taxa contributing most to between-group dissimilarity
sim <- simper(comm, group = treatment, permutations = 0)
summary(sim)

# Indicator species analysis: taxa significantly associated with each treatment
comm_pa <- (comm > 0) * 1
ind <- multipatt(comm_pa, treatment, control = how(nperm = 999))
summary(ind)

