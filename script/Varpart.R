library(vegan)

siteXspe <- data.frame(t(mags_abun_tab))
env <- metadata[, -c(1:4)]
traits <- cwm_results
mag_tab <- decostand(siteXspe, 'hellinger')

# Bray-Curtis dissimilarity
dist_bc <- vegdist(mag_tab, method = "bray")

# Full db-RDA models for each predictor group
env_full <- dbrda(dist_bc ~ ., data = env)
trait_full <- dbrda(dist_bc ~ ., data = traits)

# Check multicollinearity
vif.cca(env_full)
vif.cca(trait_full)

# Refine the full models by removing 
# predictors with high VIF values (e.g., VIF > 10)
env_refine <- env[, -which(vif.cca(env_full) > 10)]
trait_refine <- traits[, -which(vif.cca(trait_full) > 10)]

env_full_refine <- dbrda(dist_bc ~ ., 
                         data = env_refine)
trait_full_refine <- dbrda(dist_bc ~ .,
                           data = trait_refine)

# Null models
env_null <- dbrda(dist_bc ~ 1, data = env_refine)
trait_null <- dbrda(dist_bc ~ 1, data = trait_refine)

# Forward selection with adjusted R² constraint
env_sel <- ordiR2step(
  env_null,
  scope = formula(env_full_refine),
  R2scope = TRUE,
  permutations = 999
)

trait_sel <- ordiR2step(
  trait_null,
  scope = formula(trait_full_refine),
  R2scope = TRUE,
  permutations = 999
)

# Select the final factors for VPA
env_final <- env[, c("pH", "NO3_N", "Soil_moisture", 
                     "MAP", "Slope", "AGB", "Plant_richness")]
trait_final <- traits[, c(
  "Resource Acquisition:Substrate uptake:carbohydrate transport",
  "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration",
  "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation",
  "Stress Tolerance:Specific:envelope stress",
  "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization",
  "Resource Acquisition:Substrate assimilation:N compounds",
  "Stress Tolerance:Specific:oxidative stress",
  "Stress Tolerance:Specific:desiccation/osmotic/salt stress",
  "Stress Tolerance:Specific:temperature:low")]

# Variation partitioning
# Replace env_final and trait_final with
# the predictor matrices retained by forward selection.
vp <- varpart(
  dist_bc,
  env_final,
  trait_final
)

plot(vp)
