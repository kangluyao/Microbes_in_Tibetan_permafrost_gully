# Microbial profiles in thermo-erosion gullies on the Tibetan plateau

The study investigates the microbial communities in thermo-erosion gullies on the Tibetan plateau, focusing on their diversity, composition, and potential ecological functions.

## Code structure

The analysis scripts are located in the `script/` directory. Script `0.read_data.R` loads the shared input data and is sourced by the numbered analysis scripts, which are intended to be run from the project root (`E:/thermokarst_gully`).

```text
script/
├── 0.read_data.R                 # Shared data loader: reads 16S OTU/taxonomy/metadata, rarefied OTU table,
│                                 # MAG abundance/taxonomy/attribute tables; builds
│                                 # phyloseq objects used by all downstream scripts.
├── 1_SingleM_data_analysis.R     # SingleM community analysis (Collapsed vs. Un-collapsed):
│                                 # phylum/genus composition, beta-regression GLMMs of differential abundance,
│                                 # shared-genera Venn diagram, alpha-diversity GLMMs, PERMANOVA/PCoA and a
│                                 # combined multi-panel figure.
├── 2_Singlem_β_partition.R       # SingleM beta diversity: beta dispersion (betadisper with within-gully
│                                 # permutations) and ecopart partitioning of Δβ into additive/subtractive and
│                                 # differentiation/homogenization components, including per-genus contributions.
├── 3_MAGs_β_partition.R          # MAG-level beta diversity: Hellinger transformation, Bray-Curtis
│                                 # PERMANOVA/PCoA, beta dispersion, ecopart partitioning of Δβ and per-MAG
│                                 # contributions; produces the combined beta-diversity figure.
├── 4_Traits_analysis.R           # MAG trait ecology: classifies MAGs as increasing/decreasing in collapsed
│                                 # gullies, compares microtrait-based traits, correlates traits, and
│                                 # tests community-weighted mean (CWM) traits with Gaussian GLMMs.
├── 5_16S_analysis.R              # 16S rRNA community analysis: phylum-level composition, alpha-diversity
│                                 # GLMMs, Bray-Curtis PERMANOVA/PCoA, beta dispersion and combined figure.
└── 6_Analyses_supp_info.R        # Supplementary analyses: MAG taxonomy and dominant genera,
                                  # and soil environmental
                                  # variables and heterogeneity, sequencing/MAG recruitment/quality statistics,
                                  # Procrustes comparison of 16S vs. MAG communities, and diagnostic exports.
```

### Typical workflow

1. Run `script/0.read_data.R` (or let the analysis scripts source it) to load all shared inputs.
2. Run scripts `1`–`6` independently; each script starts from a clean workspace and sources `0.read_data.R`.