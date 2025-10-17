# iTOL toolkit
rm(list = ls())
library(itol.toolkit)
library(ape)
library(dplyr)
# Input the tree files
tree.bac <- "E:/thermokarst_gully/data/metagenome/MAGs/gtdb_tree/bacteria/tax.unrooted.tree"
tree.archa <- "E:/thermokarst_gully/data/metagenome/MAGs/gtdb_tree/archaea/tax.unrooted.tree"
# Read the taxonomic classification files
tax_all_mags <- read.table("E:/thermokarst_gully/data/metagenome/MAGs/annotation.txt", 
                        sep = "\t", header = TRUE)
mags.abun <- read.table("E:/thermokarst_gully/data/metagenome/MAGs/bin_abundance_table.tab", 
                        sep = "\t", header = TRUE) %>%
  mutate(MRA = rowMeans(.[, -1], na.rm = TRUE)) %>%
  select(c("ID", "MRA"))
annotations_table <- tax_all_mags %>%
  left_join(mags.abun, by = c("ID" = "ID"))

id_phylum <- as_tibble(annotations_table[, c("ID", "Phylum")])

# Set the tree color with the “range” attribute
itol.file.path <- "E:/thermokarst_gully/data/metagenome/MAGs/itol"
unit_7_bac <- create_unit(data = id_phylum, 
                      key = "E007_tree_colors_1_bac", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      tree = tree.bac)
unit_7_archa <- create_unit(data = id_phylum, 
                      key = "E007_tree_colors_1_archa", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      tree = tree.archa)
write_unit(unit_7_bac, file = itol.file.path)
write_unit(unit_7_archa, file = itol.file.path)

# Set the tree color with the “clade” attribute
unit_8_bac <- create_unit(data = id_phylum, 
                      key = "E008_tree_colors_2_bac", 
                      type = "TREE_COLORS", 
                      subtype = "clade", 
                      line_type = "normal",
                      color = "table2itol", 
                      size_factor = 3, 
                      tree = tree.bac)
unit_8_archa <- create_unit(data = id_phylum, 
                          key = "E008_tree_colors_2_archa", 
                          type = "TREE_COLORS", 
                          subtype = "clade", 
                          line_type = "normal",
                          color = "table2itol", 
                          size_factor = 3, 
                          tree = tree.archa)
write_unit(unit_8_bac, file = itol.file.path)
write_unit(unit_8_archa, file = itol.file.path)

# Set the DATASET_GRADIENT of completeness
unit_26@specific_themes$heatmap$color$min <- "#e5e5e5"
unit_26@specific_themes$heatmap$color$max <- "#e95f5c"
id_completeness <- as_tibble(annotations_table[, c("ID", "completeness")])

unit_27_bac <- create_unit(data = id_completeness,
                       key = "E027_gradient_2_completeness_bac",
                       type = "DATASET_GRADIENT",
                       tree = tree.bac)
unit_27_archa <- create_unit(data = id_completeness,
                           key = "E027_gradient_2_completeness_archa",
                           type = "DATASET_GRADIENT",
                           tree = tree.archa)

write_unit(unit_27_bac, file = itol.file.path)
write_unit(unit_27_archa, file = itol.file.path)

# Set the DATASET_GRADIENT of contamination
id_contamination <- as_tibble(annotations_table[, c("ID", "contamination")])
unit_28_bac <- create_unit(data =id_contamination,
                       key = "E028_gradient_3_contamination_bac",
                       type = "DATASET_GRADIENT",
                       tree = tree.bac)
unit_28_archa <- create_unit(data = id_contamination,
                           key = "E028_gradient_3_contamination_archa",
                           type = "DATASET_GRADIENT",
                           tree = tree.archa)
write_unit(unit_28_bac, file = itol.file.path)
write_unit(unit_28_archa, file = itol.file.path)







# Draw symbols on the tree
library(itol.toolkit)
library(data.table)
library(ape)
tree <- system.file("extdata",
                    "tree_of_itol_templates.tree",
                    package = "itol.toolkit")
data("template_groups")
data("template_parameters_count")
tree <- read.tree(tree)


template_parameters_count <- data.frame(template_parameters_count)
df_data <- cbind(template_groups,as.data.frame(rowSums(template_parameters_count)))
unit_30 <- create_unit(data = df_data,
                       key = "E030_symbol_1",
                       type = "DATASET_SYMBOL",
                       position = 1,
                       tree = tree)
