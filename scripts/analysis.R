# Influence of Synapse Type on Receptor Protein Geometry in Homo Sapiens
# Author: Ciel Baumann
# Last updated: NA

source("scripts/setup.R")
source("scripts/geometry.R")
log_info("CONGRATS! YOU GOT THIS FAR!")

CLUSTER <- makeCluster(detectCores())

csv_files <- list.files("data_csv", pattern = "\\.csv$", full.names = TRUE)
geoms <- parLapplyLB(CLUSTER, csv_files, build_geom_tibble)

stopCluster(CLUSTER)


# 
# ## Categorise Proteins by GO Annotations
# Add Inhibitory Synapse Receptor and Excitatory Synapse Receptor column to
# Geometry Descriptor Table, by iteratively getting and filtering GO annotations
# per protein.
# 
# ## Visualise Geometric Distributions
# 
# ## Get Statistical Significance and Effect Sizes
# Do a nonparametric test to determine whether excitatory synapse receptor proteins,
# inhibitory synapse receptor proteins, and other proteins have different geometric
# distributions. I am not sure yet whether to use Mann-Whitney U or Kruskal-Wallis
# or something else. I could use Cliff's Delta for effect size.
# 
# ## Create Post-Hoc Visualisations if Applicable