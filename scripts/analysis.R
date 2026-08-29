# Influence of Synapse Type on Receptor Protein Geometry in Homo Sapiens
# Author: Ciel Baumann
# Last updated: NA

source("setup.R")



# ## Build Geometry Descriptor Table
# Useful basic metrics include:
#   - Length, Width, Height : derived from PCA projections
# - Radius of Gyration (Rg) : compactness / RMS distance from centroid
# - Eccentricity : elongation
# - Sphericity
# - Coefficient of Variation : spread
# - Kurtosis : tailedness
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