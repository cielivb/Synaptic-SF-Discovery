# Make plots describing relationships between taxa and GMS clusters, including:
# - Tile plot illustrating GMS cluster presence/absence per taxa
# - Stacked bar plot illustrating proportion of proteins in each taxa assigned 
#   to each GMS cluster, with taxa side-by-side, and bars stacked per taxa
# - One bar plot per GMS cluster illustrating proportion of proteins in each taxa
#   assigned to that cluster
# - One stacked single-column bar per GMS cluster illustrating proportion of 
#   proteins in that cluster belonging to each taxa (colour by taxa)
# - Grid plots for the above two plot types for ease of cluster comparison

outdir <- "results/6_make_taxa_plots"
log_file <- "logs/6_make_taxa_plots.log"
source("scripts/setup.R")
configure_logger(log_file)

library(ggplot2)


# Read in dataframe
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv", show_col_types = FALSE
)


### Create tile plot illustrating GMS cluster presence/absence per taxa ---

# Create tibble: taxon | module | present
