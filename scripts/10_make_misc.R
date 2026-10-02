# This script is for creating miscellaneous/supplementary tables and plots that
# have not already been made

outdir <- "results/10_make_misc"
log_file <- "logs/10_make_misc.log"
source("scripts/setup.R")
configure_logger(log_file)

library(png)
library(tidyr)
library(ggplot2)

# Create stacked histogram of proteins per cluster, coloured by taxon
cluster_counts <- read_csv(
  "results/4_cluster_by_gms/cluster_counts.csv",
  show_col_types = FALSE
) |> pull(n_proteins)
svg(glue("{outdir}/proteins_per_gms_cluster.svg"), height = 8, width = 8)
hist(
  cluster_counts, 
  breaks = 5,
  xlim = c(min(cluster_counts), max(cluster_counts)),
  xlab = "Number of Proteins", ylab = "Number of GMS Clusters",
  main = "Distribution of Protein Counts Per GMS Cluster",
)
invisible(dev.off())