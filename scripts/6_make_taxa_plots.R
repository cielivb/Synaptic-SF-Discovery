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

# Read in dataframe and phylogenetic tree
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv", show_col_types = FALSE
)
tree <- read.tree("data/taxa_tree.nwk")
p <- ggtree(tree) + theme_tree2()


### Create tile plot illustrating GMS cluster presence/absence per taxa ---

# Create P/A tibble: taxon | 0 | ... | N (N GMS clusters, values = 1 or 0)
pa_df <- df |> select(taxon, gms_cluster) |>
  mutate(present = 1) |> # gms cluster is present in taxon
  distinct() |>
  pivot_wider(
    names_from = gms_cluster,
    values_from = present,
    values_fill = 0, # gms cluster absent in taxon
    names_prefix = "Cluster "
  )

# Create tile plot with attached phylogenetic tree
svg(glue("{outdir}/cluster_pa_by_taxa_tile_plot.svg"), height = 10, width = 18)
plt <- gheatmap(p, pa_df, offset = 0.02, width = 0.8, colnames = TRUE) +
  scale_fill_manual(values = c("0" = "white", "1" = "#2ECC71")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(plt)
invisible(dev.off())