# Cluster protein shapes within taxa, then detect similar clusters across taxa.

source("scripts/setup.R")
configure_logger("logs/4_cluster_shapes.log")

library(GGally)
library(uwot)
library(dbscan)




### CLUSTER WITHIN TAXA --------------------------------------------------------


# Plot UMAP embedding, colouring points by clusters
plot_umap_embedding_clusters <- function(geom_df, taxon) {
  
  plt <- ggplot(geom_df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_discrete(name = "Cluster") +
    theme_minimal(base_size = 14) +
    labs(
      title = glue("UMAP Embedding of Protein Geometry ({taxon})"),
      x = "UMAP1", y = "UMAP2"
    )
  
  umap_plot_path <- glue("results/4_cluster_shapes/umap_{taxon}.svg")
  svg(umap_plot_path, width = 8, height = 6)
  print(plt)
  invisible(dev.off())
  
}


# Read in dataframe
geom_df <- read_csv(
  "results/3_select_shape_features/selected.csv", show_col_types = FALSE)

# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))


# Find within-taxa shape clusters
taxa <- c("fly","human","mouse","rat","worm","zebrafish")
clustered <- future_lapply(taxa, future.seed = TRUE, FUN = function (taxon) {
  
  taxon_geom <- geom_df |> filter(taxon == taxon)
  
  # Create UMAP embedding
  set.seed(254)
  features <- taxon_geom |> select(starts_with("scaled_"))
  umap_embedding <- umap(
    features,
    n_neighbours = 50,
    min_dist = 0.1,
    n_components = 2,
    metric = "euclidean",
    init = "pca"
  )
  taxon_geom$UMAP1 <- umap_embedding[,1]
  taxon_geom$UMAP2 <- umap_embedding[,2]
  
  # Run HDBSCAN on UMAP embedding to discover clusters
  hdb <- hdbscan(umap_embedding, minPts = 10)
  taxon_geom$cluster <- hdb$cluster # cluster 0 = noise
  taxon_geom$cluster_prob <- hdb$membership_prob
  
  # Visualise clusters on UMAP embedding
  plot_umap_embedding_clusters(taxon_geom, taxon)
  
  return(taxon_geom)
  
})
cluster_df <- bind_rows(clustered)




### FIND GEOMETRICALLY SIMILAR WITHIN-TAXA-CLUSTERS ACROSS TAXA ----------------



