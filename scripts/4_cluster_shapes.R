# Cluster protein shapes within taxa, then detect similar clusters across taxa.

source("scripts/setup.R")
log_file <- "logs/4_cluster_shapes.log"
configure_logger(log_file)

library(cluster)
library(GGally)
library(uwot)
library(dbscan)


dir.create("results/4_cluster_shapes", showWarnings = FALSE)



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
  "results/3_select_shape_features/selected.csv", 
  col_names = TRUE, show_col_types = FALSE)

# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))


# Find within-taxa shape clusters
log_info("Computing within-taxa shape clusters ...")
taxa <- c("fly","human","mouse","rat","worm","zebrafish")
clustered <- future_lapply(taxa, future.seed = TRUE, FUN = function (taxon) {
  
  configure_logger(log_file, worker = TRUE)
  taxon_geom <- geom_df |> filter(taxon == !!taxon)
  log_info(glue("Clustering {nrow(taxon_geom)} {taxon} proteins ..."))
  set.seed(254)
  
  # Create UMAP embedding
  features <- taxon_geom |> select(starts_with("scaled_"))
  umap_embedding <- umap(
    features,
    n_neighbors = 30,
    min_dist = 0.2,
    n_components = 2,
    metric = "euclidean",
    init = "pca"
  )
  taxon_geom$UMAP1 <- umap_embedding[,1]
  taxon_geom$UMAP2 <- umap_embedding[,2]
  
  # Run HDBSCAN on UMAP embedding to discover UMAPped clusters
  hdb <- hdbscan(umap_embedding, minPts = 5)
  taxon_geom$umap_cluster <- hdb$cluster # cluster 0 = noise
  taxon_geom$umap_cluster_prob <- hdb$membership_prob
  
  # Visualise clusters on UMAP embedding
  plot_umap_embedding_clusters(taxon_geom, taxon)
  
  # Run HDBSCAN on the scaled features themselves (skipping UMAP - HDBSCAN-only)
  hdb <- hdbscan(features, minPts = 5)
  taxon_geom$cluster <- hdb$cluster # cluster 0 = noise
  taxon_geom$cluster_prob <- hdb$membership_prob
  
  # Visualise HDBSCAN-only clusters on scaled features
  plt <- ggpairs(
    taxon_geom, columns = colnames(features), aes(color = cluster)) +
    labs(title = glue("Feature Distributions with HDBSCAN Clusters ({taxon})"))
  outfile <- glue("results/4_cluster_shapes/features_clusters_{taxon}.svg")
  svg(outfile, width = 12, height = 12)
  print(plt)
  invisible(dev.off())
  
  # Visualise HDBSCAN-only cluster probabilities
  plt <- ggplot(taxon_geom, aes(cluster_prob, fill = cluster)) +
    geom_histogram(binwidth = 0.05) +
    labs(title = glue("HDBSCAN Cluster Probabilities ({taxon})"))
  
  # Get and visualise HDBSCAN-only cluster silhouette scores
  distance_matrix <- dist(features, method = "euclidean")
  sil <- silhouette(taxon_geom$cluster, distance_matrix)
  sink(glue("results/4_cluster_shapes/silhouette_summary_{taxon}.txt"))
  summary(sil)
  sink()
  plt <- plot(sil, col = taxon_geom$cluster)
  outfile <- glue("results/4_cluster_shapes/silhouette_{taxon}.svg")
  svg(outfile, width = 12, height = 12)
  print(plt)
  invisible(dev.off())
  
  # Project HDBSCAN-only clusters onto UMAP for visualisation
  umap_df <- as.data.frame(umap_embedding$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  proj_df <- bind_cols(umap_df, taxon_geom["cluster"])
  plt <- ggplot(proj_df, aes(UMAP1, UMAP2, color = factor(cluster))) +
    geom_point(size = 1.5, alpha = 0.8) +
   # scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    labs(
      title = glue("UMAP-Projected HDBSCAN-Only Clusters ({taxon})")
      color = "Cluster"
    )
  
  log_info(glue("Computed {taxon} clusters"))
  return(taxon_geom)
  
})
cluster_df <- bind_rows(clustered)





### FIND GEOMETRICALLY SIMILAR WITHIN-TAXA-CLUSTERS ACROSS TAXA ----------------



