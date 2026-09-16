# (Attempt to) Cluster protein shapes within taxa.

# "Do synaptic proteins across taxa form discrete geometric clusters in 
# geometric feature space?"

# Results: no meaningful geometric clusters could be found. The manifold appears
# continous across species, and many proteins consistently fall into cluster 0
# ("noise"), and the majority of proteins have lower membership probabilities
# in their assigned clusters. Silhouette analysis (for taxa with more than 1
# assigned cluster) also yields average silhouette widths between 0.05 and 0.24
# indicating poor cluster assignment. These results were taken from using 
# minPts = 5 in both hdbscan calls. It seems unlikely that adjusting parameters
# of either UMAP or HDBSCAN will alter these results. Therefore, I will conclude
# that I could not detect clear and meaningful cluster structure in geometric 
# feature space.

# Interesting tangents:
# - Would this be true for other protein sets beyond synaptic proteins?


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
  
  umap_plot_path <- glue("results/4_cluster_shapes/umap_{taxon}.svg")
  svg(umap_plot_path, width = 8, height = 6)
  
  plt <- ggplot(geom_df, aes(x = UMAP1, y = UMAP2, color = factor(umap_cluster))) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_discrete(name = "Cluster") +
    theme_minimal(base_size = 14) +
    labs(
      title = glue("UMAP Embedding of Protein Geometry ({taxon})"),
      x = "UMAP1", y = "UMAP2"
    )
  
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
  while (dev.cur() > 1) dev.off()
  set.seed(254)
  taxon_geom <- geom_df |> filter(taxon == !!taxon)
  log_info(glue("Clustering {nrow(taxon_geom)} {taxon} proteins ..."))
  
  
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
  taxon_geom$cluster <- as.factor(hdb$cluster) # cluster 0 = noise
  taxon_geom$cluster_prob <- hdb$membership_prob
  
  # Visualise HDBSCAN-only clusters on scaled features
  outfile <- glue("results/4_cluster_shapes/hdbscan_clusters_on_features_{taxon}.svg")
  svg(outfile, width = 12, height = 12)
  plt <- ggpairs(
    taxon_geom, columns = colnames(features), aes(color = cluster)) +
    labs(title = glue("Feature Distributions with HDBSCAN Clusters ({taxon})"))
  print(plt)
  invisible(dev.off())
  
  # Visualise HDBSCAN-only cluster probabilities
  outfile <- glue("results/4_cluster_shapes/hdbscan_cluster_probs_{taxon}.svg")
  svg(outfile, width = 12, height = 12)
  plt <- ggplot(taxon_geom, aes(cluster_prob, fill = cluster)) +
    geom_histogram(binwidth = 0.05) +
    labs(title = glue("HDBSCAN Cluster Probabilities ({taxon})"))
  print(plt)
  invisible(dev.off())
  
  # Get and visualise HDBSCAN-only cluster silhouette scores
  if (length(unique(as.integer(taxon_geom$cluster))) > 1) {
    distance_matrix <- dist(features, method = "euclidean")
    sink(glue("results/4_cluster_shapes/silhouette_summary_{taxon}.txt"))
    sil <- silhouette(as.integer(taxon_geom$cluster), distance_matrix)
    print(summary(sil))
    sink()
    outfile <- glue("results/4_cluster_shapes/hdbscan_silhouette_{taxon}.svg")
    svg(outfile, width = 12, height = 12)
    plot(sil, col = taxon_geom$cluster)
    invisible(dev.off())
  } else {
    log_info(glue("Skipping silhouette for {taxon}: only 1 cluster"))
  }
  
  # Project HDBSCAN-only clusters onto UMAP for visualisation
  umap_df <- as.data.frame(umap_embedding)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  proj_df <- bind_cols(umap_df, taxon_geom["cluster"])
  outfile <- glue("results/4_cluster_shapes/hdbscan_umap_proj_{taxon}.svg")
  svg(outfile, width = 12, height = 12)
  plt <- ggplot(proj_df, aes(UMAP1, UMAP2, color = factor(cluster))) +
    geom_point(size = 1.5, alpha = 0.8) +
    # scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    labs(
      title = glue("UMAP-Projected HDBSCAN-Only Clusters ({taxon})"),
      color = "Cluster"
    )
  print(plt)
  invisible(dev.off())
  
  log_info(glue("Computed {taxon} clusters"))
  while (dev.cur() > 1) dev.off()
  return(taxon_geom)
  
})
cluster_df <- bind_rows(clustered)





### FIND GEOMETRICALLY SIMILAR WITHIN-TAXA-CLUSTERS ACROSS TAXA ----------------



