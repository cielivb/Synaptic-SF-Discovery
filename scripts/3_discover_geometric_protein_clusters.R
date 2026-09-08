# Cluster proteins by geometry and identify similar clusters across taxa


# Create a violin plot where each violin represents a feature and the panels
# represent geometric protein clusters within an organism
plot_cluster_feature_distr_within_organism <- function(organism) {
  return() # TODO
}


# Plot UMAP embedding, colouring points by clusters
plot_umap_embedding_clusters <- function(geom_df, organism) {
  
  plt <- ggplot(geom_df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_discrete(name = "Cluster") +
    theme_minimal(base_size = 14) +
    labs(
      title = glue("UMAP Embedding of Protein Geometry ({organism})"),
      x = "UMAP1",
      y = "UMAP2"
    )
  
  umap_plot_path <- glue("results/clusters/umap/{organism}_umap.svg")
  svg(umap_plot_path, width = 8, height = 6)
  print(plt)
  dev.off()
}


# Assign proteins to geometric clusters within an organism
cluster_proteins_by_geometry <- function(organism) {
  
  # Load in geometry dataframe tibble and remove all rows with NA values.
  # NA values are most likely to be present in highly disordered proteins 
  # without stable cores. UMAP requires all features to contain no NA values.
  geom_path <- SPECIES$geom_path[match(organism, SPECIES$common_name)]
  geom_df <- read_csv(geom_path, show_col_types = FALSE) |> drop_na()
  
  # Get UMAP embeddings
  set.seed(524)
  features <- geom_df |> select(-uniprot) |> scale()
  umap_embedding <- umap(
    features,
    n_neighbors = 50,
    min_dist = 0.1,
    n_components = 2,
    metric = "euclidean",
    init = "pca"
  )
  geom_df$UMAP1 <- umap_embedding[,1]
  geom_df$UMAP2 <- umap_embedding[,2]
  
  # Run HDBSCAN on UMAP embedding columns to discover clusters
  organism_min_pts = list(
    human=30, fruit_fly=10, mouse=30, rat=25, worm=10, zebrafish=15
  )
  hdb <- hdbscan(umap_embedding, minPts = organism_min_pts[[organism]])
  geom_df$cluster <- hdb$cluster # cluster 0 = noise
  geom_df$cluster_prob <- hdb$membership_prob
  
  # Compute and write cluster level summaries
  cluster_summary_df <- geom_df |> 
    group_by(cluster) |> 
    summarise(
      across(
        high_pLDDT:kurtosis_ratio,
        list(
          mean = mean,
          min = min,
          max = max,
          mad = mad
        ),
        .names = "{.col}_{.fn}"
      )
    )
  write_csv(cluster_summary_df, glue("results/clusters/{organism}_cluster_summary.csv"))
  
  # Create UMAP embedding plot
  plot_umap_embedding_clusters(geom_df, organism)
  
  # Visualise feature distributions of each cluster side-by-side
  plot_cluster_feature_distr_within_organism(organism)
  
  # Write dataframe with UMAP embeddings and cluster IDs
  outpath <- SPECIES$clusters_path[match(organism, SPECIES$common_name)]
  write_csv(geom_df, outpath)
  log_info("Clustered {organism} proteins by geometry")
  
}
