# Cluster proteins by geometry and identify similar clusters across taxa


# Create a violin plot where each violin represents a feature and the panels
# represent geometric protein clusters within an organism
plot_cluster_feature_distr_within_organism <- function(organism) {
  return() # TODO
}


# Assign proteins to geometric clusters within an organism
cluster_proteins_by_geometry <- function(organism) {
  
  # Load in geometry dataframe tibble
  geom_path <- SPECIES$geom_path[match(organism, SPECIES$common_name)]
  geom_df <- read_csv(geom_path, show_col_types = FALSE)
  
  # Get UMAP embeddings
  set.seed(254)
  features <- geom_df |> select(-uniprot) |> scale()
  umap_embedding <- umap(
    features,
    n_neighbors = 15,
    min_dist = 0.1,
    n_components = 3,
    metric = "euclidean"
  )
  geom_df$UMAP1 <- umap_emb[,1]
  geom_df$UMAP2 <- umap_emb[,2]
  
  # Run HDBSCAN on UMAP embedding columns to discover clusters
  hdb <- hdbscan(umap_emb, minPts = 10)
  geom_df$cluster <- hdb$cluster # cluster 0 = noise
  geom_df$cluster_prob <- hdb$membership_prob
  
  # Compute and write cluster level summaries
  cluster_summary_df <- geom_df |> 
    group_by(cluster) |> 
    summarise(
      across(
        length:high_pLDDT,
        list(
          mean = mean,
          min = min,
          max = max,
          mad = mad
        ),
        .names = "{.col}_{.fn}"
      )
    )
  write_csv(cluster_summary_df, "results/clusters/{organism}_cluster_summary.csv")
  
  # Visualise feature distributions of each cluster side-by-side
  plot_cluster_feature_distr_within_organism(organism)
  
  # Write dataframe with UMAP embeddings and cluster IDs
  outpath <- SPECIES$clusters_path[match(organism, SPECIES$common_name)]
  write_csv(geom_df, outpath)
  log_info("Clustered {organism} proteins by geometry")
  
}
