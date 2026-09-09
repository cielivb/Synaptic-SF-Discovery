# Identify clusters that are not statistically significantly different
# based on length and disorder content, as well as the full feature set.

# Future enhancements could include pairwise computation of the proportions of 
# proteins in a cluster that have an ortholog in other clusters via OrthoDB.


# Load in and clean organism's geometry descriptor + clusters dataframe
load_clusters_df <- function(organism) {
  clusters_df_path <- SPECIES$clusters_path[match(organism, SPECIES$common_name)]
  clusters_df <- read_csv(clusters_df_path, show_col_types = FALSE) |>
    select(-uniprot, -UMAP1, -UMAP2, -cluster_prob) |>
    filter(cluster != 0) |> # exclude noise "cluster" 0
    mutate(
      species = organism,
      label = as.factor(glue("{organism}_{cluster}"))
    )
}


# Calculate raw permanova results for each cluster combination in dataframe
# using all available features. No multiple testing adjustment is done here.
pairwise_permanova <- function(df) {
  
  # Get list of every cluster pair
  labels <- sort(unique(df$label))
  pair_list <- combn(labels, 2, simplify = FALSE)
  
  # Perform PERMANOVA on every cluster pair
  clusterExport(CLUSTER, "df")
  results <- parLapplyLB(CLUSTER, pair_list, function(pair) {
    
    df_pair <- df |> filter(label %in% pair)
    
    permanova_pair <- tibble(adonis2(
      df_pair |> select(-species, -cluster, -label) ~ label,
      data = df, method = "euclidean"
    ))
    label_row <- permanova_pair |> filter(Term == "label")
    
    return(list(
      cluster_1 = pair[1],
      cluster_2 = pair[2],
      Df = label_row$Df,
      SumOfSqs = label_row$SumOfSqs,
      R2 = label_row$R2,
      F = label_row$F,
      PrF = label_row$`Pr(>F)`
    ))
    
  })
  
  clusterEvalQ(CLUSTER, rm(df)) # Free worker memory
  return(bind_rows(results))

}


# "Cluster" similar clusters by finding connected components of a graph
# where nodes are cluster labels and edges represent being similar (i.e.,
# not being statistically different / alpha >= 0.05)
cluster_similar_clusters <- function(name, pw_permanova_results) {
  alpha <- 0.05
  
  # Create adjacency list where edges are similar clusters
  similar_pairs <- pw_permanova_results |> 
    filter(pval_adj >= alpha) |> select(cluster_1, cluster_2)
  
  # Find connected components ("groups" / clusters of clusters)
  graph <- graph_from_data_frame(similar_pairs, directed = FALSE)
  groups <- components(graph)$membership
  groups_tbl <- tibble(cluster = names(groups), group = groups)
  
  # Plot similarity graph
  svg("results/clusters/similarity/{name}_similarity_graph.svg")
  plot(g_sub, vertex.size = 30, vertex.label.cex = 0.8)
  dev.off()
  
  return(groups_tbl)
}


discover_similar_clusters_across_species <- function() {
  
  # Read in each species' geometry + cluster dataframes and scale features
  full_df <- bind_rows(lapply(SPECIES$common_name, load_clusters_df)) |>
    mutate(
      across(
        .cols = where(is.numeric) & !cluster,
        .fns = scale
      )
    )
  
  # Run overall PERMANOVA on full feature set
  #full_permanova <- tibble(adonis2(
  #  full_df |> select(-species, -cluster, -label) ~ species / label,
  #  data = full_df, method = "euclidean"
  #))
  #write_csv(full_permanova, "results/clusters/similarity/permanova_full.csv")
  
  # Run overall PERMANOVA on feature subset
  #subset_df <- full_df |>
  #  select(species, label, length_f, length_hc, high_pLDDT, low_pLDDT)
  #subset_permanova <- tibble(adonis2(
  #  subset_df |> select(-species, -label) ~ species / label,
  #  data = subset_df, method = "euclidean"
  #))
  #write_csv(subset_permanova, "results/clusters/similarity/permanova_subset.csv")
  
  # Do pairwise PERMANOVA, adjust for multiple testing, and write outputs
  raw_pw_permanova_results_full <- pairwise_permanova(full_df)
  raw_pw_permanova_results_sub <- pairwise_permanova(subset_df)
  pw_permanova_results_full <- raw_pw_permanova_results_full |>
    mutate(pval_adj = p.adjust(`Pr(>F)`, method = "BH"))
  pw_permanova_results_sub <- raw_pw_permanova_results_sub |>
    mutate(pval_adj = p.adjust(`Pr(>F)`, method = "BH"))
  write_csv(
    pw_permanova_results_full, 
    "results/clusters/similarity/pw_permanova_full.csv"
  )
  write_csv(
    pw_permanova_results_sub, 
    "results/clusters/similarity/pw_permanova_subset.csv"
  )
  
  # Identify similar clusters (those that are not statistically significantly 
  # different after multiple testing)
  sim_clusters_full <- cluster_similar_clusters("full", pw_permanova_results_full)
  sim_clusters_sub <- cluster_similar_clusters("subset", pw_permanova_results_sub)
  write_csv(sim_clusters_full, "results/clusters/similarity/full_similar_clusters.csv")
  write_csv(sim_clusters_sub, "results/clusters/similarity/subset_similar_clusters.csv")
  
  # Visualise results with panelled radar plots
  
  log_info("Discovered similar clusters across species")
}