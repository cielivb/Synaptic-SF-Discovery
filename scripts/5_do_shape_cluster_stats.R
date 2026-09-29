# Identify protein shape clusters within GMS clusters, and identify shape
# clusters that are significantly and not significantly different in their
# shape distributions via PERMANOVA (non-parametric - no assumption of normality
# etc, critical for my non-normal shape data even after log transform)

outdir <- "results/5_do_shape_cluster_stats"
log_file <- "logs/5_do_shape_cluster_stats.log"
source("scripts/setup.R")
configure_logger(log_file)

library(cluster) # silhouette scoring
library(dbscan)


# Launch cluster for parallel processing
plan(multisession, workers = parallel::detectCores(logical = FALSE))

# Load data
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv", show_col_types = FALSE
)


### Identify shape clusters within each GMS cluster -----------------------

# Shape cluster keys are composite keys (GMS cluster + shape cluster (unique 
# only within a GMS cluster))

gms_cluster_dfs <- future_lapply(
  unique(df$go_module_cluster), function (gms_cluster) {
    
    # Get all proteins within GMS cluster
    gms_df <- df |> filter(go_module_cluster == gms_cluster)
    
    # Identify shape clusters within GMS cluster via pure HDBSCAN, using minPts
    # = 10 for bare minimum statistical robustness later
    shape_features <- gms_df |> select(starts_with("scaled_"))
    hdb <- hdbscan(shape_features, minPts = 10)
    gms_df$shape_cluster <- hdb$cluster # cluster 0 = noise
    gms_df$shape_cluster_membership_prob <- hdb$membership_prob
    
    # Visualise shape cluster membership probability distribution
    svg(
      glue("{outdir}/gms_{gms_cluster}_hdbscan_probs.svg"), 
      width = 8, height = 8
    )
    plt <- ggplot(gms_df, aes(shape_cluster_membership_prob)) +
      geom_histogram(binwidth = 0.05) +
      labs(
        title = "HDBSCAN Shape Cluster Probabilities",
        subtitle = glue("GMS Cluster {gms_cluster}")
      )
    print(plt)
    invisible(dev.off())
    
    # Visualise silhouette scores
    num_shape_clusters <- length(unique(as.integer(gms_df$shape_cluster)))
    if (num_shape_clusters > 1) {
      
      # Create the silhouette scoring plot
      distance_matrix <- dist(shape_features, method = "euclidean")
      sil <- silhouette(as.integer(gms_df$shape_cluster), distance_matrix)
      svg(
        glue("{outdir}/gms_{gms_cluster}_silhouette.svg"), 
        width = 8, height = 8
      )
      plot(sil, col = gms_df$shape_cluster)
      invisible(dev.off())
      
    } else if (num_shpae_clusters == 1) {
      log_info(glue(
        "Skipping shape cluster silhouette scoring for {gms_cluster}:" /
        "only one cluster found"
      ))
    } else {
      log_info(glue(
        "Skipping shape cluster silhouette scoring for {gms_cluster}:" /
        "no clusters found"
      ))
    }
    
  }, outdir = outdir
)
clustered_df <- bind_rows(gms_cluster_dfs) 








