# Identify protein shape clusters within GMS clusters, and identify shape
# clusters that are significantly and not significantly different in their
# shape distributions via PERMANOVA (non-parametric - no assumption of normality
# etc, critical for my non-normal shape data even after log transform)

outdir <- "results/5_do_shape_cluster_stats"
silhouette_dir <- glue("{outdir}/silhouettes")
log_file <- "logs/5_do_shape_cluster_stats.log"
source("scripts/setup.R")
configure_logger(log_file)


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
    shape_features <- gms_df |>
      select(
        scaled_mad_mean_pLDDT, scaled_mad_eccentricity_u,
        scaled_mad_coeff_var_u, scaled_mad_log_kurtosis_u,
        scaled_mad_log_length_u
      )
    hdb <- hdbscan(shape_features, minPts = 10)
    gms_df$shape_cluster <- hdb$cluster # cluster 0 = noise
    gms_df$shape_cluster_membership_prob <- hdb$membership_prob
    
    # Visualise shape cluster membership probability distribution
    
    # Run a silhouette analysis
  }
)
clustered_df <- bind_rows(gms_cluster_dfs) 








