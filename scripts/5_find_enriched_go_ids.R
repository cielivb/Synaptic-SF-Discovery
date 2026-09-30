# Find GO IDs enriched in each GMS cluster


outdir <- "results/5_find_enriched_go_ids"
log_file <- "logs/5_find_enriched_go_ids.log"
source("scripts/setup.R")
configure_logger(log_file)


# Read in and merge dataframes
cluster_df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
)
go_presence_df <- read_csv(
  "results/3_extract_functional_features/functional_features.csv"
)
df <- inner_join(cluster_df, go_presence_df, by = c("uniprot", "taxon")) |>
  select(gms_cluster, starts_with("GO:"))


# Find enriched GO IDs per GMS cluster
clusters <- unique(sort(df$gms_cluster))
results_df <- bind_rows(lapply(clusters, function (cluster_id, df) {
  
  cluster_df <- df |> filter(gms_cluster == cluster_id)
  bg_df <- df |> filter(gms_cluster != cluster_id)
  
  # For each GO ID, does the GO ID occur significantly more often in the cluster
  # compared to background? Get the unadjusted p-value and effect size.
  go_ids <- df |> select(-gms_cluster) |> names()
  results <- bind_rows(future_lapply(go_ids, function (go_id) {
    
    # Counts for contingency table
    a <- sum(cluster_df[[go_id]] == 1) # GO ID present in cluster df
    b <- sum(cluster_df[[go_id]] == 0) # GO ID absent in cluster df
    c <- sum(bg_df[[go_id]] == 1) # GO ID present in background df
    d <- sum(bg_df[[go_id]] == 0) # GO ID absent in background df
    
    # Fisher's exact test
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    return(tibble(
      gms_cluster = cluster_id,
      go_id = go_id,
      p_value = ft$p.value,
      odds_ratio = ft$estimate
    ))
  }, cluster_df = cluster_df, bg_df = bg_df))
  return(results)
  
}, df = df))


# Adjust for multiple testing via Benjamini-Hochberg
results_df <- results_df |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Write results
write_csv(results_df, glue("{outdir}/go_enrichment_results.csv"))

