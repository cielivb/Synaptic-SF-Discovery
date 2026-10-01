# Compute GO enrichment per GMS cluster, then make GO graphs illustrating GO IDs 
# that appear to characterise each cluster and illustrate relationships between 
# GO IDs for ease of functional interpretation

outdir <- "results/5_describe_gms_cluster_func"
log_file <- "logs/5_describe_gms_cluster_func.log"
source("scripts/setup.R")
configure_logger(log_file)
dir.create(glue("{outdir}/cluster_dfs"), showWarnings = FALSE)
plan(multisession, workers = parallel::detectCores(logical = FALSE))




### Statistically evaluate characteristic GO IDs per cluster --------------

log_info("Running GO ID enrichment tests ...")


# Read in and merge dataframes
cluster_df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
)
go_presence_df <- read_csv(
  "results/3_extract_go_module_scores/go_module_scores.csv",
  show_col_types = FALSE
)
df <- inner_join(cluster_df, go_presence_df, by = c("uniprot", "taxon")) |>
  select(gms_cluster, starts_with("GO:"))


# Find enriched GO IDs per GMS cluster
clusters <- unique(sort(df$gms_cluster))
results_df <- bind_rows(lapply(clusters, function (cluster_id) {
  
  cluster_df <- df |> filter(gms_cluster == cluster_id)
  bg_df <- df |> filter(gms_cluster != cluster_id)
  
  # For each GO ID, does the GO ID occur significantly more often in the cluster
  # compared to background? Get the unadjusted p-value and effect size.
  go_ids <- df |> select(-gms_cluster) |> names()
  results <- bind_rows(future_lapply(go_ids, function (go_id, cluster_df, bg_df) {
    
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
      prop_in_cluster = a / (a + b),
      prop_in_bg = c / (c + d),
      p_value = ft$p.value,
      odds_ratio = ft$estimate
    ))
  }, cluster_df = cluster_df, bg_df = bg_df))
  return(results)
  
}))


# Adjust for multiple testing via Benjamini-Hochberg on per-cluster basis and
# flag statistically enriched and high odds GO IDs
results_df <- results_df |>
  group_by(gms_cluster) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(p_adj) |>
  mutate(enriched = (p_adj < 0.1), high_odds = (odds_ratio >= 2))

# Write results
write_csv(results_df, glue("{outdir}/go_enrichment_results.csv"))
summary_df <- results_df |>
  group_by(gms_cluster) |>
  summarise(
    n_enriched = sum(enriched, na.rm = TRUE),
    n_high_odds = sum(high_odds, na.rm = TRUE),
    n_both = sum(enriched & high_odds, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_df, glue("{outdir}/go_enrichment_per_cluster.csv"))




### Create GO Graph Visual Summary per GMS Cluster ------------------------

log_info("Creating GO summary graphs ...")

gms_clusters <- unique(sort(results_df$gms_cluster))
done <- future_lapply(gms_clusters, function (gms_cluster_id, outdir) {
  cluster_df <- results_df |> filter(gms_cluster == gms_cluster_id)
  write_csv(cluster_df, glue("{outdir}/cluster_dfs/{gms_cluster_id}.csv"))
  system(glue("python scripts/make_go_graph_plot.py {gms_cluster_id}"))
}, outdir = outdir)

log_info("Done!")

