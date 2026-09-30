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
enriched <- future_lapply(clusters, function (cluster_id, df) {
  cluster_df <- df |> filter(gms_cluster == cluster_id)
  
  # For each GO ID, does the GO ID occur significantly more often in the cluster
  # compared to background?
  
}, df = df)