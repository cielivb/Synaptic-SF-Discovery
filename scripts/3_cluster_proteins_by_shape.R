# Cluster proteins by shape.

source("scripts/setup.R")
configure_logger("logs/3_cluster_proteins_by_shape.log")

library(GGally)
library(ggplot2)
library(tidyr)
library(uwot)
library(dbscan)




plot_features_species_violin <- function(df, title) {
  
  # Rotate df to get df with taxon, feature, and value columns
  geom_long <- df |>
    pivot_longer(cols = -taxon, names_to = "feature", values_to = "value")
  
  # Create violin plot
  plt <- ggplot(geom_long, aes(x = taxon, y = value, fill = taxon)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point", size = 1.5, color = "black") +
    facet_grid(rows = vars(feature), scales = "free_y") +
    labs(title = title) +
    theme(strip.text.y = element_text(angle = 0))
  
  # Export violin plot
  outfile <- glue(
    "results/3_cluster_proteins_by_shape/feature_selection/{clean_string(title)}.svg"
  )
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  svg(outfile, width = 12, height = 18)
  print(plt)
  invisible(dev.off())
  
}




dir.create(
  "results/3_cluster_proteins_by_shape", showWarnings = FALSE, recursive = TRUE
)

# Read in geometry descriptions for each taxon
taxa <- c("fly","human","mouse","rat","worm","zebrafish")
geom_df <- bind_rows(lapply(taxa, function (taxon) {
  geom_file <- glue("results/2_extract_shape_data/{taxon}_geom.csv")
  taxon_df <- read_csv(geom_file, show_col_types = FALSE) |>
    mutate(taxon = taxon)
  return(taxon_df)
})) |> select(-uniprot)


# Visualise feature distributions before corrections
geom_df_f <- geom_df |> select(taxon, ends_with("_f"))
geom_df_hc <- geom_df |> select(taxon, ends_with("_hc"))
geom_df_misc <- geom_df |> select(-ends_with("_f"), -ends_with("_hc"))
plot_features_species_violin(geom_df_f, "Feature Set 1 Distributions by Taxa (Raw)")
plot_features_species_violin(geom_df_hc, "Feature Set 2 Distributions by Taxa (Raw)")
plot_features_species_violin(geom_df_misc, "Feature Set 3 Distributions by Taxa (Raw)")

# Correct for skewness in candidate features within taxa, then visualise.
# kurtosis_f, length_f, rg_f, kurtosis_hc, length_hc, rg_hc
no_skew_f <- geom_df_f |> 
  mutate(log_kurtosis_f = log(kurtosis_f),
         log_length_f = log(length_f),
         log_rg_f = log(rg_f)) |>
  select(-kurtosis_f, -length_f, -rg_f)
no_skew_hc <- geom_df_hc |> 
  mutate(log_kurtosis_hc = log(kurtosis_hc),
         log_length_hc = log(length_hc),
         log_rg_hc = log(rg_hc)) |>
  select(-kurtosis_hc, -length_hc, -rg_hc)
no_skew_misc <- geom_df_misc
plot_features_species_violin(
  no_skew_f, "Feature Set 1 Distributions by Taxa (Skew-Corrected)")
plot_features_species_violin(
  no_skew_hc, "Feature Set 2 Distributions by Taxa (Skew-Corrected)")
plot_features_species_violin(
  no_skew_misc, "Feature Set 3 Distributions by Taxa (Skew-Corrected)")

# Scale feature distributions within taxa, then visualise
# Use MAD and median to stay robust against outliers

mad_scale <- function(v) {
  (v - median(v)) / mad(v)
}

scaled_f <- no_skew_f |> group_by(taxon) |>
  mutate( # TODO
    log_kurtosis_f = mad_scale(log_kurtosis_f)
  ) |> ungroup()
scaled_hc <- no_skew_hc |> group_by(taxon) |>
  mutate( # TODO
    ) |> ungroup()
scaled_misc <- no_skew_misc
plot_features_species_violin(
  scaled_f, "Feature Set 1 Distributions by Taxa (Skew-Corrected, MAD-Scaled)")
plot_features_species_violin(
  scaled_hc, "Feature Set 2 Distributions by Taxa (Skew-Corrected, MAD-Scaled")
plot_features_species_violin(
  scaled_misc, "Feature Set 3 Distributions by Taxa (Skew-Corrected, MAD-Scaled")


# Run per-taxon PCA on descriptor matrices per species to find variances per 
# feature, then filter out features with consistently low variance across species
# TODO


# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))

# Find within-taxa shape clusters
future_lapply(taxa, future.seed = TRUE, FUN = function (taxon) {
  
  # Create UMAP embedding
  set.seed(254)
  features <-  # TODO
  umap_embedding <- umap(
    features,
    n_neighbours = 50,
    min_dist = 0.1,
    n_components = 2,
    metric = "euclidean",
    init = "pca"
  )
  features$UMAP1 <- umap_embedding[,1]
  features$UMAP2 <- umap_embedding[,2]
  
  # Run HDBSCAN on UMAP embedding to discover clusters
  hdb <- hdbscan(umap_embedding, minPts = 10)

})




