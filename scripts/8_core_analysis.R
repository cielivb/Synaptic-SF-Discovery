# The core analysis is applied both across taxa only. Within taxa is not planned
# as many clusters will be of an inadequate sample size within taxa.

# Want to know if the distributions of geometric shapes of proteins in cluster
# A are statistically different from those in cluster B -> use raw/absolute 
# geometric features (except length and kurtosis, which will be log transformed
# due to high skew - still biologically meaningful)

# Ho: The full multivariate distribution of protein shape features in cluster X
# is identical to the full multivariate distribution (location, scale, shape, 
# tails, modality) of proteins shape features in cluster Y.

# A tile plot illustrating significance of difference between pairwise cluster 
# comparisons is also planned. It will be an nxn triangle tile plot where n is 
# the number of GMS clusters. Intensity of hue should match significance. 
# Ideally, tiles matching the same cluster should be greyed out.




outdir <- "results/8_core_analysis"
log_file <- "logs/8_core_analysis.log"
source("scripts/setup.R")
configure_logger(log_file)
plan(multisession, workers = parallel::detectCores(logical = FALSE))

library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(gghalves)
library(energy) # statistical testing
library(tidyr)

set.seed(254)




### Prepare data ----------------------------------------------------------

# Read in protein GMS cluster assignments and taxa
cluster_df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
) |> select(taxon, uniprot, gms_cluster)

# Read in geometric features (non-normalised)
geom_df <- read_csv(
  "results/2_extract_shape_features/shape_features.csv",
  show_col_types = FALSE
) |> rename(
  length = length_u, kurtosis = kurtosis_u, 
  coeff_var = coeff_var_u, eccentricity = eccentricity_u
) |>
  mutate(log_length = log(length), log_kurtosis = log(kurtosis)
) |>
  select(uniprot, mean_pLDDT, log_length, log_kurtosis, coeff_var, eccentricity) 

# Merge dataframes - taxon, uniprot, gms_cluster, mean_pLDDT, log_length, 
# log_kurtosis, coeff_var, eccentricity
df <- inner_join(cluster_df, geom_df, by = "uniprot")
write_csv(df, glue("{outdir}/core_df.csv"))




### Define analysis functions ---------------------------------------------

# Define null result dataframes (use when sample size too small or 
# to avoid duplicate analyses)
get_na_energy_result <- function(clust1, clust2) {
  tibble(clust1 = clust1, clust2 = clust2, statistic = NA, p = NA)
}

# Run pairwise energy tests between two dataframes of shape features
pairwise_energy <- function(clust1, clust2, df1, df2, feature_cols) {
  set.seed(254)
  X <- as.matrix(df1[, feature_cols])
  Y <- as.matrix(df2[, feature_cols])
  result <- eqdist.etest(rbind(X, Y), sizes = c(nrow(X), nrow(Y)), R = 199)
  return(tibble(
    clust1 = clust1,
    clust2 = clust2,
    statistic = result$statistic,
    p = result$p.value
  ))
}


# Take a results dataframe containing some NA results and adjust for multiple
# testing via Benjamini Hochberg, omitting all "results" with NA values from
# the BH total. NA values represent analyses not run either because the pair
# contained the same cluster (clustX vs clustX), or sample size was insufficient
# in at least one cluster
adjust_pvalues <- function(results_df) {

  # Separate NA and non-NA results and add p_adj column
  na_results <- results_df |> filter(is.na(p)) |> mutate(p_adj = NA)
  non_na_results <- results_df |> filter(!is.na(p)) |> 
    mutate(p_adj = p.adjust(p, method = "BH"))
  
  # Reattach NA results and return
  adjusted_df <- bind_rows(na_results, non_na_results)
  return(adjusted_df)
}




### Run analysis ----------------------------------------------------------

log_info("Running pairwise energy tests ...")

clusters <- unique(sort(df$gms_cluster))
min_size <- 10 # Require at least 10 proteins in a cluster to do a comparison

# Get unadjusted pairwise energy test results
energy_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(
    clusters,
    function (clust2, clust1, df, min_size) {
      
      # Skip duplicate tests or tests with tiny clusters
      if (clust1 <= clust2)
        return(get_na_energy_result(clust1, clust2))
      df1 <- df |> filter(gms_cluster == clust1)
      df2 <- df |> filter(gms_cluster == clust2)
      if (nrow(df1) < min_size || nrow(df2) < min_size) {
        return(get_na_energy_result(clust1, clust2))
      }
      
      # Define shape feature sets
      all_features <- c(
        "mean_pLDDT", "log_length", "log_kurtosis", "coeff_var", "eccentricity"
      )
      core_features <- c("mean_pLDDT", "log_length")
      irregularity_features <- c("log_kurtosis", "coeff_var", "eccentricity")
      
      # Get energy results for each shape feature set
      all_result <- pairwise_energy(clust1, clust2, df1, df2, all_features) |>
        mutate(features = "All Shape\nFeatures")
      core_result <- pairwise_energy(clust1, clust2, df1, df2, core_features) |>
        mutate(features = "Core Shape\nFeatures")
      irreg_result <- pairwise_energy(
        clust1, clust2, df1, df2, irregularity_features) |>
        mutate(features = "Shape Irregularity\nFeatures")
      
      # Combine results and return
      return(bind_rows(all_result, core_result, irreg_result))
    
    },
    clust1=clust1, df=df, min_size=min_size, future.seed = TRUE
  )
  return(bind_rows(results))
}))


# Adjust pairwise energy results for multiple testing, then write to file.
# Treat different feature sets as independent tests therefore requiring 
# independent adjustment.
log_info("Adjusting for multiple testing ...")
all_energy_results <- adjust_pvalues(
  energy_results |> 
    filter(features == "All Shape\nFeatures" | is.na(features)))
core_energy_results <- adjust_pvalues(
  energy_results |> 
    filter(features == "Core Shape\nFeatures" | is.na(features)))
irreg_energy_results <- adjust_pvalues(
  energy_results |> 
    filter(features == "Shape Irregularity\nFeatures" | is.na(features)))
write_csv(all_energy_results, glue("{outdir}/energy_results_all_features.csv"))
write_csv(core_energy_results, glue("{outdir}/energy_results_core_features.csv"))
write_csv(irreg_energy_results, glue("{outdir}/energy_results_irreg_features.csv"))




### Inspect p-value distributions -----------------------------------------

log_info("Inspecting p-value distributions ...")

# Log number of significant results and write quick statistical summary
log_sig_result_summary <- function (df, name) {
  num_sig <- nrow(df |> filter(p_adj < 0.05))
  log_info(glue("{name}: {num_sig} significant results found"))
  p_adj <- df$p_adj[!is.na(df$p_adj)]
  write_csv(
    tibble(
      mean = mean(p_adj), sd = sd(p_adj), min = min(p_adj),
      q1 = quantile(p_adj, 0.25), median = median(p_adj),
      q3 = quantile(p_adj, 0.75), max = max(p_adj)
    ),
    glue("{outdir}/p_adj_summary_stats_{clean_string(name)}.csv")
  )
}
log_sig_result_summary(all_energy_results, "All Shape Features")
log_sig_result_summary(core_energy_results, "Core Shape Features")
log_sig_result_summary(irreg_energy_results, "Shape Irregularity Features")

# Draw horizontal violin plots illustrating adjusted & unadjusted p-value 
# distribution
feature_p_vals <- bind_rows(
  all_energy_results, core_energy_results, irreg_energy_results
  ) |> select(features, p, p_adj)
feature_p_vals$features <- factor(
  feature_p_vals$features,
  levels = c(
    "All Shape\nFeatures",
    "Core Shape\nFeatures",
    "Shape Irregularity\nFeatures"
  )
)
feature_p_vals <- feature_p_vals |> # Get df of features, type, p_value
  pivot_longer(cols = c(p, p_adj), names_to = "type", values_to = "p_value") |>
  mutate(type = ifelse(type == "p", "Raw", "Adjusted"))
plt <- ggplot(feature_p_vals, aes(x = features, y = p_value, fill = type)) +
  geom_half_violin( # raw
    data = feature_p_vals |> filter(type == "Raw"),
    side = "l", trim = TRUE, alpha = 0.7
  ) + 
  geom_half_violin( # adjusted
    data = feature_p_vals |> filter(type == "Adjusted"),
    side = "r", trim = TRUE, alpha = 0.7
  ) +
  geom_jitter(aes(color = type), width = 0.15, height = 0, size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Raw vs Adjusted P-Value Distributions",
    fill = "",
    color = ""
  )
svg(glue("{outdir}/p_distr.svg"), height = 12, width = 8)
print(plt)
invisible(dev.off())




### Prepare result dataframes for tile plots ------------------------------

log_info("Adding aesthetics ...")

add_aesthetic_data <- function (df) {
  # Drop all rows where clust2 <= clust1 - not required for triangle tile plot
  df <- df |> filter(clust1 > clust2)
  
  # Add tile aesthetics settings
  df <- df |>
    mutate(
      is_sig = (p_adj < 0.05),
      clust1 = factor(clust1),
      clust2 = factor(clust2)
    ) |>
    mutate(alpha = ifelse(is_sig, 1 - p_adj, 0)) |> # make non-sig invisible, and sig visible
    mutate(
      alpha = ifelse(is.na(p_adj), 1, alpha), # make grey/NA visible
      label = ifelse(is_sig, sprintf("%.3f", p_adj), "")
    )
  return(df)
  
}
all_energy_results <- add_aesthetic_data(all_energy_results)
core_energy_results <- add_aesthetic_data(core_energy_results)
irreg_energy_results <- add_aesthetic_data(irreg_energy_results)




### Create tile plots (core results!!!) -----------------------------------

log_info("Creating tile plots ...")

# Create triangular tile plot
make_tile_plot <- function (df, subtitle) {
  plt <- ggplot(df, aes(x = clust1, y = clust2)) +
    geom_tile(aes(fill = is_sig, alpha = alpha), color = "grey80") +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_manual(
      values = c(
        "TRUE"  = "#2ECC71", # significant -> green
        "FALSE" = "white" # non‑significant -> white
      ),
      guide = "none"
    ) +
    scale_alpha(range = c(0.2, 1), guide = "none") +
    geom_tile( # NA tiles need to be grey — override using a second layer
      data = subset(df, is.na(p_adj)), 
      fill = "grey80", color = "grey80"
    ) +
    scale_x_discrete(limits = rev) +
    coord_fixed() +
    labs(
      x = "Cluster",
      y = "Cluster",
      title = "Pairwise Energy Test Significance Matrix",
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 24, hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 18, margin = margin(t=30)),
      axis.title.y = element_text(size = 18, margin = margin(r=30)),
      panel.grid = element_blank()
    )
  png(
    glue("{outdir}/{clean_string(subtitle)}.png"), 
    height = 12, width = 12, units = "in", res = 300
  )
  print(plt)
  invisible(dev.off())
}

# Make tile plots for each feature set
make_tile_plot(all_energy_results, "All Shape Features")
make_tile_plot(core_energy_results, "Core Shape Features")
make_tile_plot(irreg_energy_results, "Shape Irregularity Features")


# Read PNGs
img_all  <- rasterGrob(
  readPNG("{outdir}/{clean_string('All Shape Features')}.png"), 
  interpolate = TRUE
)
img_core <- rasterGrob(
  readPNG("{outdir}/{clean_string('Core Shape Features')}.png"),
  interpolate = TRUE
)
img_irreg<- rasterGrob(
  readPNG("{outdir}/{clean_string('Shape Irregularity Features')}.png"),
  interpolate = TRUE
)

# Plot 1: all three plots in one row
grid_plt_1 <- grid.arrange(grobs = c(img_all, img_core, img_irreg), ncol = 3)
png(
  glue("{outdir}/result_grid_1.png"),
  width = 36, height = 12, units = "in", res = 300
)
grid.draw(grid_plt)
invisible(dev.off())

# Plot 2: arrange PNGs such that "All" plot is above the other two plots
grid_plt_2 <- grid.arrange(
  img_all, 
  arrangeGrob(img_core, img_irreg, ncol = 2),
  nrow = 2, heights = c(2, 1) # top row twice as tall as bottom row
)
png(
  glue("{outdir}/result_grid_2.png"),
  width = 20, height = 20, units = "in", res = 300
)
grid.draw(grid_plt_2)
invisible(dev.off())




log_info("Done!")


