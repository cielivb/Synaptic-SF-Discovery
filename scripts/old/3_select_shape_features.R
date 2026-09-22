# Inspect correlations and variances of candidate features for UMAP + HDBSCAN
# clustering, and prepare the highest signal features for clustering.


source("scripts/setup.R")
configure_logger("logs/3_select_shape_features.log")

library(GGally)
library(ggplot2)
library(tidyr)



dir.create(
  "results/3_select_shape_features", showWarnings = FALSE, recursive = TRUE
)



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
    "results/3_select_shape_features/{clean_string(title)}.svg"
  )
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  svg(outfile, width = 12, height = 18)
  print(plt)
  invisible(dev.off())
  
}


# Read in geometry descriptions for each taxon. Keep an unmodified geom tibble
# for merging cluster results against uniprot and raw geometric features later.
# Create a copy geom_df without uniprot column to update during feature selection.
taxa <- c("fly","human","mouse","rat","worm","zebrafish")
og_geom_df <- bind_rows(lapply(taxa, function (taxon) {
  geom_file <- glue("results/2_extract_shape_data/{taxon}_geom.csv")
  taxon_df <- read_csv(geom_file, show_col_types = FALSE) |>
    mutate(taxon = taxon)
  return(taxon_df)
}))
geom_df <- og_geom_df |> select(-uniprot)


# Visualise feature distributions before corrections
geom_df_f <- geom_df |> select(taxon, ends_with("_f"))
geom_df_hc <- geom_df |> select(taxon, ends_with("_hc"))
geom_df_misc <- geom_df |> select(-ends_with("_f"), -ends_with("_hc"))
plot_features_species_violin(geom_df_f, "Feature Set 1 Distributions by Taxa (Raw)")
plot_features_species_violin(geom_df_hc, "Feature Set 2 Distributions by Taxa (Raw)")
plot_features_species_violin(geom_df_misc, "Feature Set 3 Distributions by Taxa (Raw)")


# Correct for skewness in candidate features across taxa, then visualise.
# kurtosis_f, eccentricity_f, length_f, rg_f, kurtosis_hc, eccentricity_hc,
# length_hc, rg_hc
log_info("Correcting for skew ...")
no_skew_f <- geom_df_f |> 
  mutate(log_kurtosis_f = log(kurtosis_f),
         log_eccentricity_f = log(eccentricity_f),
         log_length_f = log(length_f),
         log_rg_f = log(rg_f)) |>
  select(-kurtosis_f, -eccentricity_f, -length_f, -rg_f)
no_skew_hc <- geom_df_hc |> 
  mutate(log_kurtosis_hc = log(kurtosis_hc),
         log_eccentricity_hc = log(eccentricity_hc),
         log_length_hc = log(length_hc),
         log_rg_hc = log(rg_hc)) |>
  select(-kurtosis_hc, -eccentricity_hc, -length_hc, -rg_hc)
no_skew_misc <- geom_df_misc
plot_features_species_violin(
  no_skew_f, "Feature Set 1 Distributions by Taxa (Skew-Corrected)")
plot_features_species_violin(
  no_skew_hc, "Feature Set 2 Distributions by Taxa (Skew-Corrected)")
plot_features_species_violin(
  no_skew_misc, "Feature Set 3 Distributions by Taxa (Skew-Corrected)")


# Use median-MAD-standardisation on geometric features to normalise
# variance while staying robust against outliers, so that PCA-based variance
# filtering can be used to eliminate low variance features
log_info("Normalising variances ...")
mad_scale <- function(v) {(v - median(v)) / mad(v)}
mad_scaled_f <- no_skew_f |>
  mutate(
    scaled_coeff_var_f = mad_scale(coeff_var_f),
    scaled_log_eccentricity_f = mad_scale(log_eccentricity_f),
    scaled_log_kurtosis_f = mad_scale(log_kurtosis_f),
    scaled_log_length_f = mad_scale(log_length_f),
    scaled_log_rg_f = mad_scale(log_rg_f)
  ) |> 
  select(
    -coeff_var_f, -log_eccentricity_f, -log_kurtosis_f, 
    -log_length_f, -log_rg_f
  )
mad_scaled_hc <- no_skew_hc |>
  mutate(
    scaled_coeff_var_hc = mad_scale(coeff_var_hc),
    scaled_log_eccentricity_hc = mad_scale(log_eccentricity_hc),
    scaled_log_kurtosis_hc = mad_scale(log_kurtosis_hc),
    scaled_log_length_hc = mad_scale(log_length_hc),
    scaled_log_rg_hc = mad_scale(log_rg_hc)
  ) |> 
  select(
    -coeff_var_hc, -log_eccentricity_hc, -log_kurtosis_hc,
    -log_length_hc, -log_rg_hc
  )
mad_scaled_misc <- no_skew_misc |>
  mutate(
    scaled_low_pLDDT = mad_scale(low_pLDDT),
    scaled_high_pLDDT = mad_scale(high_pLDDT)
  ) |>
  select(-low_pLDDT, -high_pLDDT)
plot_features_species_violin(
  mad_scaled_f, "Feature Set 1 Distributions by Taxa (Skew-Corrected & MAD-Scaled)")
plot_features_species_violin(
  mad_scaled_hc, "Feature Set 2 Distributions by Taxa (Skew-Corrected & MAD-Scaled")
plot_features_species_violin(
  mad_scaled_misc, "Feature Set 3 Distributions by Taxa (Skew-Corrected & MAD-Scaled")


# Get correlations and filter features to reduce redundancy
log_info("Decorrelating features ...")
mad_scaled_hc <- mad_scaled_hc |> select(-taxon) # taxon column already in scaled_f
mad_scaled_misc <- mad_scaled_misc |> select(-taxon)
mad_geom_df <- bind_cols(mad_scaled_f, mad_scaled_hc, mad_scaled_misc)
mad_geom_df_numeric <- mad_geom_df |> select(where(is.numeric))
plt <- ggcorr(mad_geom_df_numeric, nbreaks = 6, label = TRUE) + 
  labs(title = "MAD-Scaled & Skew-Corrected Feature Correlations Before Decorrelation") +
  theme(plot.title = element_text(size = 18, hjust = 0.2))
outfile <- "results/3_select_shape_features/corr_before.svg"
svg(outfile, width = 16, height = 16)
print(plt)
invisible(dev.off())
mad_geom_df <- mad_geom_df |> # Remove features based on correlation plot
  select(
    -scaled_low_pLDDT, # high correlation with high_pLDDT
    -scaled_log_rg_hc, # high corr with length_f
    -scaled_log_rg_f, # high corr with length_f
    -scaled_log_length_hc # high corr with lenght_f
  )
mad_geom_df_numeric <- mad_geom_df |> select(where(is.numeric))
plt <- ggcorr(mad_geom_df_numeric, nbreaks = 6, label = TRUE) + 
  labs(title = "MAD-Scaled & Skew-Corrected Feature Correlations After Decorrelation") +
  theme(plot.title = element_text(size = 18, hjust = 0.2))
outfile <- "results/3_select_shape_features/corr_after.svg"
svg(outfile, width = 16, height = 16)
print(plt)
invisible(dev.off())


# Run PCA on full geometric description dataframe to find variances per 
# feature, then filter out low variance features
log_info("Filtering out low variance features ...")
# Manual centering and scaling already done above
pca <- prcomp(mad_geom_df_numeric, center = FALSE, scale. = FALSE)
sink("results/3_select_shape_features/pca_summary.txt")
print(summary(pca)) # First 5 PCs explain ~90% of variance
sink()
# Get feature contributions to first 5 PCs (explain ~90% of variance)
feature_contrib <- rowSums((pca$rotation[, 1:5])^2)
contrib_df <- data.frame(
  feature = names(feature_contrib), contrib = feature_contrib)
plt <- contrib_df |> arrange(contrib) |> # sort by contrib value
  mutate(feature = factor(feature, levels=feature)) |> # update factor levels
  ggplot(aes(x = feature, levels = feature)) +
  geom_segment(aes(x=feature, xend = feature, yend=0, y=contrib)) +
  geom_point(aes(y=contrib), size = 4, color = "orange") +
  coord_flip() + theme_bw() + 
  labs(title = "Feature Contributions to First 5 Principal Components") +
  theme(plot.title = element_text(size = 18, hjust = 0.2))
outfile <- "results/3_select_shape_features/pca_feature_contrib.svg"
svg(outfile, width = 12, height = 8)
print(plt)
invisible(dev.off())
sink("results/3_select_shape_features/pca_feature_contrib.txt")
sort(feature_contrib, decreasing = TRUE)
sink()
# Keep only important features, determined by evaluating feature contribs and
# remaining correlations between features
mad_geom_df <- mad_geom_df |>
  select(
    taxon, scaled_high_pLDDT, scaled_log_eccentricity_f, 
    scaled_log_length_f, scaled_log_kurtosis_hc, scaled_coeff_var_hc
  )


# Restrict feature ranges to between 0 and 1 via min-max normalisation so they 
# are treated with equal importance in UMAP neighbourhood graph construction.
log_info("Min-max normalising filtered features ...")
min_max_norm <- function(v) {(v - min(v)) / (max(v) - min(v))}
mad_geom_df <- mad_geom_df |>
  mutate(across(where(is.numeric), min_max_norm)) |>
  select(-taxon) # Already present in og_geom_df (merging below)


# Merge scaled features with raw og_geom_df, select useful feature columns and
# metadata columns, then write to csv
log_info("Writing dataframe containing selected feature columns ...")
geom_df <- bind_cols(mad_geom_df, og_geom_df) |>
  select(
    taxon, uniprot, high_pLDDT, eccentricity_f, length_f,
    kurtosis_hc, coeff_var_hc, scaled_high_pLDDT, scaled_log_eccentricity_f,
    scaled_log_length_f, scaled_log_kurtosis_hc, scaled_coeff_var_hc
  )
write_csv(geom_df, "results/3_select_shape_features/selected.csv")



log_info("Shape feature selection script complete")