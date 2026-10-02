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

library(energy)

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
  result <- eqdist.etest(rbind(X, Y), sizes = c(nrow(X), nrow(Y)), R = 9)
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

clusters <- unique(sort(df$gms_cluster))
min_size <- 10 # Require at least 10 proteins in a cluster to do a comparison
feature_cols <- c("mean_pLDDT", "log_length", "log_kurtosis", "coeff_var", "eccentricity")

# Get unadjusted pairwise energy test results
log_info("Running pairwise energy tests ...")
energy_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(
    clusters,
    function (clust2, clust1, df, min_size, feature_cols) {
      
      if (clust2 <= clust1)
        return(get_na_energy_result(clust1, clust2))
      
      df1 <- df |> filter(gms_cluster == clust1)
      df2 <- df |> filter(gms_cluster == clust2)
      if (nrow(df1) < min_size || nrow(df2) < min_size) {
        return(get_na_energy_result(clust1, clust2))
      }
      
      return(pairwise_energy(clust1, clust2, df1, df2, feature_cols))
    
    },
    clust1=clust1, df=df, min_size=min_size, feature_cols=feature_cols,
    future.seed = TRUE
  )
  return(bind_rows(results))
}))


# Adjust pairwise energy results for multiple testing, then write to file
log_info("Adjusting for multiple testing ...")
energy_results <- adjust_pvalues(energy_results)
write_csv(energy_results, glue("{outdir}/energy_results.csv"))




### Create tile plot (core result!!!) -------------------------------------





log_info("Done!")




