# The core analysis is applied both across taxa only, and consists of two 
# distinct stages:
# 1. Pairwise PERMANOVA between GMS clusters on protein shapes (centroid diffs)
# 2. Pairwise betadisper between GMS clusters on protein shapes (dispersion diffs)

# Within taxa is not planned as many clusters will be of an inadequate sample
# size within taxa.

# Want to know if the geometric shapes of proteins in cluster
# A are statistically different from those in cluster B -> use raw/absolute 
# geometric features (except length and kurtosis, which will be log transformed
# due to high skew - still biologically meaningful)

# A tile plot illustrating significance of difference between pairwise cluster 
# comparisons is also planned. It will be an nxn tile plot where n is the number of 
# GMS clusters. The lower triangle should contain betadisper results which should be
# coloured with a blue hue. The upper triangle should contain PERMANOVA results
# which should be coloured with a red hue. Intensity of hue should match
# significance. Ideally, tiles matching the same cluster should be greyed out.




outdir <- "results/8_core_analysis"
log_file <- "logs/8_core_analysis.log"
source("scripts/setup.R")
configure_logger(log_file)
plan(multisession, workers = parallel::detectCores(logical = FALSE))

library(vegan) # adonis2 (permanova) + betadisper

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
get_null_pw_permanova_result <- function (id_1, id_2) {
  tibble(id_1, id_2, F = NA, R2 = NA, p = NA)
}
get_null_pw_betadisper_result <- function(id_1, id_2) {
  tibble(id_1, id_2, F = NA, p = NA)
}


# Run pairwise permanova between two dataframes of shape features.
# Assumes valid inputs (i.e., adequate size)
pw_permanova <- function(id_1, id_2, df, full_dist_matrix) {
  set.seed(254)
  
  # Get indices for each protein in dataframe assigned to either id_1 or id_2,
  # then extract corresponding rows and columns from full distance matrix,
  # and subset the dataframe at the row/protein indices
  row_indices <- which(df$gms_cluster %in% c(id_1, id_2))
  distance_matrix <- full_dist_matrix[row_indices, row_indices]
  df_pair <- df[idx, ]
  
  # Run PERMANOVA
  result <- adonis2(
    distance_matrix ~ gms_cluster, data = df_pair, permutations = 9 # 199 later
  )
  
  return(tibble(
    id_1 = id_1, id_2 = id_2,
    F = result$F[1],
    R2 = result$R2[1],
    p = result$`Pr(>F)`[1]
  ))

}


# Run pairwise betadisper between two dataframes of shape features.
# Assumes valid inputs.
pw_betadisper <- function(id_1, id_2, df_1, df_2, full_dist_matrix) {
  set.seed(254)
  
  # Get indices for each protein in dataframe assigned to either id_1 or id_2,
  # then extract corresponding rows and columns from full distance matrix,
  # and subset the dataframe at the row/protein indices
  row_indices <- which(df$gms_cluster %in% c(id_1, id_2))
  distance_matrix <- full_dist_matrix[row_indices, row_indices]
  df_pair <- df[idx, ]
  
  # Run betadisper
  bd <- betadisper(dist_sub, df_pair$gms_cluster)
  anova_result <- anova(bd)
  
  return(tibble(
    id_1 = id_1, id_2 = id_2,
    F = anova_result$`F value`[1],
    p = anova_result$`Pr(>F)`[1]
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

log_info("Computing full distance matrix ...")
feature_cols <- c("mean_pLDDT", "log_length", "eccentricity", "coeff_var", "log_kurtosis")
full_dist_matrix <- as.matrix(dist(df[, feature_cols], method = "euclidean"))


# Get unadjusted pairwise PERMANOVA results
log_info("Running pairwise PERMANOVA ...")
permanova_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(
    clusters, 
    function (clust2, clust1, df, min_size, full_dist_matrix) {
    
      # Return NA result if appropriate
      if (clust2 <= clust1)
        return(get_null_pw_permanova_result(clust1, clust2))
      if (sum(df$gms_cluster == clust1) < min_size ||
          sum(df$gms_cluster == clust2) < min_size)
        return(get_null_pw_permanova_result(clust1, clust2))
      
      return(pw_permanova(clust1, clust2, df, full_dist_mat))
  
    }, 
    clust1=clust1, df=df, min_size=min_size, 
    full_dist_matrix=full_distance_matrix, 
    future.seed = TRUE
  )
  return(bind_rows(results))
}))


# Get unadjusted pairwise betadisper results
log_info("Running pairwise betadisper ...")
betadisper_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(
    clusters, 
    function (clust2, clust1, df, min_size, full_dist_matrix) {
    
      # Return NA result if appropriate
      if (clust2 <= clust1)
        return(get_null_pw_betadisper_result(clust1, clust2))
      if (sum(df$gms_cluster == clust1) < min_size ||
          sum(df$gms_cluster == clust2) < min_size)
        return(get_null_pw_permanova_result(clust1, clust2))
      
      return(pw_betadisper(clust1, clust2, df, full_dist_matrix))
      
    },
    clust1=clust1, df=df, min_size, 
    full_distance_matrix=full_distance_matrix, 
    future.seed = TRUE
  )
  return(bind_rows(results))
}))


# Adjust pairwise PERMANOVA and betadisper results for multiple testing,
# then write to file
log_info("Adjusting for multiple testing ...")
permanova_results <- adjust_pvalues(permanova_results)
betadisper_results <- adjust_pvalues(betadisper_results)
write_csv(permanova_results, glue("{outdir}/permanova_results.csv"))
write_csv(betadisper_results, glue("{outdir}/betadisper_results.csv"))




### Create tile plot (core result!!!) -------------------------------------










