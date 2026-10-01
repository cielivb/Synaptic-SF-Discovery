# The core analysis is applied both across and within taxa, and consists of two 
# distinct stages:
# 1. Pairwise PERMANOVA between GMS clusters on protein shapes (centroid diffs)
# 2. Pairwise betadisper between GMS clusters on protein shapes (dispersion diffs)

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

# The above tile plot should also be generated for within taxa analyses, and
# tiles that match GMS clusters that are not present in the taxa should be
# greyed out.


outdir <- "results/8_core_analysis"
log_file <- "logs/8_core_analysis.log"
source("scripts/setup.R")
configure_logger(log_file)
plan(multisession, workers = parallel::detectCores(logical = FALSE))

library(vegan) # adonis2 (permanova) + betadisper




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
  select(mean_pLDDT, log_length, log_kurtosis, coeff_var, eccentricity) 

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
pw_permanova <- function(id_1, id_2, df_1, df_2) {

  # Merge dataframes then compute distance matrix
  df_pair <- bind_rows(df_1, df_2)
  dist_mat <- dist(df_pair |> select(-gms_cluster), method = "euclidean")
  
  result <- adonis2( # run PERMANOVA
    dist_mat ~ gms_cluster, data = df_pair, permutations = 9 # use 999 for real
  )
  
  return(tibble(
    id_1 = id_1, id_2 = id_2,
    F = result$F[1], R2 = result$R2[1], p = result$`Pr(>F)`[1]
  ))

}


# Run pairwise betadisper between two dataframes of shape features.
# Assumes valid inputs.
pw_betadisper <- function(df_1, df_2, id_1, id_2) {
  
  # Merge dataframes then compute distance matrix
  df_pair <- bind_rows(df_1, df_2)
  dist_mat <- dist(df_pair |> select(-gms_cluster), method = "euclidean")
  
  bd <- betadisper(dist_mat, df_pair$cluster) # Run betadisper
  anova_result <- anova(bd) # Run ANOVA on the betadisper object
  
  return(tibble(
    id_1 = id_1, id_2 = id_2,
    F = anova_result$`F value`[1], p = anova_result$`Pr(>F)`[1]
  ))
}





### Run analysis ----------------------------------------------------------

clusters <- unique(sort(df$gms_cluster))
min_size <- 10 # Require at least 10 proteins in a cluster to do a comparison


# Get unadjusted pairwise PERMANOVA results
permanova_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(clusters, function (clust2, clust1, df, min_size) {
    
    if (clust2 <= clust1) { 
      return(get_null_pw_permanova_result(clust1, clust2)) 
    }
    
    clust1_df <- df |> filter(gms_cluster == clust1) |> 
      select(gms_cluster, mean_pLDDT, log_length, eccentricity, coeff_var, log_kurtosis)
    clust2_df <- df |> filter(gms_cluster == clust2) |>
      select(gms_cluster, mean_pLDDT, log_length, eccentricity, coeff_var, log_kurtosis)
    
    if (nrow(clust1_df) < min_size || nrow(clust2_df) < min_size) {
      return(get_null_pw_permanova_result(clust1, clust2)) 
    }
    
    return(pw_permanova(clust_1, clust_2, clust1_df, clust2_df))
  
  }, clust1=clust1, df=df, min_size=min_size)
  return(bind_rows(results))
}))


# Get unadjusted pairwise betadisper results
betadisper_results <- bind_rows(lapply(clusters, function (clust1) {
  results <- future_lapply(clusters, function (clust2, clust1, df, min_size) {
    
    if (clust2 <= clust1) { 
      return(get_null_pw_betadisper_result(clust1, clust2)) 
    }
    
    clust1_df <- df |> filter(gms_cluster == clust1) |> 
      select(gms_cluster, mean_pLDDT, log_length, eccentricity, coeff_var, log_kurtosis)
    clust2_df <- df |> filter(gms_cluster == clust2) |>
      select(gms_cluster, mean_pLDDT, log_length, eccentricity, coeff_var, log_kurtosis)
    
    if (nrow(clust1_df) < min_size || nrow(clust2_df) < min_size) {
      return(get_null_pw_permanova_result(clust1, clust2)) 
    }
    
    return(pw_betadisper(clust_1, clust_2, clust1_df, clust2_df))
    
  }, clust1=clust1, df=df, min_size)
  return(bind_rows(results))
}))

