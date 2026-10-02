# Group GO IDs into functional modules using a WGCNA-like method, write mapping
# between module label and member GO IDs, then
# create a feature dataframe containing taxon, uniprot, geometric features,
# and module scores (representing proportion of module GO IDs present in that
# protein).


source("scripts/setup.R")
log_file <- "logs/5_get_functional_modules.log"
configure_logger(log_file)

suppressPackageStartupMessages(library(dendextend))
library(dynamicTreeCut)
library(ggplot2)
library(tidyr)

dir.create("results/5_get_functional_modules", showWarnings = FALSE)
plan(multisession, workers = parallel::detectCores(logical = FALSE))
taxa <- c("fly","human","mouse","rat","worm","zebrafish")


# Read in function data (uniprot, go_id)
function_df <- bind_rows(lapply(taxa, function(taxon) {
  func_file <- glue("results/1_filter_by_func/{taxon}.csv")
  taxon_df <- read_csv(func_file, show_col_types = FALSE) |>
    mutate(taxon = !!taxon)
  return(taxon_df)
}))

str(function_df)

# Transform function data into sparse binary matrix with rows = proteins and
# columns = taxon, uniprot, and a column per GO ID
go_presence_matrix <- function_df |>
  select(-taxon) |>
  mutate(value = 1L) |> # all existing GO IDs are associated with 1 (present)
  pivot_wider(
    id_cols = uniprot, # keep
    names_from = go_id, # each go_id value becomes a column
    values_from = value,
    values_fill = list(value = 0L) # absent in protein
  )


# Filter out low variance GO IDs (i.e., GO IDs that are almost always present or
# almost always absent)
log_info("Filtering out low variance GO IDs ...")
prevalence <- go_presence_matrix |> select(-uniprot) |>
  summarise(across(everything(), ~ mean(.x))) |> # get presence proportion
  pivot_longer(everything(), names_to = "go_id", values_to = "p") # go_id, p cols
outfile <- "results/5_get_functional_modules/go_prevalence_before.svg"
svg(outfile, width = 8, height = 8)
hist(prevalence$p, xlab = "GO ID Prevalence (Proportion of Proteins)",
     main = "GO ID Prevalence Before Filtering Out Low-Variance GO IDs",
     breaks = seq(0, 1, by = 0.01))
invisible(dev.off())
keep_go_ids <- prevalence |>
  # remove go terms present in < 2% of proteins and > 99% of proteins
  filter(p >= 0.02, p <= 0.99) |> # based on histogram above
  pull(go_id) # transform to vector

filtered_go_presence_matrix <- go_presence_matrix |>
  select(uniprot, all_of(keep_go_ids)) # discard low variance GO IDs
filtered_prevalence <- filtered_go_presence_matrix |> select(-uniprot) |>
  summarise(across(everything(), ~ mean(.x))) |> # get presence proportion
  pivot_longer(everything(), names_to = "go_id", values_to = "p") # go_id, p cols
outfile <- "results/5_get_functional_modules/go_prevalence_after.svg"
svg(outfile, width = 8, height = 8)
hist(filtered_prevalence$p, xlab = "GO ID Prevalence (Proportion of Proteins)",
     main = "GO ID Prevalence After Filtering Out Low-Variance GO IDs",
     breaks = seq(0, 1, by = 0.01), xlim = c(0,1))
invisible(dev.off())




# Create a co-presence matrix, where each row and column is a GO ID, and the
# value is the Jaccard similarity score between the respective GO ID columns
# in the filtered GO presence matrix containing binary values.
# Use Jaccard because it treats co-absence as uninformative, and the GO ID
# presence vectors are likely to be sparse based on the above histograms.
log_info("Creating GO ID co-presence matrix ...")

jaccard <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union <- sum(x == 1 | y == 1)
  return(intersection / union)
}

# Get tibble containing Jaccard scores between each GO ID
j_scores <- future_lapply(colnames(filtered_go_presence_matrix), function (go) {
  # Skip uniprot column
  if (go == "uniprot") { return(tibble(go_1 = NA, go_2 = NA, jaccard = NA)) }
  
  # Compute Jaccard similarity between this GO ID's raw prevalence and other GO IDs
  go_raw <- filtered_go_presence_matrix[[go]] # vector of 0 and 1
  scores <- lapply(colnames(filtered_go_presence_matrix), function(go_target) {
    if (go_target == "uniprot") { 
      return(list(go_1 = NA, go_2 = NA, jaccard = NA)) 
    } else {
      j_score <- jaccard(go_raw, filtered_go_presence_matrix[[go_target]])
      list(go_1 = go, go_2 = go_target, jaccard = j_score)
    }
  })
  
  return(tibble(
    go_1 = sapply(scores, `[[`, "go_1"),
    go_2 = sapply(scores, `[[`, "go_2"),
    jaccard = sapply(scores, `[[`, "jaccard")
  ))
  
})
j_scores <- bind_rows(j_scores) |> drop_na() |> distinct()

# Visualise similarity score distribution
outfile <- "results/5_get_functional_modules/jaccard_score_distribution.svg"
svg(outfile, width = 8, height = 8)
hist(j_scores$jaccard, xlab = "Jaccard Similarity Score",
     main = "Distribution of Pairwise Jaccard Similarity Scores",
     breaks = seq(0, 1, by = 0.01), xlim = c(0,1))
invisible(dev.off())

# Transform tibble (go_1, go_2, jaccard) into matrix
jaccard_matrix <- j_scores |> 
  pivot_wider(names_from = go_2, values_from = jaccard) |> 
  tibble::column_to_rownames("go_1") |> as.matrix()



# Apply weighted TOM-like transformation, adapted from the set-theoretic 1-M 
# formula presented in Yip & Horvath (2007) to a weighted jaccard matrix.
log_info("Applying weighted TOM transformation ...")
all_go_ids <- rownames(jaccard_matrix)

# Returns vector of similarities between GO ID and all other GO IDs
get_go_pair_similarities <- function (go_id) {
  return(
    vapply(all_go_ids, FUN.VALUE = double(1), FUN = function(other_go_id) {
      jaccard_matrix[go_id, other_go_id]
    }
  ))
}
similarity_vectors <- future_lapply(all_go_ids, get_go_pair_similarities)
names(similarity_vectors) <- all_go_ids

# Calculate the weighted TOM score for the ij entry in a TOM matrix
wTOM <- function (go_1, go_2) {
  # Get vector of similarities between go_1 and all other GO IDs, likewise for
  # go_2
  v1 <- similarity_vectors[[go_1]]
  v2 <- similarity_vectors[[go_2]]
  
  # Multiply elementwise and then sum result (i.e., get dot product of similarity
  # score vectors) to get numerator
  numerator <- sum(v1 * v2)
  
  # Choose mean over min or max for denominator to avoid artificially inflating
  # or suppressing similiarities and mitigate large "cliffs" in distribution
  denominator <- mean(sum(v1), sum(v2)) + 1e-6
  
  return(numerator / denominator)
}

# Apply wTOM transformation to jaccard matrix
wTOM_matrix <- future_sapply(all_go_ids, function(go_i) {
  return(sapply(all_go_ids, function(go_j) { wTOM(go_i, go_j) }))
})
wTOM_matrix <- t(wTOM_matrix) # want row cols to be go_i values (for my sanity)
rownames(wTOM_matrix) <- all_go_ids
colnames(wTOM_matrix) <- all_go_ids

# Collapse wTOM matrix into long tibble for wTOM score distribution visualisation
wTOM_tibble <- wTOM_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("go_1") |>
  pivot_longer(cols = -go_1, names_to = "go_2", values_to = "wTOM_score") |>
  filter(go_1 < go_2) # Drop diagonal and duplicates

# Visualise wTOM score distribution
outfile <- "results/5_get_functional_modules/wTOM_score_distribution.svg"
svg(outfile, width = 8, height = 8)
hist(wTOM_tibble$wTOM_score, xlab = "wTOM score",
     main = "Distribution of wTOM scores",
     breaks = seq(0, 1, by = 0.01), xlim = c(0,1),)
invisible(dev.off())

# Visualise wTOM matrix heatmaps
outfile <- "results/5_get_functional_modules/wTOM_score_heatmap.svg"
svg(outfile, width = 8, height = 8)
heatmap(wTOM_matrix, Rowv = NA, Colv = NA, scale = "none")
invisible(dev.off())
outfile <- "results/5_get_functional_modules/wTOM_score_heatmap_2.svg"
svg(outfile, width = 8, height = 8)
heatmap(wTOM_matrix, scale = "none")
invisible(dev.off())




# Cluster GO IDs into functional modules based off wTOM matrix

log_info("Finding GO ID clusters ...")
dissimilarity_matrix <- 1 - wTOM_matrix
tree <- hclust(as.dist(dissimilarity_matrix), method = "average")

# Cut into modules and plot dendrogram with associated modules
modules <- cutreeDynamic(
  dendro = tree,
  distM = dissimilarity_matrix,
  deepSplit = 3, # Want to err on the side of more fine-grained modules
  pamRespectsDendro = TRUE, # re-assign borderline GO IDs if helpful
  minClusterSize = 5 # Produces a relatively coherent clustering
)

# For visualisation, create dendrogram, assigning empty modules to a dedicated 
# NA module
modules_ints <- as.integer(factor(modules))
module_NA <- max(modules_ints, na.rm = TRUE) + 1
modules_ints[is.na(modules_ints)] <- module_NA
names(modules_ints) <- rownames(wTOM_matrix)
dend <- as.dendrogram(tree)
dend <- set_labels(dend, rep("", length(labels(dend))))

palette <- viridis::viridis(length(unique(modules_ints)))
module_colours <- palette[modules_ints]

# Plot
outfile <- "results/5_get_functional_modules/go_modules_dendrogram.svg"
svg(outfile, width = 16, height = 12)
par(mar = c(16, 4, 4, 4))
plot(dend, main = "GO ID Dendrogram")
colored_bars(module_colours, dend, rowLabels = "Module", y_scale = 0.2)
invisible(dev.off())

# Get and export mapping between module ID and member GO IDs
module_go_map <- tibble(module = modules_ints, go_id = names(modules_ints))
write_csv(module_go_map, "results/5_get_functional_modules/module_go_map.csv")
module_go_map <- module_go_map |>
  group_by(module) |> summarise(go_ids = list(go_id), .groups = "drop")

# For each module in each protein, get the percentage of GO IDs in that module
# that are linked to that protein ('module strength')

log_info("Computing module strengths per protein ...")

# For each module, take the GO ID columns from filtered_go_presence_matrix, then
# compute the proportion of GO IDs present
module_strengths <- sapply(module_go_map$go_ids, function(module_go_ids) {
  module_go_cols <- filtered_go_presence_matrix[, module_go_ids, drop = FALSE]
  return(rowMeans(module_go_cols)) # mean of binary values -> proportion present
})
colnames(module_strengths) <- paste0("module_", module_go_map$module)

# Merge module_strengths into filtered_go_presence_matrix to get a complete
# functional feature set with raw GO presence/absence and module strengths
functional_features <- inner_join(
  function_df |> distinct(uniprot, taxon), 
  bind_cols(
    filtered_go_presence_matrix, 
    tibble(as.data.frame(module_strengths))
  )
)
write_csv(
  functional_features, "results/5_get_functional_modules/functional_features.csv")


# Visualise module strength distributions over entire dataset
plot_module_strengths <- function(df, title) {
  # Code adapted from plot_features_species_violin in select_shape_features.R
  df_long <- df |> select(-starts_with("GO:"), -uniprot) |>
    pivot_longer(cols = -taxon, names_to = "feature", values_to="value")
  plt <- ggplot(df_long, aes(x=taxon, y=value, fill=taxon)) +
    geom_violin(trim=FALSE) +
    stat_summary(fun=median, geom="point", size=1.5, color="black") +
    facet_grid(rows=vars(feature), scales="free_y") +
    labs(title = title) +
    theme(strip.text.y = element_text(angle=0))
  outfile <- glue("results/5_get_functional_modules/{clean_string(title)}.svg")
  svg(outfile, width=12, height = 18)
  print(plt)
  invisible(dev.off())
}
plot_module_strengths(functional_features,  "Module Strength Distributions by Taxa")


# Report degree of functional feature space reduction
num_go_ids_bef_filt <- ncol(go_presence_matrix |> select(starts_with("GO:"))) 
num_go_ids_aft_filt <- ncol(functional_features |> select(starts_with("GO:")))
num_modules <- ncol(functional_features |> select(starts_with("module_")))
log_info(
  glue("Functional feature space: \
  {num_go_ids_bef_filt} GO IDs before variance filtering; \
  {num_go_ids_aft_filt} GO IDs after variance filtering; \
  {num_modules} functional modules detected")
)

head(functional_features |> select(starts_with("module_")))
log_info("Script complete")



