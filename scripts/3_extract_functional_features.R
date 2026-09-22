# Extract core synaptic protein functional features.

# Group GO IDs into functional modules using a WGCNA-like method, write mapping
# between module label and member GO IDs, then create a feature dataframe 
# containing taxon, uniprot, and module scores representing proportion of module
# GO IDs associated with that protein.

outdir <- "results/3_extract_functional_features"
log_file <- "logs/3_extract_functional_features.log"
source("scripts/setup.R")
configure_logger(log_file)

suppressPackageStartupMessages(library(dendextend))
library(dynamicTreeCut)
library(ggplot2)
library(tidyr)

taxa <- c("fly","human","mouse","rat","worm","zebrafish")
plan(multisession, workers = parallel::detectCores(logical = FALSE))




### Define visualisation functions ----------------------------------------

draw_go_prevalence_props_hist <- function (df, title) {
  prevalence <- df |> select(-uniprot, -taxon) |>
    summarise(across(everything(), ~ mean(.x))) |> # get GO ID presence proportions
    pivot_longer(everything(), names_to = "go_id", values_to = "p") # go_id, p cols
  
  svg(glue("{outdir}/{clean_string(title)}.svg", width = 8, height = 8))
  hist(prevalence$p, xlab = "GO ID Prevalence (Proportion of Proteins)",
       main = title, breaks = seq(0, 1, by = 0.01), xlim = c(0,1))
  invisible(dev.off())
  
  return(prevalence)
}


visualise_matrix_value_distr <- function (mat, title, x_label) {
  
  # Collapse matrix into long tibble
  long_tibble <- mat |> as.data.frame() |>
    tibble::rownames_to_column("go_1") |>
    pivot_longer(cols = -go_1, names_to = "go_2", values_to = "score") |>
    filter(go_1 < go_2) # Drop diagonals and duplicates
  
  # Visualise score distribution
  svg(glue("{outdir}/{clean_string(title)}.svg"), width = 8, height = 8)
  hist(long_tibble$score, xlab = x_label, main = title,
       breaks = seq(0, 1, by = 0.01), xlim = c(0, 1))
  invisible(dev.off())
  
}




### Load data into matrix -------------------------------------------------

log_info("Loading data ...")

# Read in taxon GO annotation data (uniprot, go_id)
ann_df <- bind_rows(lapply(taxa, function(taxon) {
  func_file <- glue("results/1_filter_proteins/{taxon}.csv")
  taxon_df <- read_csv(func_file, show_col_types = FALSE) |>
    mutate(taxon = !!taxon)
  return(taxon_df)
}))

# Transform annnotation dataframe into sparse binary matrix with rows = proteins
# and columns = uniprot and a column per GO ID
go_presence_matrix <- ann_df |>
  mutate(value = 1L) |> # all existing GO IDs are assigned value = 1 (present)
  pivot_wider(
    id_cols = c(taxon, uniprot), # keep
    names_from = go_id, # each go_id value becomes a column
    values_from = value,
    values_fill = list(value = 0L) # value 0 = absent in protein
  )




### Filter out low variance GO IDs (those almost always present or absent)

log_info("Filtering out low-variance GO IDs ...")

prevalence <- draw_go_prevalence_props_hist(
  go_presence_matrix, "GO ID Prevalence Distribution Before Filtering"
)

# Remove go terms present in < 2% of proteins and > 99% of proteins
keep_go_ids <- prevalence |>
  filter(p >= 0.02, p <= 0.99) |> # based on histogram above
  pull(go_id) # transform to vector
go_presence_matrix_f <- go_presence_matrix |>
  select(uniprot, taxon, all_of(keep_go_ids))

prevalence <- draw_go_prevalence_props_hist(
  go_presence_matrix_f, "GO ID Prevalence Distribution After Filtering"
)




### Create a taxa-normalised GO ID Jaccard Similarity matrix --------------

log_info("Creating taxa-normalised GO ID Jaccard similarity matrix ...")

jaccard <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union <- sum(x == 1 | y == 1)
  return(intersection / union)
}

# Create raw Jaccard similarity matrices per taxa.
# Per taxa similarity matrices are needed to average similarities across taxa.
# An aggregate cross-taxa Jaccard similarity matrix constructed from the entire
# dataset would bias the similarity scores towards trends present in taxa with
# more proteins. I am aiming to detect GO modules that are informative
# across taxa without biasing towards any particular taxon. The per taxa matrices
# will be averaged to give each taxon even weight.
# The choice to use Jaccard similarity stems from the idea that GO IDs that
# tend to appear in the same proteins are correlated and likely to be part of a
# functional module, and one can use Jaccard similarity to evaluate how freque-
# ntly GO IDs tend to appear together. Jaccard similarity also treats co-absence
# as uninformative, which is handy for sparse matrices.
similarity_matrices <- future_lapply(taxa, function (taxon) {
  
  taxon_matrix <- go_presence_matrix_f |> 
    select(-taxon, -uniprot) |> filter(taxon == !!taxon)
  
  # Get tibble containing Jaccard scores between each GO ID
  j_scores <- lapply(colnames(taxon_matrix), function (go) {
    # Compute Jaccard similarity between this GO ID's prevalence and other GO IDs
    go_vec <- taxon_matrix[[go]] # Vector of 0 (absent) and 1 (present)
    scores <- lapply(colnames(taxon_matrix), function (go_target) {
      j_score <- jaccard(go_vec, taxon_matrix[[go_target]])
      list(go_1 = go, go_2 = go_target, jaccard = j_score)
    })
    return(tibble(
      go_1 = sapply(scores, `[[`, "go_1"),
      go_2 = sapply(scores, `[[`, "go_2"),
      jaccard = sapply(scores, `[[`, "jaccard")
    ))
  })
  j_scores <- bind_rows(j_scores)
  
  # Transform tibble (go_1, go_2, jaccard) into matrix
  similarity_matrix <- j_scores |> 
    pivot_wider(names_from = go_2, values_from = jaccard) |> 
    tibble::column_to_rownames("go_1") |> as.matrix()
  
  return(similarity_matrix)

})

# Average the taxa-specific similarity matrices into a unified similarity matrix
similarity_matrix <- Reduce("+", similarity_matrices) / length(similarity_matrices)

# Visualise similarity score distribution and heatmap
visualise_matrix_value_distr(
  similarity_matrix, x_label = "Jaccard Similarity Score",
  title = "Distribution of Pairwise Jaccard Similarity Scores (GO ID Presence)"
)
svg(glue("{outdir}/heatmap_jaccard.svg"), width = 8, height = 8)
heatmap(similarity_matrix, scale = "none")
invisible(dev.off())




### Apply wTOM transformation to similarity matrix ------------------------

# This transformation is adapted from the set-theoretic 1-M formula presented
# in Yip & Horvath (2007) and applied to a weighted Jaccard similarity matrix.

log_info("Applying wTOM transformation to similarity matrix ...")


# Returns vector of similarities between GO ID and all other GO IDs
get_go_pair_similarities <- function (go_id) {
  return(
    vapply(all_go_ids, FUN.VALUE = double(1), FUN = function(other_go_id) {
      similarity_matrix[go_id, other_go_id]
    }
  ))
}

# Calculate the weighted TOM score for the ij entry in a TOM matrix
wTOM <- function (go_1, go_2) {
  # Get vector of similarities between go_1 and all other GO IDs, same for go_2
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


all_go_ids <- rownames(similarity_matrix)
similarity_vectors <- future_lapply(all_go_ids, get_go_pair_similarities)
names(similarity_vectors) <- all_go_ids

# Apply wTOM transformation to jaccard matrix
wTOM_matrix <- future_sapply(all_go_ids, function(go_i) {
  return(sapply(all_go_ids, function(go_j) { wTOM(go_i, go_j) }))
})
wTOM_matrix <- t(wTOM_matrix) # want row cols to be go_i values (for my sanity)
rownames(wTOM_matrix) <- all_go_ids
colnames(wTOM_matrix) <- all_go_ids


# Visualise wTOM score distribution and heatmap
visualise_matrix_value_distr(
  wTOM_matrix, x_label = "wTOM Score",
  title = "Distribution of wTOM Scores"
)
svg(glue("{outdir}/heatmap_wtom.svg"), width = 8, height = 8)
heatmap(wTOM_matrix, scale = "none")
invisible(dev.off())




### Find GO ID modules ----------------------------------------------------

log_info("Finding GO modules ...")

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
svg(glue("{outdir}/go_modules_dendrogram.svg"), width = 16, height = 12)
par(mar = c(16, 4, 4, 4))
plot(dend, main = "GO ID Dendrogram")
colored_bars(module_colours, dend, rowLabels = "Module", y_scale = 0.2)
invisible(dev.off())


# Get and export mapping between module ID and member GO IDs
go_module_members <- tibble(module = modules_ints, go_id = names(modules_ints))
write_csv(go_module_members, glue("{outdir}/go_module_members.csv"))
go_module_members <- go_module_members |>
  group_by(module) |> summarise(go_ids = list(go_id), .groups = "drop")




### Compute module scores per protein -------------------------------------

# The module score represents the proportion of GO IDs associated with a given
# module that are also associated with a given protein.

# For each module, take the GO ID columns from go_presence_matrix_f, then
# compute the proportion of GO IDs present
module_scores <- sapply(go_module_members$go_ids, function(module_go_ids) {
  module_go_cols <- go_presence_matrix_f[, module_go_ids, drop = FALSE]
  return(rowMeans(module_go_cols)) # mean of binary values -> proportion present
})
colnames(module_scores) <- paste0("module_", go_module_members$module)

# Merge module_scores into go_presence_matrix_f to get a complete
# functional feature set with raw GO presence/absence and module score data
functional_features <- inner_join( # join on uniprot
  ann_df |> distinct(uniprot, taxon), 
  bind_cols(go_presence_matrix_f, tibble(as.data.frame(module_scores)))
)
write_csv(functional_features, glue("{outdir}/functional_features.csv"))




### Report degree of functional feature space reduction -------------------

num_go_ids_bef_filt <- ncol(go_presence_matrix |> select(starts_with("GO:"))) 
num_go_ids_aft_filt <- ncol(functional_features |> select(starts_with("GO:")))
num_modules <- ncol(functional_features |> select(starts_with("module_")))
log_info(
  glue("Functional feature space: \
  {num_go_ids_bef_filt} GO IDs before variance filtering; \
  {num_go_ids_aft_filt} GO IDs after variance filtering; \
  {num_modules} functional modules detected")
)




### Create GO Network graph coloured by module ----------------------------

# The objective is to understand how GO modules sit relative to each other - 
# e.g., are any modules nested within other modules, or scattered?

# TODO



### Create GO network graphs per GO module --------------------------------

# What GO terms make up each module, and how do they sit relative to each other?

# TODO








