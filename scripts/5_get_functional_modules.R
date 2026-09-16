# Group GO IDs into functional modules using a WGCNA-like method, write mapping
# between module label and member GO IDs, then
# create a feature dataframe containing taxon, uniprot, geometric features,
# and module scores (representing proportion of module GO IDs present in that
# protein).


source("scripts/setup.R")
log_file <- "logs/5_get_functional_modules.log"
configure_logger(log_file)

library(tidyr)


dir.create("results/5_get_functional_modules", showWarnings = FALSE)


# Read in shape data (taxon, uniprot, geom-etc)
taxa <- c("fly", "human", "mouse", "rat", "worm", "zebrafish")
shape_df <- read_csv(
  "results/3_select_shape_features/selected.csv", show_col_types = FALSE
)

# Read in function data (uniprot, go_id)
function_df <- bind_rows(lapply(taxa, function(taxon) {
  func_file <- glue("results/1_filter_by_func/{taxon}.csv")
  taxon_df <- read_csv(func_file, show_col_types = FALSE)
  return(taxon_df)
}))

# Transform function data into sparse binary matrix with rows = proteins and
# columns = taxon, uniprot, and a column per GO ID
go_presence_matrix <- function_df |>
  mutate(value = 1L) |> # all existing GO IDs are associated with 1 (present)
  pivot_wider(
    id_cols = uniprot, # keep
    names_from = go_id, # each go_id value becomes a column
    values_from = value,
    values_fill = list(value = 0L) # absent in protein
  )


# Filter out low variance GO IDs (i.e., GO IDs that are almost always present or
# almost always absent)
prevalence <- go_presence_matrix |> select(-uniprot) |>
  summarise(across(everything(), ~ mean(.x))) |> # get presence proportion
  pivot_longer(everything(), names_to = "go_id", values_to = "p") # go_id, p cols
outfile <- "results/5_get_functional_modules/go_prevalence_before.svg"
svg(outfile, width = 8, height = 8)
hist(prevalence$p, xlab = "GO ID Prevalence (Proportion of Proteins)",
     main = "GO ID Prevalence Before Filtering Out Low-Variance GO IDs",
     breaks <- seq(0, 1, by = 0.01))
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
     breaks <- seq(0, 1, by = 0.01), xlim = c(0,1))
invisible(dev.off())




# Create a co-presence matrix, where each row and column is a GO ID, and the
# value is the Jaccard similarity score between the respective GO ID columns
# in the filtered GO presence matrix containing binary values.
# Use Jaccard because it treats co-absence as uninformative, and the GO ID
# presence vectors are likely to be sparse based on the above histograms.

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
  go_raw <- filtered_go_presence_matrix$go # vector of 0 and 1
  scores <- lapply(colnames(filtered_go_presence_matrix), function(go_target) {
    if (go_target == "uniprot") { 
      return(list(go_1 = NA, go_2 = NA, jaccard = NA)) 
    } else {
      jscore <- jaccard(go_raw, filtered_go_presence_matrix$go_target)
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

# Transform tibble (go_1, go_2, jaccard) into matrix
jaccard_matrix <- j_scores |> 
  pivot_wider(names_from = go_2, values_from = jaccard) |> 
  tibble::column_to_rownames("go_1") |> as.matrix()









