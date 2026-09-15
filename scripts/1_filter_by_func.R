# Filter candidate proteins by function.

# For each taxon, curate a set of proteins present at neuron-neuron chemical 
# synapses as inferred via GO IDs and their ancestors.

# This script requires GO Slim Synapse, GO Basic Ontology, and GO 
# Taxon-Specific Datasets.


source("scripts/setup.R")
log_file <- "logs/1_filter_by_func.log"
configure_logger(log_file)

library(jsonlite)
library(tidyr)


dir.create("results/1_filter_by_func", showWarnings = FALSE, recursive = TRUE)


# Load and clean GO Slim Synapse
log_info("Cleaning GO Slim Synapse data ...")
goslim_synapse <- read_tsv(
  "data/goslim_synapse.tsv", col_names = TRUE, show_col_types = FALSE
) |>
  select(go_id = `?x`) |>
  mutate(go_id = sub("_", ":", basename(gsub("[<>]", "", go_id)))) |>
  arrange(go_id)
write_csv(goslim_synapse, "results/1_filter_by_func/goslim_synapse_go_ids.csv")


# Map synaptic GO IDs to ancestor GO IDs
log_info("Mapping synaptic GO IDs to ancestor GO IDs ...")
system("python scripts/map_synaptic_go_ancestors.py")
ancestors <- fromJSON("results/1_filter_by_func/synaptic_go_ancestors.json")
ancestor_df <- tibble::enframe( # unpack named list ancestors into tibble
  ancestors, name = "go_id", value = "ancestor_ids"
)
ancestor_df$ancestor_ids <- lapply(ancestor_df$ancestor_ids, function(x) {
    if (is.null(x)) {character()} else if (is.list(x)) {unlist(x)} else {x}
})


# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))




# Filter proteins by function for each taxon in parallel
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
go_taxon_meta <- download_meta |> filter(category == "go_taxon")
log_info("Filtering each taxon's proteins ...")
done <- future_lapply(go_taxon_meta$destination, function (go_taxon_file) {
  configure_logger(log_file, worker = TRUE)
  
  taxon <- sub("\\.gaf.gz$", "", basename(go_taxon_file)) # e.g., worm
  dest <- glue("results/1_filter_by_func/{taxon}.csv")
  dest_uniprot_only <- glue("results/1_filter_by_func/{taxon}_uniprot_only.csv")
  
  # Read in taxon annotations
  ann <- read_tsv(
    go_taxon_file, comment = "!", col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession
      X5 = col_character() # GO ID
    ),
    show_col_types = FALSE)
  names(ann) <- c("uniprot", "go_id")
  
  # Expand annotations to include ancestors
  anc_ann <- ann |>
    left_join(ancestor_df, by = "go_id") |>
    select(-go_id) |>
    unnest(ancestor_ids) |>
    rename(go_id = ancestor_ids)
  all_ann <- bind_rows(anc_ann, ann) |> arrange(uniprot, go_id)
  
  # Include only annotations involving those proteins with both the chemical 
  # synaptic transmission GO ID and the neuron-to-neuron synapse GO ID
  chemsyn_go <- "GO:0007268"
  nn_go <- "GO:0098984"
  synaptic_proteins <- all_ann |> group_by(uniprot) |>
    filter(any(go_id == chemsyn_go) & any(go_id == nn_go)) |> distinct(uniprot)
  filtered_ann <- semi_join(all_ann, synaptic_proteins, by = "uniprot")
  num_proteins <- filtered_ann |> distinct(uniprot) |> nrow()
  log_info(glue("{taxon}: {nrow(filtered_ann)} annotations & {num_proteins} proteins after filter"))
  
  # Export filtered synaptic protein annotation data
  write_csv(filtered_ann, dest, col_names = TRUE)
  write_csv(filtered_ann |> distinct(uniprot), dest_uniprot_only, col_names = TRUE)
  log_info(glue("Filtered {taxon} proteins"))

})

log_info("Protein filtration script done")