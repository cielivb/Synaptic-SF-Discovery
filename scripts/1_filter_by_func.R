# Filter candidate proteins by function.

# For each taxon, curate a set of proteins present at neuron-neuron chemical 
# synapses as inferred via GO IDs and their ancestors.

# This script requires GO Slim Synapse, GO Basic Ontology, and GO 
# Taxon-Specific Datasets.


source("scripts/setup.R")
configure_logger("logs/1_filter_by_func.log")

library(jsonlite)
library(tidyr)



# Load and clean GO Slim Synapse
log_info("Cleaning GO Slim Synapse data ...")
goslim_synapse <- read_tsv(
  "data/goslim_synapse.tsv", col_names = TRUE, show_col_types = FALSE
) |>
  select(go_id = `?x`) |>
  mutate(go_id = sub("_", ":", basename(gsub("[<>]", "", go_id)))) |>
  arrange(go_id)
write_csv(goslim_synapse, "results/1_filter_by_func/goslim_synapse_go_ids.csv")


# Map GO Slim Synapse GO IDs to ancestor GO IDs
log_info("Mapping GO Slim Synapse GO IDs to their ancestor GO IDs ...")
system("python scripts/map_synaptic_go_ancestors.py")
ancestors <- fromJSON("results/1_filter_by_func/synaptic_go_ancestors.json")
ancestor_df <- enframe( # unpack named list ancestors into tibble
  ancestors, name = "go_id", value = "ancestor_ids"
)


# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))




# Filter proteins by function for each taxon in parallel
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
go_taxon_meta <- download_meta |> filter(category == "go_taxon")
log_info("Filtering each taxon's proteins ...")
future_lapply(go_taxon_meta$dest, function (go_taxon_file) {
  
  taxon <- sub("\\.gaf.gz$", "", basename(go_taxon_file)) # e.g., worm
  dest <- glue("results/1_filter_by_func/{taxon}.csv")
  dest_uniprot_only <- glue("results/1_filter_by_func/{taxon}_uniprot_only.csv")
  
  # Read in taxon annotations
  all_ann <- read_tsv(
    go_taxon_file, comment = "!", col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession
      X5 = col_character() # GO ID
    ),
    show_col_types = FALSE)
  names(all_ann) <- c("uniprot", "go_id")
  
  # Keep only those proteins with at least one GO ID in GO Slim Synapse
  synapse_ann <- all_ann |>
    group_by(uniprot) |>
    filter(any(!is.na(match(go_id, goslim_synapse$go_id)))) |>
    ungroup()
  
  
  # Merge ancestor GO IDs into synapse_ann df -> (uniprot, go_id, ancestor_ids)
  synapse_ann <- synapse_ann |> left_join(ancestor_df, by = "go_id")
  
  # Filter synapse_ann again to include only annotations involving those 
  # proteins with both the chemical synaptic transmission GO ID and the 
  # neuron-to-neuron synapse GO ID
  chemsyn_go <- "GO:0007268"
  nn_go <- "GO:0098984"
  chemsyn_proteins <- synapse_ann |> group_by(uniprot) |>
    filter(any(go_id == chemsyn_go | chemsyn_go %in% ancestor_ids)) |>
    distinct(uniprot)
  nn_proteins <- synapse_ann |> group_by(uniprot) |>
    filter(any(go_id == nn_go | nn_go %in% ancestor_ids)) |>
    distinct(uniprot)
  intersection <- inner_join(chemsyn_proteins, nn_proteins, by = "uniprot")
  filtered <- semi_join(synapse_ann, intersection, by = "uniprot")
  
  
  # Export filtered synaptic protein annotation data
  write_csv(filtered, dest, col_names = TRUE)
  write_csv(filtered |> distinct(uniprot), dest_uniprot_only, col_names = TRUE)
  log_info(glue("Filtered {taxon} proteins"))

})

log_info("Protein filtration script done")