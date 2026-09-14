# Filter candidate proteins by function.

# For each taxon, curate a set of proteins present at neuron-neuron chemical 
# synapses as inferred via GO IDs and their ancestors.

# This script requires GO Slim Synapse, GO Basic Ontology, and GO 
# Taxon-Specific Datasets.


source("scripts/setup.R")
configure_logger("logs/1_filter_by_func.log")

library(jsonlite)


# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))


# Load and clean GO Slim Synapse
#go_basic <- fromJSON("data/go-basic.json")
goslim_synapse <- read_tsv(
  "data/goslim_synapse.tsv", col_names = TRUE, show_col_types = FALSE
) |>
  select(go_id = `?x`) |>
  mutate(go_id = sub("_", ":", basename(gsub("[<>]", "", go_id)))) |>
  arrange(go_id)
write_csv(goslim_synpase, "results/1_filter_by_func/goslim_synapse_go_ids.csv")


# Create JSON map between GO IDs in GO Slim Synapse and ancestor GO IDs
system("python scripts/map_synaptic_go_ancestors.py")


# Filter proteins by function for each taxon in parallel
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
go_taxon_meta <- download_meta |> select(category = "go_taxon")
future_lapply(go_taxon_meta$dest, function (go_taxon_file) {
  
  taxon <- sub(".gaf.gz", "", basename(go_taxon_file)) # e.g., worm
  dest <- glue("results/1_filter_by_func/{taxon}.json")
  
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
  
  # Add ancestor GO IDs to synapse_ann dataframe
  
  # Filter synapse_ann again to include only those proteins with both the 
  # chemical synapse GO ID and the neuron-to-neuron synapse GO ID
  
})






