# Filter candidate proteins by function.

# For each taxon, curate a set of proteins present at neuron-neuron chemical 
# synapses as inferred via GO IDs and their ancestors.

# This script requires GO Slim Synapse, GO Basic Ontology, and GO 
# Taxon-Specific Datasets.


outdir <- "results/1_filter_proteins"
log_file <- "logs/1_filter_proteins.log"
source("scripts/setup.R")
configure_logger(log_file)

library(jsonlite)
library(tidyr)




# Launch cluster for parallel processing
# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))




# Load and clean GO Slim Synapse
log_info("Cleaning GO Slim Synapse data ...")
goslim_synapse <- read_tsv(
  "data/goslim_synapse.tsv", col_names = TRUE, show_col_types = FALSE
) |>
  select(go_id = `?x`) |> # go id column is literally called ?x
  mutate(go_id = sub( # extract go id from URI format
    "_", ":", basename(gsub("[<>]", "", go_id)))
  ) |>
  arrange(go_id) # Sort GO IDs
write_csv(goslim_synapse, glue("{outdir}/goslim_synapse_go_ids.csv"))


# Map synaptic GO IDs to ancestor GO IDs
log_info("Mapping synaptic GO IDs to ancestor GO IDs ...")
system("python scripts/map_synaptic_go_ancestors.py")
ancestors <- fromJSON(glue("{outdir}/synaptic_go_ancestors.json"))
ancestor_df <- tibble::enframe( # unpack named list ancestors into tibble
  ancestors, name = "go_id", value = "ancestor_ids"
)
# Some ancestor ID values may be null, single items, or stored in a list.
# Coerce these different data types into characters.
ancestor_df$ancestor_ids <- lapply(ancestor_df$ancestor_ids, function(x) {
    if (is.null(x)) {character()} else if (is.list(x)) {unlist(x)} else {x}
})




# Filter proteins for each taxon in parallel, keeping only proteins associated
# with chemical neuron-neuron synapses.

download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
go_taxon_meta <- download_meta |> filter(category == "go_taxon")

log_info("Finding proteins associated with neuron-neuron synapses ...")
done <- future_lapply(go_taxon_meta$destination, function (go_taxon_file, outdir) {
  configure_logger(log_file, worker = TRUE)
  
  # Extract taxon name from taxon-specific GO file and define output filenames
  taxon <- sub("\\.gaf.gz$", "", basename(go_taxon_file)) # e.g., worm
  dest <- glue("{outdir}/{taxon}.csv")
  dest_uniprot_only <- glue("{outdir}/{taxon}_uniprot_only.csv")
  
  # Read in taxon GO annotations present in taxon-specific GO file
  ann <- read_tsv(
    go_taxon_file, comment = "!", col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession
      X5 = col_character() # GO ID
    ),
    show_col_types = FALSE)
  names(ann) <- c("uniprot", "go_id")
  
  # Expand GO annotations to include ancestor GO IDs
  ancestor_ann <- ann |>
    left_join(ancestor_df, by = "go_id") |>
    select(-go_id) |>
    unnest(ancestor_ids) |>
    rename(go_id = ancestor_ids)
  all_ann <- bind_rows(ancestor_ann, ann) |> arrange(uniprot, go_id)
  
  # Include only annotations involving those proteins with both the chemical 
  # synaptic transmission GO ID and the neuron-to-neuron synapse GO ID
  num_proteins <- all_ann |> distinct(uniprot) |> nrow()
  log_info(glue(
    "{taxon}: {nrow(all_ann)} annotations & {num_proteins} proteins before filter"
  ))
  chemsyn_go <- "GO:0007268" # Chemical synaptic transmission
  nn_go <- "GO:0098984" # Neuron-to-neuron synapse
  target_proteins <- all_ann |> group_by(uniprot) |>
    filter(any(go_id == chemsyn_go) & any(go_id == nn_go)) |> distinct(uniprot)
  filtered_ann <- semi_join(all_ann, target_proteins, by = "uniprot") |> 
    distinct(uniprot, go_id)
  num_proteins <- filtered_ann |> distinct(uniprot) |> nrow()
  log_info(glue(
    "{taxon}: {nrow(filtered_ann)} annotations & {num_proteins} proteins after filter"
  ))
  
  # Export filtered synaptic protein annotation data
  write_csv(filtered_ann, dest, col_names = TRUE)
  write_csv(filtered_ann |> distinct(uniprot), dest_uniprot_only, col_names = TRUE)

}, outdir = outdir)

log_info("Protein filtration script done")