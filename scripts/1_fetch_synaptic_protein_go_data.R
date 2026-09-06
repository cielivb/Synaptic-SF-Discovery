# Map synaptic proteins to their GO terms across organism


# Download and clean GO slim synapse data
fetch_go_slim_synapse_data <- function() {
  goslim_path <- "data/goslim_synapse.tsv"
  endpoint <- "http://release.geneontology.org/2026-08-05/ontology/subsets/goslim_synapse.tsv"
  request(endpoint) |> req_perform(goslim_path)
  
  raw_slim <- read_tsv(goslim_path, col_names = TRUE, show_col_types = FALSE)
  
  clean_slim <- raw_slim |>
    select(go_id = `?x`) |>
    mutate(go_id = sub("_", ":", basename(gsub("[<>]", "", go_id)))) |>
    arrange(go_id)
  
  write_csv(clean_slim, "results/goa_cleaning/goslim_synapse_clean.csv")
}


# Download GO IDs for each protein in organism
download_go_metadata <- function(organism) {
  goa_request <- SPECIES$goa_request[match(organism, SPECIES$organism)]
  goa_path <- SPEICES$goa_path[match(organism, SPECIES$organism)]
  req_perform(goa_request, goa_path)
}


# Make table of total number of proteins and total number of GO annotations
get_protein_and_go_counts <- function(df) {
  counts <- df |>
    summarise(
      total_proteins = n_distinct(uniprot),
      total_annotations = n()
    )
  return(counts)
}


# Read protein accessions and GO IDs into tibble
read_go_metadata <- function(goa_path) {
  df <- read_tsv(
    goa_path, comment = "!", col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession 
      X5 = col_character() # GO ID
    ),
    show_col_types = FALSE
  )
  names(df) <- c("uniprot", "go_id")
  return(df)
}


# Keep entries only where the protein has at least one synapse-associated GO ID
keep_synaptic_proteins <- function(df) {
  df_filtered <- df |>
    group_by(uniprot) |>
    filter(any(!is.na(match(go_id, SYNAPSE_GO_IDS$go_id)))) |>
    ungroup()
  return(df_filtered)
}


# Export csv of synaptic proteins and associated GO IDs
fetch_organism_go_metadata <- function(organism) {
  
  # Download and read in unfiltered GO annotations
  download_go_metadata(organism)
  raw_df <- read_go_metadata(goa_path)
  
  # Export number of proteins and total number of GO annotations
  init_counts <- get_protein_and_go_counts(raw_df)
  write_csv(init_counts, glue("results/goa_cleaning/{organism}_0_init_counts.csv"))

  # Filter to keep only synaptic proteins
  filtered_df <- keep_synaptic_proteins(raw_df)
  
  # Export number of proteins and total number of GO annotations again
  filtered_counts <- get_protein_and_go_counts(filtered_df)
  write_csv(filtered_counts, glue("results/goa_cleaning/{organism}_1_filtered_counts.csv"))
  
  # Export filtered dataframe
  write_csv(filtered_df, glue("resuts/goa_clean/{organism}.csv"))
  
}


# Create CSV per organism mapping synaptic proteins to their GO terms
system("mkdir -p data/goa results/goa_cleaning results/goa_clean")
log_info("Fetching GO Slim synapse data ...")
fetch_go_slim_synapse_data()
SYNAPSE_GO_IDS <- read_csv("results/goa_cleaning/goslim_synapse_clean.csv")
log_info("Fetching GO metadata (this may take a while) ...")
parLapplyLB(CLUSTER, as.vector(SPECIES$common_name), fetch_organism_go_metadata)
remove(SYNAPSE_GO_IDS)
