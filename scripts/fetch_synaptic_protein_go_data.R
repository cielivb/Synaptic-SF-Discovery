# Map synaptic proteins to their GO terms across organism


# Download and clean GO slim synapse data
fetch_go_slim_synapse_data <- function() {
  goslim_path <- "data/goslim_synapse.tsv"
  endpoint <- "http://release.geneontology.org/2026-08-05/ontology/subsets/goslim_synapse.tsv"
  request(endpoint) |> req_perform(goslim_path)
  
  raw_slim <- read_tsv(goslim_path, col_names = TRUE)
  
  clean_slim <- raw_slim |>
    select(go_id = `?x`) |>
    mutate(
      go_id = go_id |> 
        gsub("[<>]", "") |>
        basename() |>
        sub("_", ":")
    ) |>
    arrange(go_id)
  
  write_csv(clean_slim, "results/goa_cleaning/goslim_synapse_clean.csv")
}


# Download bulk GO metadata files
download_go_metadata <- function() {
  req_perform_parallel(
    as.list(SPECIES$goa_request),
    as.vector(SPECIES$goa_path),
    on_error="continue"
  )
}


# Read protein accessions and GO IDs into tibble
read_go_metadata <- function(goa_path) {
  df <- read_tsv(
    goa_path, comment = "!", col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession 
      X5 = col_character() # GO ID
    )
  )
  names(df) <- c("uniprot", "go_id")
  return(df)
}


# Keep entries only where the protein has at least one synapse-associated GO ID
keep_synaptic_proteins <- function(df) {
  synapse_go_ids <- read_csv("results/goa_cleaning/goslim_synapse_clean.csv")
  df_filtered <- df |>
    group_by(uniprot) |>
    filter(any(!is.na(match(go_id, synapse_go_ids$go_id)))) |>
    ungroup()
  return(df_filtered)
}


# Export csv of synaptic proteins and associated GO IDs
clean_go_metadata <- function(goa_path) {
  raw_df <- read_go_metadata(goa_path)
  organism <- SPECIES$common_name[match(goa_path, SPECIES$goa_path)]
  
  # Export number of proteins and total number of GO annotations
  init_counts <- raw_df |>
    summarise(
      total_proteins = n_distinct(uniprot),
      total_annotations = n()
    )
  write_csv(init_counts, glue("results/goa_cleaning/{organism}_0_init_counts.csv"))

  # Filter to keep only synaptic proteins
  filtered_df <- keep_synaptic_proteins(raw_df)
  
  # Export number of proteins and total number of GO annotations again
  filtered_counts <- filtered_df |>
    summarise(
      total_proteins = n_distinct(uniprot),
      total_annotations = n()
    )
  write_csv(init_counts, glue("results/goa_cleaning/{organism}_1_filtered_counts.csv"))
  
  # Export filtered dataframe
  write_csv(filtered_df, glue("resuts/goa_clean/{organism}.csv"))
}


# Create CSV per organism mapping synaptic proteins to their GO terms
fetch_synaptic_protein_go_data <- {
  system("mkdir -p data/goa results/goa_cleaning results/goa_clean")
  fetch_go_slim_synapse_data()
  download_go_metadata()
  parLapplyLB(CLUSTER, as.vector(SPECIES$goa_path), clean_go_metadata)
}