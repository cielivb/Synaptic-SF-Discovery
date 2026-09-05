

# Download bulk Homo sapiens GO metadata file

download_go_metadata <- function() {
  log_info("Downloading Homo sapiens GO metadata ...")
  endpoint <- "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
  system(glue("curl -o data/goa_human.gaf.gz {endpoint}"))
}


# Read metadata into dataframe containing protein accession and assoc. GO IDs

read_go_metadata <- function() {
  log_info("Reading GO metadata ...")
  df <- read_tsv(
    "data/goa_human.gaf.gz",
    comment = "!",
    col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(), # UniProt accession
      X5 = col_character(), # GO ID
    )
  )
  names(df) <- c("uniprot", "go_id")
  return(df)
}


# Manipulate merged GO JSON into a GO metadata table

clean_go_metadata <- function(df) {
  log_info("Cleaning GO metadata ...")
  
  # Export number of proteins and total number of GO annotations
  init_counts <- df |>
    summarise(
      total_proteins = n_distinct(uniprot),
      total_annotations = n()
    )
  write_csv(init_counts, "results/go_cleaning/0_init_counts.csv")
  
  # Define GO IDs
  psra_go <- "GO:0098960" # post-synaptic receptor activity
  excit_go <- "GO:0060076" # in excitatory synapse
  inhib_go <- "GO:0060077" # in inhibitory synapse
  
  # Aggregate tibble, indicating whether each protein has post-synaptic receptor
  # activity, and whether each is in an excitatory or inhibitory synapse
  agg_df <- df |>
    group_by(uniprot) |>
    summarise(
      psra = any(go_id == psra_go),
      excit = any(go_id == excit_go),
      inhib = any(go_id == inhib_go),
      .groups = "drop"
    )
  
  # Export percentages of TRUE values per boolean column in tibble
  pct_true <- agg_df |>
    summarise(
      postsyn_recept_pct = mean(psra),
      excit_syn_pct = mean(excit),
      inhib_syn_pct = mean(inhib)
    )
  write_csv(pct_true, "results/go_cleaning/1_pct_true_before_filter.csv")
  
  # Keep only proteins that have post-synaptic receptor activity AND are 
  # associated with (an excitatory synapse or an inhibitory synapse). It is 
  # permissible for a given protein to be associated with both synapse types.
  filtered_df <- agg_df |> filter(psra, excit | inhib)
  write.csv(filtered_df, "results/cleaned_go_metadata.csv")
  
  # Export percentages of TRUE values per bool col after filtering
  pct_true_2 <- filtered_df |>
    summarise(
      postsyn_recept_pct = mean(psra),
      excit_syn_pct = mean(excit),
      inhib_syn_pct = mean(inhib)
    )
  write_csv(pct_true_2, "results/go_cleaning/2_pct_true_after_filter.csv")
  
  # Export synapse exclusivity overview, i.e., percentage of proteins associated
  # with excitatory and not inhibitory synapses, or both excitatory and inhibi-
  # tory synapse, etc.
  df_exclusivity <- filtered_df |>
    mutate(
      excit_true_inhib_false = excit & !inhib,
      excit_false_inhib_true = !excit & inhib,
      excit_false_inhib_false = !excit & !inhib,
      excit_true_inhib_true = excit & inhib
    ) |>
    summarise(
      excit_true_inhib_false_pct = mean(excit_true_inhib_false),
      excit_false_inhib_true = mean(excit_false_inhib_true),
      excit_false_inhib_false = mean(excit_false_inhib_false),
      excit_true_inhib_true = mean(excit_true_inhib_true)
    )
  write_csv(df_exclusivity, "results/go_cleaning/3_pct_exclusive.csv")
  
}





# Build a base metadata table containing protein accession (str), excit_syn 
# (bool), and inhib_syn (bool) columns. Each entry/row has post-synaptic neuro-
# transmitter receptor activity (GO:0098960), and is associated with either an
# excitatory synapse (GO:0060076) or an inhibitory synapse (GO:0060077).

fetch_go_metadata <- function() {
  #download_go_metadata()
  df <- read_go_metadata()
  clean_go_metadata(df)
}


fetch_go_metadata()



