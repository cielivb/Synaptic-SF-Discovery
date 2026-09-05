

# Download bulk Homo sapiens GO metadata file

download_go_metadata <- function() {
  log_info("Downloading Homo sapiens GO metadata")
  endpoint <- "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
  system(glue("curl -o data/goa_human.gaf.gz {endpoint}"))
}


# Read metadata into dataframe containing protein accession and assoc. GO IDs

read_go_metadata <- function() {
  df <- read_tsv(
    "data/goa_human.gaf.gz",
    comment = "!",
    col_names = FALSE,
    col_types = cols_only(
      X2 = col_character(),
      X5 = col_character(),
      X9 = col_character()
    )
  )
  names(df) <- c("uniprot", "go_id", "aspect")
  return(df)
}


# Manipulate merged GO JSON into a GO metadata table

clean_go_metadata <- function(df) {
  
  # Remove proteins that do not have post-synaptic receptor activity (i.e., 
  # remove proteins that do not have GO:0098960 as a molecular function)
  
  # Group tibble by protein accession, and indicate whether each protein is 
  # associated with an excitatory synapse (GO:0060076) or an inhibitory synapse 
  # (GO:0060077).
  
  # Remove proteins that are not associated with excitatory nor inhibitory
  # synapses. It is permissible for a given protein to be associated with both.
  
  write.csv(df, "results/cleaned_go_metadata.csv")
}




# 1.
# Build a base metadata table containing protein accession (str), excit_syn 
# (bool), and inhib_syn (bool) columns. Each entry/row has post-synaptic neuro-
# transmitter receptor activity (GO:0098960), and is associated with either an
# excitatory synapse (GO:0060076) or an inhibitory synapse (GO:0060077).
# Requires the QuickGO API.

fetch_go_metadata <- function() {
  #download_go_metadata()
  df <- read_go_metadata()
  clean_go_metadata(df)
}


fetch_go_metadata()



