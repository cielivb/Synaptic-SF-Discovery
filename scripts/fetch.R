

##### 1. FILTER ACCESSIONS BY GO METADATA ======================================


# Download GO JSON metadata files for each protein accession

download_go_metadata <- function() {
  
  logger_info("Downloading Homo sapiens GO metadata ...")
  taxon_id <- 9606 # Homo sapiens
  base_endpoint <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
  num_entries_per_page <- 2000
  
  # QuickGO paginates its entries. Determine the number of pages required.
  req <- request(base_endpoint) |>
    req_url_query(taxonId = taxon_id, limit = 1)
  res <- req_perform(req)
  hits <- resp_body_json(res)$numberOfHits
  num_pages <- ceiling(hits / num_entries_per_page)
  pages <- 1:num_pages
  
  # Generate QuickGO URL and result paths for each page
  reqs <- lapply(pages, function(p) {
    request(base_endpoint) |>
      req_url_query(taxonId = taxon_id, limit = num_entries_per_page, page = p) |>
      req_throttle(capacity = 1, fill_time_s = 5)
  })
  paths <- lapply(pages, function(p) {glue("data/go/go_{p}.json")})
  
  # Download GO metadata JSON files
  results <- req_perform_parallel(
    reqs, paths, on_error = "continue", max_active = 10
  )

}


# Merge many small JSON files into one gzipped jsonl file

merge_go_metadata <- function() {
  logger_info("Merging GO metadata ...")
  system("cat data/go/go_*.json > data/all_annotations.jsonl")
  system("gzip data/all_annotations.jsonl")
  system("rm -r data/go/")
}


# Manipulate merged GO JSON into a GO metadata table

wrangle_go_metadata <- function() {
  
  logger_info("Wrangling GO metadata ...")

  # Create tibble containing protein accession, go accession, go type columns
  con <- gzfile("data/all_annotations.jsonl.gz", open = "rt")
  go_data <- stream_in(con)
  close(con)
  
  # Remove proteins that do not have post-synaptic receptor activity (i.e., 
  # remove proteins that do not have GO:0098960 as a molecular function)
  
  # Group tibble by protein accession, and indicate whether each protein is 
  # associated with an excitatory synapse (GO:0060076) or an inhibitory synapse 
  # (GO:0060077).
  
  # Remove proteins that are not associated with excitatory nor inhibitory
  # synapses. It is permissible for a given protein to be associated with both.
  
}


# 1. Shortlist proteins
# Build a base metadata table containing protein accession (str), excit_syn 
# (bool), and inhib_syn (bool) columns. Each entry/row has post-synaptic neuro-
# transmitter receptor activity (GO:0098960), and is associated with either an
# excitatory synapse (GO:0060076) or an inhibitory synapse (GO:0060077).
# Requires the QuickGO API.

filter_by_go_metadata() <- function() {
  
  download_go_metadata()
  merge_go_metadata()
  wrangle_go_metadata()
  
}



