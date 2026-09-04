# Influence of Synapse Type on Receptor Protein Geometry in Homo Sapiens
# Author: Ciel Baumann
# Last updated: NA

# Both reviewed and unreviewed protein accessions are used to increase power.
# Using reviewed proteins alone introduces human selection and experimental
# accessibility biases, and greatly reduces statistical power. Combining both
# reviewed and unreviewed proteins mitigates such biases and increases statist-
# ical power, but at the risk of including fragments and lower-confidence and
# not yet experimentally validated structures. Later steps reduce the severity
# of this risk through strict GO-annotation filtering, outlier stripping, 
# and confidence-aware analysis.

source("scripts/setup.R")

# 0. Download .tsv of all protein accessions in the human proteome.
proteome <- "UP000005640"
log_info(glue("Checking proteome accession {proteome} is Homo sapiens ..."))
req <- request("https://rest.uniprot.org/proteomes/UP000005640") |>
  req_url_query(fields = "organism_name")
species <- resp_body_json(req_perform(req))$organism$scientificName
if (species != "Homo sapiens") {
  log_failure(glue("Expected Homo sapiens, got {species}"))
}

log_info("Downloading protein accessions (may take a few minutes) ...")
req <- request("https://rest.uniprot.org/uniprotkb/stream") |>
  req_url_query(
    query  = "proteome:UP000005640",
    fields = "accession",
    format = "tsv"
  )
res <- req_perform(req, path="data/accessions.tsv")
if (res$status != 200) {
  log_failure(glue("Failed to download protein accessions ({res$status})"))
}