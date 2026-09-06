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

source("scripts/0_setup.R")
CLUSTER <- makeCluster(detectCores())


# Enter species of interest here
SPECIES <- tibble(
  common_name = c("Human", "Chicken"),
  goa_endpoint = c(
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/CHICKEN/goa_chicken.gaf.gz"
  )
)

SPECIES <- SPECIES |>
  mutate(
    goa_path = file.path("data/goa", basename(goa_endpoint)),
    goa_request = request(goa_endpoint)
  )


source("scripts/1_fetch_synaptic_protein_go_data.R")

stopCluster(CLUSTER)
log_info("End!")

#source("scripts/2_fetch_structure_data.R")
#source("scripts/3_cluster_by_geometry.R")
#source("scripts/4_get_go_by_cluster.R")
#source("scripts/5_visualise_go_by_cluster")
