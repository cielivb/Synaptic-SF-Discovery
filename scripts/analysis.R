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

# Fetch and clean synaptic protein GO data, creating one CSV per organism 
# mapping synaptic proteins to their GO terms
log_info("Fetching GO Slim synapse data ...")
fetch_go_slim_synapse_data()
log_info("Fetching GO metadata ...")
done <- parLapplyLB(CLUSTER, SPECIES$common_name, fetch_organism_go_metadata)


# Analysis complete!
stopCluster(CLUSTER)
log_info("End!")
