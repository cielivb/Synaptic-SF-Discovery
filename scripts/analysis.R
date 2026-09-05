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
source("scripts/fetch_go_metadata.R")
source("scripts/fetch_structure_data.R")



# 