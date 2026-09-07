# Cross-Taxa Discovery of Synaptic Protein Structural Archetypes and Function
# Author: Ciel Baumann
# Last updated: NA

# What structural archetypes exist among synaptic proteins across species, and 
# what biological functions do they correspond to?



source("scripts/setup.R")

# Fetch and clean synaptic protein GO data, creating one CSV per organism 
# mapping synaptic proteins to their GO terms
log_info("Fetching GO Slim synapse data ...")
fetch_go_slim_synapse_data()
log_info("Fetching GO metadata ...")
done <- lapply(SPECIES$common_name, fetch_organism_go_metadata)

# Download synaptic protein structure files and generate geometric descriptors,
# creating a CSV per organism where each row represents the geometric properties
# of a protein in that organism's proteome
log_info("Fetching AFDB predicted structure data ...")
done <- lapply(CLUSTER, SPECIES$common_name, download_organism_structure_data)
log_info("Extracting atom site data ...")
done <- parLapplyLB(CLUSTER, SPECIES$common_name, create_atom_site_csvs)
log_info("Computing protein geometries ...")
done <- parLapplyLB(CLUSTER, SPECIES$common_name, compute_protein_geometries)

# Cluster proteins by geometric properties within each organism

# Discover geometrically similar clusters across species

# TP: visualise geometrically similar clusters

# Calculate functional enrichment for each cluster

# Visualise functional enrichment for each cluster

# TP: determine whether different clusters have statistically different 
# functional enrichment

# Analysis complete!
stopCluster(CLUSTER)
log_info("End!")
