# Retrieve PDB structure file from AFDB for each synaptic protein in organism
download_organism_structure_data <- function(organism) {
  endpoint <- SPECIES$afdb_endpoint[match(organism, SPECIES$common_name)]
  outpath <- SPECIES$afdb_path[match(organism, SPECIES$common_name)]
  req_perform(request(endpoint), outpath)
}


# Create a CSV where each row represents a synaptic protein and contains data
# about that protein's geometry and overall structure confidence
compute_protein_geometries <- function(organism) {
  
  # Extract residue-level coordinates and confidences from PDB to CSV
  
  # Compute confidence-weighted eigenvectors and eigenvalues
  
  # Compute other geometric descriptors
  
  # Export to CSV
  
}