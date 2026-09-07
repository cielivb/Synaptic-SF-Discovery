# Retrieve PDB structure file from AFDB for each synaptic protein in organism
download_organism_structure_data <- function(organism) {
  endpoint <- SPECIES$afdb_endpoint[match(organism, SPECIES$common_name)]
  outpath <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  req_perform(request(endpoint), outpath)
}


# Extract residue-level coordinates and confidences from PDB to CSV
create_residue_site_csv <- function(organism) {
  
  # Extract atom site data for each protein to csv
  afdb_tar_path <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  atom_sites_dir <- SPECIES$atom_sites_dir[match(organism, SPECIES$common_name)]
  system(glue("scripts/extract_atom_sites.sh {afdb_tar_path} {atom_sites_dir}"))
  
  # Reduce atom site data to residue-level data for each protein and make csv
  file_paths <- list.files(atom_sites_dir, full.names = TRUE)
  residue_sites_dir <- SPECIES$residue_sites_dir[match(organism), SPECIES$common_name]
  done <- lapply(file_paths, function(path){
    atom_site_df <- read_csv(path)
    residue_site_df <- atom_site_df |>
      group_by(residue) |>
      summarise(x = mean(x), y = mean(y), z = mean(z), pLDDT = mean(pLDDT))
    write_csv(residue_site_df, glue("{residue_sites_dir}/{basename(path)}"))
  })
  
}


# Create a CSV where each row represents a synaptic protein and contains data
# about that protein's geometry and overall structure confidence
compute_protein_geometries <- function(organism) {
  create_residue_site_csv(organism)
  
  # Compute confidence-weighted eigenvectors and eigenvalues
  
  # Compute geometric descriptors
  # length
  # radius of gyration
  # eccentricity
  # coefficient of variation
  # kurtosis
  # low pLDDT fraction (disorder proxy)
  
  # Export to CSV
  
}