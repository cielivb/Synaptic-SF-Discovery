# Get basic geometric descriptions of synaptic proteins across organisms


# Retrieve PDB structure file from AFDB for each synaptic protein in organism
download_organism_structure_data <- function(organism) {
  endpoint <- SPECIES$afdb_endpoint[match(organism, SPECIES$common_name)]
  outpath <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  req_perform(request(endpoint), outpath)
}


# Extract atom site data for each protein to csv
create_atom_site_csvs <- function(organism) {
  afdb_tar_path <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  atom_sites_dir <- SPECIES$atom_sites_dir[match(organism, SPECIES$common_name)]
  system(glue("scripts/extract_atom_sites.sh {afdb_tar_path} {atom_sites_dir}"))
}


# Get a list containing synaptic protein structure file paths
get_synaptic_protein_structure_paths <- function(organism) {
  
  # Get synaptic proteins in organism
  go_df <- read_csv(glue("data/goa/{organism}_clean.csv"), show_col_types = FALSE)
  synaptic_proteins <- go_df |> distinct(uniprot)
  
  # Get Uniprot accessions from AFDB structure data file names and keep only
  # those accessions and filepaths associated with synaptic proteins
  atom_sites_dir <- SPECIES$atom_sites_dir[match(organism, SPECIES$common_name)]
  df <- tibble(files = list.files(atom_sites_dir, full.names = TRUE))
  df <- df |>
    mutate(uniprot = sub(".csv", "", basename(files))) |>
    filter(!is.na(match(uniprot, synaptic_proteins$uniprot)))
  return(as.list(df$files))
  
}


# Load atom site data from csv, reduce to residue-level, and return dataframe
load_residue_level_structure <- function(path) {
  atom_site_df <- read_csv(path, show_col_types = FALSE)
  residue_site_df <- atom_site_df |>
    group_by(residue) |>
    summarise(x = mean(x), y = mean(y), z = mean(z), pLDDT = mean(pLDDT))
  return(residue_site_df)
}


# Create a CSV where each row represents a synaptic protein and contains data
# about that protein's geometry and overall structure confidence
compute_protein_geometries <- function(organism) {
  
  geom_file <- SPECIES$geom_path[match(organism, SPECIES$common_name)]
  syn_prot_structure_paths <- get_synaptic_protein_structure_paths(organism)
  
  done <- lapply(syn_prot_structure_paths, function(path) {
    
    residue_df <- load_residue_level_structure(path)
    geom_df <- residue_df
    
    # Compute confidence-weighted eigenvectors and eigenvalues
    
    # Compute geometric descriptors
    # length
    # radius of gyration
    # eccentricity
    # coefficient of variation
    # kurtosis
    # low pLDDT fraction (disorder proxy)
    
    # Append line to CSV
    write_csv(geom_df, geom_file, append = TRUE)
    
  })
  
}