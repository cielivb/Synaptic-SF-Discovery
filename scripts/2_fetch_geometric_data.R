# Get basic geometric descriptions of synaptic proteins across organisms


# Retrieve PDB structure file from AFDB for each synaptic protein in organism
download_organism_structure_data <- function(organism) {
  endpoint <- SPECIES$afdb_endpoint[match(organism, SPECIES$common_name)]
  outpath <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  req_perform(request(endpoint), outpath)
  log_info(glue("Downloaded {organism} AFDB structure data"))
}


# Extract atom site data for each protein to csv
create_atom_site_csvs <- function(organism) {
  afdb_tar_path <- SPECIES$afdb_tar_path[match(organism, SPECIES$common_name)]
  atom_sites_dir <- SPECIES$atom_sites_dir[match(organism, SPECIES$common_name)]
  system(glue("scripts/extract_atom_sites.sh {afdb_tar_path} {atom_sites_dir}"))
  log_info(glue("Created {organism} atom site csv files"))
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


# Get PC scores, eigenvectors, and eigenvalues, weighting by pLDDT
weighted_pca <- function(residue_df) {
  
  # Extract coordinate matrix and weights
  coords <- residue_df |> select(x, y, z)
  weights <- residue_df$pLDDT
  weights <- weights / sum(weights) # Normalise weights
  
  # Centre data using weighted mean
  weighted_mean <- colSums(coords * weights)
  centred_coords <- sweep(coords, 2, weighted_mean)
  
  # Compute weighted covariance matrix
  cov_w <- t(centred_coords) %*% (centred_coords * weights)
  
  # Do an eigendecomposition
  eig <- eigen(cov_w)
  eig_vectors <- eig$vectors # 3x3, each col is a principal axis
  eig_values  <- eig$values
  
  # Compute principal component scores
  pc_scores <- as.matrix(centred_coords) %*% eig_vectors
  
  return(list(
    eigenvectors = eig_vectors,
    eigenvalues = eig_values,
    pc_scores = pc_scores,
    weighted_mean = weighted_mean,
    centred_coords = centred_coords,
    weights = weights
  ))
}


# Get distances to each residue from weighted mean (centre)
get_radial_distances <- function(centred_coords) {
  return(sqrt(rowSums(centred_coords^2)))
}


# Get the length between the most distant residues in a protein
get_length_along_pc1 <- function(pc_scores) {
  pc1 <- pca$pc_scores[, 1]
  return(max(pc1) - min(pc1))
}


# Get descriptor of overall compactness
get_radius_of_gyration <- function(centred_coords, weights) {
  return(sqrt(sum(weights * rowSums(centred_coords^2))))
}


# Get descriptor of rod-like-ness vs sphere-like-ness
get_eccentricity <- function(eigenvalues) {
  eig <- sort(eigenvalues, decreasing = TRUE)
  return(sqrt(1 - (eig[3] / eig[1])))
}


# Get spread of radial distances
get_coefficient_of_variation <- function(radial_distances) {
  return(sd(radial_distances) / mean(radial_distances))
}


# Get descriptor of tailedness vs uniformity
get_kurtosis <- function(radial_distances) {
  d <- radial_distances
  return(mean((d - mean(d))^4) / sd(d)^4)
}


# Return proportion of pLDDT scores less than 50 (disorder proxy)
get_low_pLDDT_fraction <- function(pLDDT_vector) {
  return(mean(pLDDT_vector < 50))
}


# Return proportion of pLDDT scores greater than or equal to 70 (confident core)
get_high_pLDDT_fraction <- function(pLDDT_vector) {
  return(mean(pLDDT_vector >= 70))
}


# Create a CSV where each row represents a synaptic protein and contains data
# about that protein's geometry and overall structure confidence
compute_protein_geometries <- function(organism) {
  
  geom_file <- SPECIES$geom_path[match(organism, SPECIES$common_name)]
  syn_prot_structure_paths <- get_synaptic_protein_structure_paths(organism)
  
  # Write geometric description of each synaptic protein to geom_file CSV
  done <- lapply(syn_prot_structure_paths, function(path) {
    residue_df <- load_residue_level_structure(path)

    # Do confidence-weighted PCA and get radial distances
    wpca <- weighted_pca(residue_df)
    radial_distances <- get_radial_distances(centred_coords)

    # Compute geometric descriptors and append to csv
    geom_df <- tibble(
      uniprot <- sub(".csv", "", basename(path)),
      length = get_length_along_pc1(wpca$pc_scores),
      rg = get_radius_of_gyration(wpca$centred_coords, wpca$weights),
      eccentricity = get_eccentricity(wpca$eigenvalues),
      coeff_var = get_coefficient_of_variation(radial_distances),
      kurtosis = get_kurtosis(radial_distances),
      low_pLDDT = get_low_pLDDT_fraction(residue_df$pLDDT), # disorder proxy
      high_pLDDT = get_high_pLDDT_fraction(residue_df$pLDDT) # confident core
    )
    write_csv(geom_df, geom_file, append = TRUE)
    
  })
  log_info(glue("Computed {organism} protein geometries"))
  
}