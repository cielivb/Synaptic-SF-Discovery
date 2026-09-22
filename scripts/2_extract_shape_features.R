# Extract core synaptic protein shape features.

outdir <- "results/2_extract_shape_features/"
log_file = "logs/2_extract_shape_data.log"
source("scripts/setup.R")
configure_logger(log_file)

library(R.utils) # unzip gzipped files




# Launch cluster for parallel processing
# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))

download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
afdb_taxon_meta <- download_meta |> filter(category == "afdb_taxon")




### Define initial data wrangling functions -------------------------------


# Extract atomic coordinates for each synaptic protein from AFDB .tar archive
extract_atom_site_data <- function(taxon, afdb_taxon_tar, atom_site_csv_dir) {
  # Extract .tar archive
  extract_tar_dir <- glue("{outdir}/{taxon}_extracted_tar")
  untar(afdb_taxon_tar, exdir = extract_tar_dir)
  
  # Make tibble containing path and uniprot for each synaptic protein in taxon
  synaptic_uniprot <- read_csv(
    glue("results/1_filter_by_func/{taxon}_uniprot_only.csv"), 
    show_col_types = FALSE
  )
  pdb_tibble <- tibble(path = list.files(
    extract_tar_dir, pattern = "\\.pdb\\.gz$", full.names = TRUE
  )) |> 
    mutate(uniprot = strsplit(basename(path), "-") |> sapply(`[[`, 2)) |>
    semi_join(synaptic_uniprot, by = "uniprot")
  
  # Extract synaptic protein .pdb files
  extract_pdb_dir <- glue("{outdir}/{taxon}_extracted_pdb")
  dir.create(extract_pdb_dir, showWarnings = FALSE)
  lapply(pdb_tibble$path, function(gzipped_path) {
    dest <- file.path(extract_pdb_dir, sub("\\.gz$", "", basename(gzipped_path)))
    gunzip(gzipped_path, dest, remove = FALSE)
  })
  
  # Extract atom site data for each synaptic protein in each species
  dir.create(atom_site_csv_dir, showWarnings = FALSE)
  system(glue("scripts/extract_atom_sites.sh {extract_pdb_dir} {atom_site_csv_dir}"))
}


# Load atom site data from csv, reduce to residue-level, and return dataframe
load_residue_level_structure <- function(path) {
  atom_site_df <- read_csv(path, show_col_types = FALSE)
  residue_site_df <- atom_site_df |>
    group_by(residue) |>
    summarise( # Collapse atomic coordinates into representative residue coord
      x = mean(x), y = mean(y), z = mean(z), 
      pLDDT = mean(pLDDT), .groups = "drop"
    )
  return(residue_site_df)
}




### Define weighted PCA and shape descriptor functions --------------------

# Get PC scores, eigenvectors, and eigenvalues, weighting by weights column
weighted_pca <- function(residue_df) {
  
  # Extract coordinate matrix and weights
  coords <- residue_df |> select(x, y, z) |> as.matrix()
  weights <- residue_df$weights
  
  # Centre data using weighted mean
  weighted_mean <- colSums(coords * weights) / sum(weights)
  centred_coords <- sweep(coords, 2, weighted_mean)
  
  # Compute weighted covariance matrix
  cov_w <- t(centred_coords) %*% (centred_coords * weights) / sum(weights)
  
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
    centred_coords = centred_coords,
    weights = weights
  ))
}


# Get distances to each residue from weighted mean (centre)
get_radial_distances <- function(centred_coords) {
  return(sqrt(rowSums(centred_coords^2)))
}


# Get the length between the most distant residues in a protein (maximum span)
get_length_along_pc1 <- function(pc_scores) {
  pc1 <- pc_scores[, 1]
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




### Extract candidate shape features --------------------------------------


# Extract shape data for each taxon in parallel
done <- future_lapply(afdb_taxon_meta$destination, function (afdb_taxon_tar) {
  configure_logger(log_file, worker = TRUE)
  
  
  taxon <- sub("\\.tar$", "", basename(afdb_taxon_tar)) # e.g., worm
  geom_file <- glue("{outdir}/{taxon}_geom.csv")
  
  # Get .csv file per protein containing atom site coordinate + confidence data
  log_info(glue("Extracting {taxon} atom site data ..."))
  atom_site_csv_dir <- glue("{outdir}/{taxon}_atom_sites_csv")
  extract_atom_site_data(taxon, afdb_taxon_tar, atom_site_csv_dir)
  atom_site_paths <- list.files(atom_site_csv_dir, full.names = TRUE)
  log_info(glue("Extracted {taxon} atom site data"))
  
  
  # Compute values for candidate shape features for each protein
  log_info("Computing {taxon} protein shape data ...")
  lapply(atom_site_paths, function(path) {
    
    # Load structure data and run both weighted and unweighted PCA.
    # I am trying both weighted and unweighted PCA as I do not yet know which 
    # will produce the features that explain the most geometric variance in the
    # dataset.
    residue_df <- load_residue_level_structure(path)
    weighted_residue_df <- residue_df |> mutate(weights = pLDDT)
    unweighted_residue_df <- residue_df |> mutate(weights = 1)
    pca_w <- weighted_pca(weighted_residue_df)
    pca_u <- weighted_pca(unweighted_residue_df)
    
    # The below radial distances differ due to different weighted means
    radial_distances_w <- get_radial_distances(pca_w$centred_coords)
    radial_distances_u <- get_radial_distances(pca_u$centred_coords)
    
    # Create tibble of protein shape metadata
    geom_df <- tibble(
      uniprot = sub("\\.csv$", "", basename(path)),
      mean_pLDDT = mean(residue_df$pLDDT),
      length_u = get_length_along_pc1(pca_u$pc_scores),
      length_w = get_length_along_pc1(pca_w$pc_scores),
      rg_u = get_radius_of_gyration(pca_u$centred_coords, pca_u$weights),
      rg_w = get_radius_of_gyration(pca_w$centred_coords, pca_w$weights),
      eccentricity_u = get_eccentricity(pca_u$eigenvalues),
      eccentricity_w = get_eccentricity(pca_w$eigenvalues),
      coeff_var_u = get_coefficient_of_variation(radial_distances_u),
      coeff_var_w = get_coefficient_of_variation(radial_distances_w),
      kurtosis_u = get_kurtosis(radial_distances_u),
      kurtosis_w = get_kurtosis(radial_distances_w)
    )
    
    # Append protein shape metadata to csv
    if (!file.exists(geom_file)) {
      write_csv(geom_df, geom_file, append = FALSE, col_names = TRUE)
    } else {
      write_csv(geom_df, geom_file, append = TRUE, col_names = FALSE)
    }
    
  })
  log_info(glue("Computed {taxon} protein shape data"))

})
