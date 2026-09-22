# Extract synaptic protein shape properties.


source("scripts/setup.R")
log_file = "logs/2_extract_shape_data.log"
configure_logger(log_file)

library(R.utils) # unzip gzipped files


# Extract atomic coordinates for each synaptic protein from AFDB .tar archive
extract_atom_site_data <- function(taxon, afdb_taxon_tar, atom_site_csv_dir) {
  # Extract .tar archive
  extract_tar_dir <- glue("results/2_extract_shape_data/{taxon}_extracted_tar")
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
  extract_pdb_dir <- glue("results/2_extract_shape_data/{taxon}_extracted_pdb")
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
    summarise(
      x = mean(x), y = mean(y), z = mean(z), 
      pLDDT = mean(pLDDT), .groups = "drop"
    )
  return(residue_site_df)
}






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


# Get the length between the most distant residues in a protein
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


# Return proportion of pLDDT scores less than 50 (disorder proxy)
get_low_pLDDT_fraction <- function(pLDDT_vector) {
  return(mean(pLDDT_vector < 50))
}


# Return proportion of pLDDT scores greater than or equal to 70 (confident core)
get_high_pLDDT_fraction <- function(pLDDT_vector) {
  return(mean(pLDDT_vector >= 70))
}






# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))

download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
afdb_taxon_meta <- download_meta |> filter(category == "afdb_taxon")




# Extract shape data for each taxon in parallel
log_info("Extracting shape data for each taxon ...")
done <- future_lapply(afdb_taxon_meta$destination, function (afdb_taxon_tar) {
  configure_logger(log_file, worker = TRUE)
  
  taxon <- sub("\\.tar$", "", basename(afdb_taxon_tar)) # e.g., worm
  geom_file <- glue("results/2_extract_shape_data/{taxon}_geom.csv")
  
  # Get .csv file per protein containing atom site coordinate + confidence data
  log_info(glue("Extracting {taxon} atom site data ..."))
  atom_site_csv_dir <- glue("results/2_extract_shape_data/{taxon}_atom_sites_csv")
  extract_atom_site_data(taxon, afdb_taxon_tar, atom_site_csv_dir)
  atom_site_paths <- list.files(atom_site_csv_dir, full.names = TRUE)
  log_info(glue("Extracted {taxon} atom site data"))
  
  
  
  # Write full geometric description of each synaptic protein to geom_file csv
  log_info(glue("Computing {taxon} protein geometry ..."))
  lapply(atom_site_paths, function(path) {
    
    # Load structure data and create evenly weighted weights column, then
    # create geometric descriptor tibble for full structure
    residue_df_f <- load_residue_level_structure(path) |>
      mutate(weights = 1)
    wpca_f <- weighted_pca(residue_df_f)
    radial_distances_f <- get_radial_distances(wpca_f$centred_coords)
    geom_df_f <- tibble(
      uniprot = sub("\\.csv$", "", basename(path)),
      high_pLDDT = get_high_pLDDT_fraction(residue_df_f$pLDDT), # confident core
      low_pLDDT = get_low_pLDDT_fraction(residue_df_f$pLDDT), # disorder proxy
      length_f = get_length_along_pc1(wpca_f$pc_scores),
      rg_f = get_radius_of_gyration(wpca_f$centred_coords, wpca_f$weights),
      eccentricity_f = get_eccentricity(wpca_f$eigenvalues),
      coeff_var_f = get_coefficient_of_variation(radial_distances_f),
      kurtosis_f = get_kurtosis(radial_distances_f)
    )
    
    # Create a high confidence subset and normalise pLDDT into true weights,
    # then create geometric descriptor tibble for core/high confidence subset
    residue_df_hc <- residue_df_f |> 
      filter(pLDDT >= 70)
    if (nrow(residue_df_hc) < 3) { # too few residues for PCA
      geom_df_hc <- tibble(
        length_hc = NA, rg_hc = NA, eccentricity_hc = NA,
        coeff_var_hc = NA, kurtosis_hc = NA
      )
    } else {
      residue_df_hc <- residue_df_hc |> mutate(weights = pLDDT / sum(pLDDT))
      wpca_hc <- weighted_pca(residue_df_hc)
      radial_distances_hc <- get_radial_distances(wpca_hc$centred_coords)
      geom_df_hc <- tibble(
        length_hc = get_length_along_pc1(wpca_hc$pc_scores),
        rg_hc = get_radius_of_gyration(wpca_hc$centred_coords, wpca_hc$weights),
        eccentricity_hc = get_eccentricity(wpca_hc$eigenvalues),
        coeff_var_hc = get_coefficient_of_variation(radial_distances_hc),
        kurtosis_hc = get_kurtosis(radial_distances_hc)
      )
    }
    
    # Merge full tibble with high confidence tibble and append protein geometric
    # description to csv
    geom_df <- bind_cols(geom_df_f, geom_df_hc)
    if (!file.exists(geom_file)) {
      write_csv(geom_df, geom_file, append = FALSE, col_names = TRUE)
    } else {
      write_csv(geom_df, geom_file, append = TRUE, col_names = FALSE)
    }
    
    
  })
  log_info(glue("Computed {taxon} protein geometry"))
  
})


log_info("Extracted shape data for each taxon")
log_info("Shape data extraction script done")
