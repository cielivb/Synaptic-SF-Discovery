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
  system(glue("scripts/extract_atom_sites.sh {afdb_tar_path} {organism}"))
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
    mutate(uniprot = strsplit(basename(files), "-") |> sapply(`[[`, 2)) |>
    filter(uniprot %in% synaptic_proteins$uniprot)
  return(as.list(df$files))
  
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


# Create a CSV where each row represents a synaptic protein and contains data
# about that protein's geometry and overall structure confidence
compute_protein_geometries <- function(organism) {
  
  geom_file <- SPECIES$geom_path[match(organism, SPECIES$common_name)]
  syn_prot_structure_paths <- get_synaptic_protein_structure_paths(organism)

  # Write geometric description of each synaptic protein to geom_file CSV
  lapply(syn_prot_structure_paths, function(path) {
    
    # Load structure data and create evenly weighted weights column, then
    # create geometric descriptor tibble for full structure
    residue_df_f <- load_residue_level_structure(path) |>
      mutate(weights = 1)
    wpca_f <- weighted_pca(residue_df_f)
    radial_distances_f <- get_radial_distances(wpca_f$centred_coords)
    geom_df_f <- tibble(
      uniprot = sub(".csv", "", basename(path)),
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
    
    # Merge full structure geometric descriptor tibble with subset tibble,
    # then compute ratios and write to csv
    geom_df <- bind_cols(geom_df_f, geom_df_hc) |>
      mutate(
        length_ratio = ifelse(is.na(length_hc), NA, length_f / length_hc),
        rg_ratio = ifelse(is.na(rg_hc), NA, rg_f / rg_hc),
        eccentricity_ratio = ifelse(
          is.na(eccentricity_hc), NA, eccentricity_f / eccentricity_hc
        ),
        coeff_var_ratio = ifelse(
          is.na(coeff_var_hc), NA, coeff_var_f / coeff_var_hc
        ),
        kurtosis_ratio = ifelse(is.na(kurtosis_hc), NA, kurtosis_f / kurtosis_hc)
      )
    write_csv(geom_df, geom_file, append = TRUE)
    
  })
  
  log_info(glue("Computed {organism} protein geometries"))
  
}