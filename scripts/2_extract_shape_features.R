# Extract core synaptic protein shape features.

outdir <- "results/2_extract_shape_features"
log_file = "logs/2_extract_shape_data.log"
source("scripts/setup.R")
configure_logger(log_file)

library(R.utils) # unzip gzipped files
library(ggplot2)
library(GGally)
library(tidyr) # pivot_longer



# Launch cluster for parallel processing
# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))

download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
afdb_taxon_meta <- download_meta |> filter(category == "afdb_taxon")




### Define initial data wrangling functions -------------------------------


# Extract atomic coordinates for each synaptic protein from AFDB .tar archive
extract_atom_site_data <- function(taxon, afdb_taxon_tar, atom_site_csv_dir, outdir) {
  # Extract .tar archive
  extract_tar_dir <- glue("{outdir}/{taxon}_extracted_tar")
  untar(afdb_taxon_tar, exdir = extract_tar_dir)
  
  # Make tibble containing path and uniprot for each synaptic protein in taxon
  synaptic_uniprot <- read_csv(
    glue("results/1_filter_proteins/{taxon}_uniprot_only.csv"), 
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




### Define visualisation functions ----------------------------------------

plot_features_by_taxa_violin <- function(df, title) {
  
  # Rotate df to get df with taxon, feature, and value columns
  df_long <- df |>
    pivot_longer(cols = -taxon, names_to = "feature", values_to = "value")
  
  # Create violin plot
  plt <- ggplot(df_long, aes(x = taxon, y = value, fill = taxon)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point", size = 1.5, color = "black") +
    facet_grid(rows = vars(feature), scales = "free_y") +
    labs(title = title) +
    theme(strip.text.y = element_text(angle = 0))
  
  # Export violin plot
  outfile <- glue("{outdir}/{clean_string(title)}.svg")
  svg(outfile, width = 12, height = 18)
  print(plt)
  invisible(dev.off())

}




### Extract raw candidate shape features -----------------------------------

log_info("Extracting raw candidate shape features ...")

# Extract shape data for each taxon in parallel
done <- future_lapply(afdb_taxon_meta$destination, function (afdb_taxon_tar, outdir) {
  configure_logger(log_file, worker = TRUE)
  
  
  taxon <- sub("\\.tar$", "", basename(afdb_taxon_tar)) # e.g., worm
  shape_file <- glue("{outdir}/{taxon}_shape.csv")
  if (file.exists(shape_file)) { 
    return(NULL) # Skip extraction if features already extracted 
  } 

  # Get .csv file per protein containing atom site coordinate + confidence data
  log_info(glue("Extracting {taxon} atom site data ..."))
  atom_site_csv_dir <- glue("{outdir}/{taxon}_atom_sites_csv")
  extract_atom_site_data(taxon, afdb_taxon_tar, atom_site_csv_dir, outdir)
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
    shape_df <- tibble(
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
    if (!file.exists(shape_file)) {
      write_csv(shape_df, shape_file, append = FALSE, col_names = TRUE)
    } else {
      write_csv(shape_df, shape_file, append = TRUE, col_names = FALSE)
    }
    
  })
  log_info(glue("Computed {taxon} protein shape data"))

}, outdir = outdir)




### Skew-correct and MAD-scale candidate shape feature distributions ------

log_info("Skew-correcting and MAD-scaling candidate shape feature distributions ...")

# Read in shape files for each taxon
taxa <- c("fly", "human", "mouse", "rat", "worm", "zebrafish")
og_shape_df <- bind_rows(lapply(taxa, function (taxon, outdir) {
  shape_file <- glue("{outdir}/{taxon}_shape.csv")
  taxon_df <- read_csv(shape_file, show_col_types = FALSE) |>
    mutate(taxon = taxon)
  return(taxon_df)
}, outdir = outdir))
# Temporarily remove uniprot column
shape_df <- og_shape_df |> select(-uniprot)

# Visualise candidate shape feature distributions before corrections
plot_features_by_taxa_violin(
  shape_df, "Candidate Shape Feature Distributions by Taxa (Raw)"
)

# Correct for skew in candidate shape features via log transform, then visualise.
# Specifically, correct skew in weighted and unweighted kurtosis, length, and rg.
# Attempting to log transform eccentricity makes the eccentricity skew worse.
skewed_cols <- c(
  "kurtosis_w", "kurtosis_u", "length_w", "length_u", "rg_w", "rg_u"
)
shape_df <- shape_df |>
  mutate(across(all_of(skewed_cols), ~ log(.x + 1e-8), .names = "log_{.col}")) |>
  select(-all_of(skewed_cols))
plot_features_by_taxa_violin(
  shape_df, "Candidate Shape Feature Distributions by Taxa (Skew-Corrected)"
)

# Use median-MAD-standardisation on geometric features to normalise
# variance while staying robust against outliers, so that PCA-based variance
# filtering can be used to eliminate low variance features. MAD-scaling gives
# each distribution MAD = 1.
mad_scale <- function(v) {(v - median(v)) / mad(v)}
mad_df <- shape_df |>
  mutate(across(where(is.numeric), ~ mad_scale(.x), .names = "mad_{.col}")) |>
  select(taxon, starts_with("mad_"))
plot_features_by_taxa_violin(mad_df, 
  "Candidate Shape Feature Distributions by Taxa (Skew-Corrected & MAD-Scaled)"
)




### De-correlate candidate shape features to reduce redundancy ------------

log_info("De-correlating candidate shape features ...")

# Visualise correlations between candidate shape features
done <- future_lapply(taxa, function (taxon, outdir) {
  taxon_df <- mad_df |> filter(taxon == !!taxon)
  
  plt <- ggcorr( # default correlation computation = pairwise, pearson
    mad_df |> select(where(is.numeric)), nbreaks = 6, label = TRUE) +
    labs(title = stringr::str_wrap(
      glue("Shape Feature Correlations Before Decorrelation ({taxon})"),
      width = 50)
    ) +
    theme(plot.title = element_text(size = 18, hjust = 0.2))
  
  svg(glue("{outdir}/corr_{taxon}_before.svg", width = 18, height = 16))
  print(plt)
  invisible(dev.off())

  }, outdir = outdir)

# Remove features based on above correlation plot
# Naturally, unweighted and weighted versions of the same metric are highly
# correlated. The differences in their distributions are negligible, so I will
# choose the unweighted versions of these metrics as it is simpler.
# Additionally, there is high correlation (~0.8 R^2) between radius of gyration
# and length. Length is more intuitive to interpret, so I will keep length and
# drop radius of gyration.
mad_df <- mad_df |> select(-ends_with("_w"), -mad_log_rg_u)




### Min-max normalise selected shape features and write to file -----------

# Scale distributions such that minimum values equal 0 and maximum values 
# equal 1.
log_info("Min-max normalising selected shape features ...")
min_max_norm <- function(v) {(v - min(v)) / (max(v) - min(v))}
mad_df <- mad_df |> select(-taxon) |>
  mutate(across(everything(), min_max_norm, .names = "scaled_{.col}"))

# Merge back onto og_shape_df and write to file
shape_df <- bind_cols(og_shape_df, mad_df)
write_csv(shape_df, glue("{outdir}/shape_features.csv"), col_names = TRUE)




log_info("Done!")

