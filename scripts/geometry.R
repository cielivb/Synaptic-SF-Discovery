# Geometric Descriptor Table Construction Operations

# Build Geometry Descriptor tibble
build_geom_tibble <- function(csv_file) {
  protein = basename(csv_file) # UniProt accession
  
  # Load point cloud and get eigenvectors + eigenvalues
  cloud <- load_point_cloud(f)
  eigvecs, eigvals <- get_eigen_properties(cloud)
  
  # Calculate geometric descriptors
  length, width, height <- get_lwh()
  geom_tibble <- tibble(
    protein = protein,
    length = length, width = width, height = height,
    rg = get_rg(), # how far points are from centroid on average
    eccentricity = get_eccentricity() # elongation
    coeff_of_var = get_coeff_of_var() # unevenness of pt distr from centroid
    kurtosis = get_kurtosis() # tailedness
  )
  return(geom_tibble)
}


# Read x, y, z coordinates into tibble
load_point_cloud <- function(file) {
  ...
}


# Get ordered eigenvectors and eigenvalues from x,y,z tibble
get_eigen_properties <- function(cloud) {
  ...
}


# Get length, width, and height via projection
get_lwh <- function(eigvecs, eigvals, cloud) {
  ...
}


# Calculate radius of gyration (compactness indicator)
get_rg <- function(cloud) {
  ...
}


# Calculate eccentricity (elongation indicator)
get_eccentricity <- function(cloud) {
  ...
}


# Calculate coefficient of variation (spread variation indicator)
get_coeff_of_var <- function(cloud) {
  ...
}


# Calculate kurtosis (tailedness indicator)
get_kurtosis <- function(cloud) {
  ...
}