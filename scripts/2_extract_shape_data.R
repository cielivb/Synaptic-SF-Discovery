# Extract synaptic protein shape properties.


source("scripts/setup.R")
configure_logger("logs/extract_shape_data.log")


# Use 1 worker per taxon, or if less than 6 cores available, just use one
# worker per core
available_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = min(available_cores, 6))


# Extract shape data for each taxon in parallel
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)
afdb_taxon_meta <- download_meta |> filter(category == "afdb_taxon")
future_lapply(afdb_taxon_meta$dest, function (afdb_taxon_tar) {
  
  taxon <- sub(".tar", "", basename(afdb_taxon_tar)) # e.g., worm
  
  
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
  unzip(files = pdb_tibble$path, exdir = extract_pdb_dir)

  # Extract atom site data for each synaptic protein in each species
  atom_site_csv_dir <- glue("results/2_extract_shape_data/{taxon}_atom_sites_csv")
  dir.create(atom_site_csv_dir, showWarnings = FALSE)
  system("scripts/extract_atom_sites.sh {extract_pdb_dir} {atom_site_csv_dir}")
  
})
