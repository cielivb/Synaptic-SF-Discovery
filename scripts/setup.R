# Restore packages from lockfile, log execution metadata, create species table, 
# and set up cluster

### MAKE DIRECTORIES -----------------------------------------------------------

system("mkdir -p data/goa results/goa_cleaning")
system("mkdir -p data/afdb_tar data/atom_sites data/geometry")


### RESTORE PACKAGES FROM LOCKFILE ---------------------------------------------

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  print("Installing renv ...")
  install.packages("renv", verbose = FALSE)
  print("Installed renv")
}

# Restore R environment from lockfile
renv::restore(lockfile = "renv.lock", clean = TRUE, prompt = FALSE)

# Load packages
library(dplyr)
library(ggplot2)
library(glue)
library(gridExtra)
library(httr2)
library(jsonlite)
library(logger)
library(parallel)
library(readr)


### SET UP LOGGER AND LOG EXECUTION METADATA -----------------------------------

# Initialise logger
LOGFILE <- "results/analysis.log"
if (file.exists(LOGFILE)) {
  file.remove(LOGFILE)
}
log_appender(appender_console)
log_appender(appender_file(LOGFILE), index = 2)

# Log R version
log_info(glue("Using R version {getRversion()}"))

# Log general metadata
log_info(glue("Running on platform: {R.version$platform}"))
log_info(glue("Locale info: {Sys.getlocale()}"))
log_info(glue("Kernel release: {Sys.info()['release']}"))
log_info(glue("Kernel version: {Sys.info()['version']}"))

# Log environment package metadata
pkgs <- c("dplyr", "ggplot2", "glue", "gridExtra", "httr2", "jsonlite", "logger")
log_info("Environment packages:")
for (pkg in pkgs) {
  log_info(glue("{pkg}: {as.character(packageVersion(pkg))}"))
}


### CREATE GLOBAL SPECIES TABLE ------------------------------------------------

SPECIES <- tibble(
  common_name = c("human", "fruit_fly", "mouse", "rat", "worm", "zebrafish"),
  species = c(
    "Homo sapiens", 
    "Drosophila melanogaster", 
    "Mus musculus", 
    "Rattus norvegicus",
    "Caenorhabditis elegans",
    "Danio rerio"
  ),
  goa_endpoint = c(
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/FLY/goa_fly.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/RAT/goa_rat.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/WORM/goa_worm.gaf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/ZEBRAFISH/goa_zebrafish.gaf.gz"
  ),
  afdb_endpoint = c(
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000005640_9606_HUMAN_v6.tar",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000000803_7227_DROME_v6.tar",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000000589_10090_MOUSE_v6.tar",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000002494_10116_RAT_v6.tar",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000001940_6239_CAEEL_v6.tar",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000000437_7955_DANRE_v6.tar"
  )
)

SPECIES <- SPECIES |>
  mutate(
    goa_path = file.path("data/goa", basename(goa_endpoint)),
    afdb_tar_path = file.path("data/afdb_tar", basename(afdb_endpoint)),
    atom_sites_dir = file.path("data/atom_sites", common_name),
    geom_path = file.path("data/geometry", glue("{common_name}_geometry.csv"))
  )


### SET UP AND START CLUSTER ---------------------------------------------------

CLUSTER <- makeCluster(detectCores())

clusterExport(CLUSTER, "LOGFILE")
clusterEvalQ(CLUSTER, {
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(gridExtra)
  library(httr2)
  library(jsonlite)
  library(logger)
  library(parallel)
  library(readr)
  log_appender(appender_console)
  log_appender(appender_file(LOGFILE))
})

source("scripts/1_fetch_synaptic_protein_go_data.R")
source("scripts/2_fetch_geometric_data.R")

clusterExport(CLUSTER, c(
  "SPECIES",
  
  # 1_fetch_synaptic_protein_go_data.R functions
  "read_go_metadata",
  "get_protein_and_go_counts",
  "keep_synaptic_proteins",
  "filter_organism_go_metadata",
  
  # 2_fetch_geometric_data.R functions
  "create_atom_site_csvs",
  "get_synaptic_protein_structure_paths",
  "load_residue_level_structure",
  "weighted_pca",
  "get_radial_distances",
  "get_length_along_pc1",
  "get_radius_of_gyration",
  "get_eccentricity",
  "get_coefficient_of_variation",
  "get_kurtosis",
  "get_low_pLDDT_fraction",
  "get_high_pLDDT_fraction",
  "compute_protein_geometries"
  
  
))

