# Set up environment


### Restore packages from lockfile ----------------------------------------

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  print("Installing renv ...")
  install.packages("renv", verbose = FALSE)
  print("Installed renv")
}

# Restore R environment from lockfile
renv::restore(lockfile = "renv.lock", clean = TRUE, prompt = FALSE)

# Load packages used in most scripts
library(dplyr)
library(future)
library(future.apply)
library(glue)
library(logger)
library(readr)

options(dplyr.summarise.inform = FALSE)


### Create project directories --------------------------------------------

dir.create("logs", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
result_dirs <- list(
  "1_filter_proteins", "2_extract_shape_features", 
  "3_extract_functional_features", "4_visualise_manifolds", 
  "5_cluster_and_test_cross_taxa_sf_manifold",
  "6_characterise_clusters", "7_bonus_analyses"
)
lapply(result_dirs, function (d) {
  dir.create(glue("results/{d}"), showWarnings = FALSE)
})


### Define utility functions ----------------------------------------------

configure_logger <- function(logfile, worker=FALSE) {
  if (file.exists(logfile) & !worker) {
    file.remove(logfile)
  }
  log_appender(appender_console)
  log_appender(appender_file(logfile), index = 2)
}

# Remove punctuation and whitespace etc for consistent filename formatting
clean_string <- function(s) {
  cleaned <- s |>
    gsub(" ", "_", x = _) |> 
    gsub("[\\%()&#:/.,-]", "", x = _) |> 
    tolower()
  return(cleaned)
}
