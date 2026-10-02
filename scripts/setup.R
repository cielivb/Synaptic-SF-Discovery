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
  "01_filter_proteins", "02_extract_shape_features", 
  "03_extract_go_module_scores", "04_cluster_by_gms",
  "05_describe_gms_cluster_func", "06_make_taxa_plots",
  "07_make_shape_plots", "08_core_analysis",
  "09_make_gms_cluster_summaries", "10_make_misc"
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
