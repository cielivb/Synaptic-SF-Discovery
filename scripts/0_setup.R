# Restore packages from lockfile and log execution metadata

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

# Initialise logger
logfile <- "results/analysis.log"
if (file.exists(logfile)) {
  file.remove(logfile)
}
log_appender(appender_console)
log_appender(appender_file(logfile), index = 2)

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