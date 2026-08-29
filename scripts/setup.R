# Restore packages from lockfile

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  print("Installing renv ...")
  install.packages("renv", verbose = FALSE)
  print("Installed renv")
}

# Restore R environment from lockfile and printing dependency versions
print("Restoring R environment from lockfile ...")
renv::restore(lockfile = "renv.lock", verbose = FALSE)
print("Restored R environment from lockfile")

# Load packages
print("Loading libraries ...")
library(logger)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(gridExtra)
library(glue)
library(parallel)
print("Libraries loaded")

# Initialise logger
print("Initialising logger ...")
logfile <- "results/analysis.log"
if (file.exists(logfile)) {
  file.remove(logfile)
}
log_appender(appender_console)
log_appender(appender_file(logfile), index = 2)
print("Logger initialised")

# Log R version
log_info(glue("Using R version {getRversion()}"))

# Log general metadata
log_info(glue("Running on platform: {R.version$platform}"))
log_info(glue("Locale info: {Sys.getlocale()}"))
log_info(glue("Kernel release: {Sys.info()['release']}"))
log_info(glue("Kernel version: {Sys.info()['version']}"))

# Log environment package metadata
pkgs <- c("logger", "dplyr", "ggplot2", "readr", "stringr", "gridExtra", "glue")
log_info("Environment packages:")
for (pkg in pkgs) {
  log_info(glue("{pkg}: {as.character(packageVersion(pkg))}"))
}