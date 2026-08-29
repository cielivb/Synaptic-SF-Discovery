# Restore packages from lockfile

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Restore R environment from lockfile and printing dependency versions
renv::restore(lockfile = "../renv.lock")

# Load packages
library(logger)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(gridExtra)
library(glue)

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
pkgs <- c("logger", "dplyr", "ggplot2", "readr", "stringr", "gridExtra", "glue")
log_info("Environment packages:")
for (pkg in pkgs) {
  log_info(glue("{pkg}: {as.character(packageVersion(pkg))}"))
}