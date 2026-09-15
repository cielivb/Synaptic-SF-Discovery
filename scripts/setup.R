# Set up environment

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)


### RESTORE PACKAGES FROM LOCKFILE ---------------------------------------------

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



### CONFIGURE LOGGER -----------------------------------------------------------

configure_logger <- function(logfile) {
  if (file.exists(logfile)) {
    file.remove(logfile)
  }
  log_appender(appender_console)
  log_appender(appender_file(logfile), index = 2)
}



### OTHER UTILITY FUNCTIONS ---------------------------------------------------

clean_string <- function(s) {
  cleaned <- gsub(" ", "_", x) |> gsub("[[:punct:]]", "", _) |> tolower()
  return(cleaned)
}