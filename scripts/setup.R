# Restore packages from lockfile

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Restore R environment from lockfile and printing dependency versions
renv::restore(lockfile = "../renv.lock")

# Reset logfile
logfile <- "results/analysis.log"
if (file.exists(logfile)) {
  file.remove(logfile)
}

# Initialise logger
library(logger)
log_appender(appender_console)
log_appender(appender_file(logfile), index = 2)

# Log R version
log_info(sprintf("Using R version %s", getRversion()))

# Log general metadata
log_info(sprintf("Running on platform: %s", R.version$platform))
log_info(sprintf("Locale info: %s", Sys.getlocale()))
log_info(sprintf("Kernel release: %s", Sys.info()["release"]))
log_info(sprintf("Kernel version: %s", Sys.info()["version"]))

# Log environment package metadata
log_info("Environment packages:")
for (pkg in c("tidyverse", "logger", "kableExtra")) {
  log_info(sprintf("%s: %s", pkg, as.character(packageVersion(pkg))))
}