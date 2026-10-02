# This script is for creating miscellaneous/supplementary tables and plots that
# have not already been made

outdir <- "results/10_make_misc"
log_file <- "logs/10_make_misc.log"
source("scripts/setup.R")
configure_logger(log_file)

library(png)
library(tidyr)
library(ggplot2)
