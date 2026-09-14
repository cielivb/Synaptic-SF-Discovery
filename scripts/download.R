# Download all files specified in data/downloads.csv
# Most files belong on same connection- and rate-limited server, hence
# must be patient and download files sequentially.

source("scripts/setup.R")
configure_logger("logs/download.log")

library(readr)
library(httr2)

log_info("Downloading multi-GB dataset - this will take a while")
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)

for (i in seq_len(nrow(download_meta))) {
  
  endpoint <- download_meta$endpoint[i]
  dest <- download_meta$destination[i]
  size_c <- download_meta$size_compressed[i]
  size_u <- download_meta$size_uncompressed[i]
  download_size <- ifelse(!is.na(size_c), size_c, size_u)
  
  req <- request(endpoint) |> req_retry(max_tries = 5)
  log_info(glue("Downloading {endpoint} into {dest} ({download_size}) ..."))
  req_perform(req, path = dest)
  downloaded_bytes <- file.size(dest)
  downloaded_mb <- downloaded_bytes / 1000000
  downloaded_gb <- downloaded_mb / 1000
  log_info(glue(
    "Downloaded {round(downloaded_gb, 2)} GB ({round(downloaded_mb, 2)} MB) into {dest}"
  ))

}