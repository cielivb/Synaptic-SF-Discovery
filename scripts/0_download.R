# Download all files specified in data/downloads.csv
# Most files belong on same connection- and rate-limited server, hence
# must be patient and download files sequentially.

source("scripts/setup.R")
configure_logger("logs/0_download.log")


log_info("Downloading multi-GB dataset - this will take a while")
download_meta <- read_csv("data/downloads.csv", show_col_types = FALSE)

for (i in seq_len(nrow(download_meta))) {
  
  endpoint <- download_meta$endpoint[i]
  dest <- download_meta$destination[i]
  size_c <- download_meta$size_compressed[i]
  size_u <- download_meta$size_uncompressed[i]
  download_size <- ifelse(!is.na(size_c), size_c, size_u)
  
  log_info(glue("Downloading {endpoint} into {dest} ({download_size}) ..."))
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  system(glue("wget --tries=10 --continue -O {dest} {endpoint}"))
  downloaded_bytes <- file.size(dest)
  downloaded_mib <- downloaded_bytes / 1024^2
  downloaded_gib <- downloaded_bytes / 1024^3
  
  log_info(glue(
    "Downloaded {round(downloaded_gib, 2)} GiB ({round(downloaded_mib, 2)} MiB) into {dest}"
  ))

}