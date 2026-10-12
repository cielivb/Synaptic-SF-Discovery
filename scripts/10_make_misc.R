# This script is for creating miscellaneous/supplementary tables and plots that
# have not already been made

outdir <- "results/10_make_misc"
log_file <- "logs/10_make_misc.log"
source("scripts/setup.R")
configure_logger(log_file)
log_info("Running 10_make_misc.R ...")

library(png)
library(tidyr)
library(ggplot2)
library(scatterplot3d)
library(grid)
library(gridExtra)


# Create stacked histogram of proteins per cluster, coloured by taxon
cluster_counts <- read_csv(
  "results/04_cluster_by_gms/cluster_counts.csv",
  show_col_types = FALSE
) |> pull(n_proteins)
svg(glue("{outdir}/proteins_per_gms_cluster.svg"), height = 8, width = 8)
hist(
  cluster_counts, 
  breaks = 5,
  xlim = c(min(cluster_counts), max(cluster_counts)),
  xlab = "Number of Proteins", ylab = "Number of GMS Clusters",
  main = "Distribution of Protein Counts Per GMS Cluster",
)
invisible(dev.off())



# Create shape feature visual intuition plots
shape_features_df <- read_csv(
  "results/02_extract_shape_features/shape_features.csv",
  show_col_types = FALSE
)
features <- c("length_u", "mean_pLDDT", "kurtosis_u", "coeff_var_u", "eccentricity_u")
lapply(features, function (feature) {
  
  # Get uniprot accessions for extreme values of feature
  by_feature <- shape_features_df |> arrange(!!sym(feature)) |> select(uniprot, taxon) |> distinct(uniprot, taxon)
  low_feature <- head(by_feature, 3)
  high_feature <- tail(by_feature, 3)
  
  # Create and save a plot for each uniprot
  make_plot <- function (uniprot, taxon, feature, direction) {
    coords_file <- glue("results/02_extract_shape_features/{taxon}_atom_sites_csv/{uniprot}.csv")
    
    # Load residue-level coordinate data for uniprot protein
    # This code block taken from load_residue_level_data in 02_extract_shape_Features.R
    atom_site_df <- read_csv(coords_file, show_col_types = FALSE)
    residue_site_df <- atom_site_df |>
      group_by(residue) |>
      summarise( # Collapse atomic coordinates into representative residue coord
        x = mean(x), y = mean(y), z = mean(z), 
        pLDDT = mean(pLDDT), .groups = "drop"
      )

    # Create 3D plot and save
    png(filename=glue("{outdir}/{direction}_{feature}_{uniprot}.png"), 
        height = 8, width = 8, res = 300, units = "in")
    scatterplot3d(
      residue_site_df$x,
      residue_site_df$y,
      residue_site_df$z,
      color = residue_site_df$pLDDT,
      pch = 16
    )
    invisible(dev.off())
  }
  lapply(seq_len(nrow(low_feature)), function(i) { # low feature value proteins
    make_plot(
      uniprot = low_feature$uniprot[i],
      taxon = low_feature$taxon[i],
      feature = feature,
      direction = "low"
    )
  })
  lapply(seq_len(nrow(high_feature)), function(i) { # high feature value proteins
    make_plot(
      uniprot = high_feature$uniprot[i],
      taxon = high_feature$taxon[i],
      feature = feature,
      direction = "high"
    )
  })
  
  # Combine feature plots into single row for each extreme
  combine_plots <- function (uniprots, feature, direction) {
    imgs <- lapply(uniprots, function (uniprot) {
      img <- readPNG(glue("{outdir}/{direction}_{feature}_{uniprot}.png"))
      rasterGrob(img, interpolate = TRUE)
    })
    png(
      filename = glue("{outdir}/{direction}_{feature}_row.png"),
      width = 12, height = 4, res = 300, units = "in"
    )
    grid.arrange(grobs = imgs, ncol = 3)
    invisible(dev.off())
  }
  combine_plots(low_feature$uniprot, feature, "low")
  combine_plots(high_feature$uniprot, feature, "high")
})





log_info("Ended 10_make_misc.R ...")
