# Draw one parallel coordinates plot per GMS cluster, illustrating distribution of
# shape features in addition to individual lines representing individual proteins
# in that GMS cluster. 
# Additionally, create a grid plot containing the above plots for ease of comparison.

outdir <- "results/7_make_shape_plots"
log_file <- "logs/7_make_shape_plots.log"
source("scripts/setup.R")
configure_logger(log_file)
plan(multisession, workers = parallel::detectCores(logical = FALSE))


library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(GGally)




# Create parallel coordinates plot showing shape features of proteins in a GMS 
# cluster, where lines are styled by taxon
make_shape_plot <- function(cluster_id, cluster_df, outdir, styles) {
  
  shape_features <- cluster_df |>
    rename(
      "Mean pLDDT" = scaled_mad_mean_pLDDT,
      "Eccentricity" = scaled_mad_eccentricity_u,
      "Coefficient of\nVariation" = scaled_mad_coeff_var_u,
      "Log(Kurtosis)" = scaled_mad_log_kurtosis_u,
      "Log(Length)" = scaled_mad_log_length_u
    ) |>
    select( # in order in which cols should appear in plot
      taxon,
      `Mean pLDDT`,
      `Log(Length)`,
      `Eccentricity`,
      `Coefficient of\nVariation`,
      `Log(Kurtosis)`
    )
  
  present_taxa <- unique(cluster_df$taxon)
  line_map <- styles$taxa_line[present_taxa]
  colour_map <- styles$taxa_colour[present_taxa]
  
  # Build full plot (standalone figure to embed in GMS cluster summaries)
  full_plt <- ggparcoord(shape_features,
                         columns = 2:6, groupColumn = "taxon",
                         scale = "uniminmax", showPoints = TRUE, alphaLines = 0.4,
                         mapping = aes(linetype = taxon), boxplot = TRUE,
  ) +
    scale_linetype_manual(values = line_map) +
    scale_colour_manual(values = colour_map) +
    labs(
      title = glue("Protein Shapes (GMS Cluster {cluster_id})"), 
      x = "", y = "Min-Max Scaled MAD-adjusted Value",
      linetype = "Taxon", color = "Taxon"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  # Write full plot to file
  png(
    glue("{outdir}/gms_{cluster_id}_protein_shapes_full.png"),
    width=12, height=8, units="in", res=300, type="cairo"
  )
  print(full_plt)
  invisible(dev.off())
  
  # Create base plot for use in grid plot
  base_plt <- ggparcoord(shape_features,
    columns = 2:6, groupColumn = "taxon", scale = "uniminmax", 
    showPoints = FALSE, alphaLines = 0, # Hide GGAlly points and lines
    boxplot = TRUE
  ) +
    geom_point(aes(color = taxon), size = 2.5) + # add custom large points
    geom_path( # add custom line style
      aes(color = taxon, linetype = taxon),
      linewidth = 1.2, alpha = 0.7
    ) +
    scale_linetype_manual(values = line_map) +
    scale_colour_manual(values = colour_map) +
    labs(title = glue("Cluster {cluster_id}"), x = "", y = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(size = 28, hjust = 0.5, margin = margin(b = 4))
    )
  
  # Write base plot to file
  png(
    glue("{outdir}/gms_{cluster_id}_protein_shapes.png"), 
    width=10, height=8, units="in", res=300, type="cairo"
  )
  print(base_plt)
  invisible(dev.off())
  
}




# Define taxon line style map and colour map
taxa_line_style_map <- c(
  Human = "dashed", Mouse = "longdash", Rat = "solid",
  Fly = "dotted", Worm = "dotdash", Zebrafish = "twodash"
)
taxa_colour_map <- c(
  Human = "#D55E00", Mouse = "#E69F00", Rat = "#F0E442",
  Fly = "#009E73", Worm = "#56B4E9", Zebrafish = "#CC79A7"
)

STYLES = list(
  taxa_line = taxa_line_style_map,
  taxa_colour = taxa_colour_map
)




### Create individual shape plots per GMS cluster -------------------------

log_info("Creating shape plots per GMS cluster ...")

df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
) |> mutate(taxon = stringr::str_to_title(taxon))
gms_clusters <- unique(sort(df$gms_cluster))

done <- future_lapply(gms_clusters, function (cluster_id, df, outdir, styles) {
  cluster_df <- df |> filter(gms_cluster == cluster_id)
  make_shape_plot(cluster_id, cluster_df, outdir, styles)
}, df = df, outdir = outdir, styles = STYLES)




# Create grid shape plot --------------------------------------------------

log_info("Creating grid shape plot ...")

imgs <- lapply(gms_clusters, function(cluster_id) {
  img <- readPNG(glue("{outdir}/gms_{cluster_id}_protein_shapes.png"))
  rasterGrob(img, interpolate = TRUE)
})
grid_plt <- grid.arrange(grobs = imgs, ncol = 5)
png(
  glue("{outdir}/shape_grid.png"),
  width = 20, height = 16, units = "in", res = 300
)
grid.draw(grid_plt)
dev.off()

