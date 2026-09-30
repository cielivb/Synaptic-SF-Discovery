# Visualise protein shapes, core GO IDs, and taxa proportions per GMS cluster,
# and visualise evolutionary graph of GMS cluster conservation

outdir <- "results/5_make_visual_summaries"
log_file <- "logs/5_make_visual_summaries.log"
source("scripts/setup.R")
configure_logger(log_file)

library(cowplot)
library(ggplot2)
library(GGally)
library(tidyr)
set.seed(254)


plan(multisession, workers = parallel::detectCores(logical = FALSE))
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
)




### Define global plot styling --------------------------------------------

# Define GMS cluster colour map
gms_clusters <- levels(factor(sort(df$gms_cluster)))
if (0 %in% gms_clusters) { # grey out noise points only if present
  gms_cluster_colour_map <- c( "0" = "grey", # noise GMS "cluster"
    setNames(
      scales::hue_pal()(length(gms_clusters) - 1),
      gms_clusters[gms_clusters != 0]
    )
  )
} else {
  gms_cluster_colour_map <- setNames(
    scales::hue_pal()(length(gms_clusters)),
    gms_clusters
  )
}

# Define GO module colour map
go_modules <- df |> select(starts_with("module_")) |> names()
go_module_colour_map <- setNames(
  viridis::viridis(length(go_modules)), go_modules
)

# Define taxon line style map and colour map
taxa_line_style_map <- c(
  human = "dashed", mouse = "longdash", rat = "solid",
  fly = "dotted", worm = "dotdash", zebrafish = "twodash"
)
taxa_colour_map <- c(
  human = "#D55E00", mouse = "#E69F00", rat = "#F0E442",
  fly = "#009E73", worm = "#56B4E9", zebrafish = "#CC79A7"
)

STYLES = list(
  gms_cluster = gms_cluster_colour_map,
  go_module = go_module_colour_map,
  taxa_line = taxa_line_style_map,
  taxa_colour = taxa_colour_map
)




### Define GMS cluster evolutionary conservation plot function (core result)

plot_gms_cluster_evo <- function () {
  NULL # TODO: Call the python script that actually creates the plot lol
}




### Define core visual summary functions ----------------------------------

# Create parallel coordinates plot showing shape features of proteins in a GMS 
# cluster, where lines are styled by taxon
make_shape_plot <- function(cluster_id, cluster_df, outdir, styles) {
  
  shape_features <- cluster_df |>
    select(taxon, starts_with("scaled_")) |>
    rename(
      "Mean pLDDT" = scaled_mad_mean_pLDDT,
      "Eccentricity" = scaled_mad_eccentricity_u,
      "Coefficient of Variation" = scaled_mad_coeff_var_u,
      "Log(Kurtosis)" = scaled_mad_log_kurtosis_u,
      "Log(Length)" = scaled_mad_log_length_u
    )
  present_taxa <- unique(cluster_df$taxon)
  line_map <- styles$taxa_line[present_taxa]
  colour_map <- styles$taxa_colour[present_taxa]
  
  # Build full plot
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
  svg(glue("{outdir}/gms_{cluster_id}_protein_shapes.svg"), width=12, height=8)
  print(full_plt)
  invisible(dev.off())
  
  # Create and return base plot (renamed title) for comparative patchwork grid
  base_plt <- full_plt +
    labs(title = glue("GMS Cluster {cluster_id}"))
  return(base_plt)
  
}


# Create GO graph showing most common GO IDs in GMS cluster, where GO IDs are
# coloured by their associated GO module
make_go_graph_plot <- function (cluster_id, cluster_df, outdir, styles) {
  NULL # TODO
}


# Create horizontal single-column stacked bar plot showing proportion of proteins 
# represented by each taxon present in GMS cluster (e.g., 10% of proteins in 
# GMS cluster are human proteins)
make_prop_repr_plot <- function(cluster_id, cluster_df, outdir, styles) {
  
  # Get proportion of proteins in cluster associated with each taxon
  prop_df <- cluster_df |> count(taxon) |> mutate(prop = n / sum(n))
  
  # Build full plot - aes(x=1) makes stacked
  full_plt <- ggplot(prop_df, aes(x = 1, y = prop, fill = taxon)) +
    geom_col(width = 0.6, color = "black") +
    coord_flip() + # make horizontal
    scale_fill_manual(values = styles$taxa_colour) +
    labs(
      title = glue("Taxon Representation (GMS Cluster {cluster_id})"),
      x = NULL, y = "Proportion", fill = "Taxon"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  
  # Write full plot to file
  svg(glue("{outdir}/gms_{cluster_id}_prop_repr.svg"), width = 12, height = 8)
  print(full_plt)
  invisible(dev.off())
  
  # Create and return base plot (renamed title) for comparative patchwork grid
  base_plt <- full_plt +
    labs(title = glue("GMS Cluster {cluster_id}"))
  return(base_plt)

}


# Create vertical bar plot showing percentage of proteins within each taxa
# that are assigned to this GMS cluster
make_perc_ded_plot <- function(cluster_id, cluster_df, df, outdir, styles) {
  
  # Compute percentage of proteins in each taxa present in this GMS cluster
  total_per_taxon <- df |> count(taxon) |> rename(total = n)
  cluster_total_per_taxon <- cluster_df |> 
    count(taxon) |> rename(cluster_total = n)
  perc_df <- left_join( # use left join as some taxa may not be in cluster
    total_per_taxon, cluster_total_per_taxon, by = "taxon") |>
    mutate(cluster_total = replace_na(cluster_total, 0)) |>
    mutate(perc = (cluster_total / total) * 100)
  
  # Build full plot
  full_plt <- ggplot(perc_df, aes(x = taxon, y = perc, fill = taxon)) +
    geom_text(
      aes(label = sprintf("%d%%", round(perc))), # print rounded perc
      size = 10, vjust = -0.5, # move % text above the bar
    ) +
    geom_col(width = 0.7, color = "black") +
    scale_fill_manual(values = styles$taxa_colour) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      title = glue(
        "Percentage of Taxon Proteins Assigned to GMS Cluster {cluster_id}"
      ),
      x = "Taxon", y = "Percentage of Taxon Proteins", fill = "Taxon"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Write full plot to file
  svg(glue("{outdir}/gms_{cluster_id}_perc_ded.svg"), width = 12, height = 8)
  print(full_plt)
  invisible(dev.off())
  
  # Create and return base plot (renamed title) for comparative patchwork grid
  base_plt <- full_plt +
    labs(title = glue("GMS Cluster {cluster_id}"))
  return(base_plt)
  
}




### Create "dictionaries" of cluster summary plots, keyed by GMS cluster ID,
# then by plot type
get_plots <- function(df) { 
  gms_clusters <- unique(sort(df$gms_cluster))
  plots_list <- future_lapply(
    gms_clusters,
    function(cluster_id, df, outdir, styles) {
      cluster_df <- df |> filter(gms_cluster == cluster_id)
      list(
        shape_plot = make_shape_plot(cluster_id, cluster_df, outdir, styles),
        go_graph = make_go_graph_plot(cluster_id, cluster_df, outdir, styles),
        prop_repr_plot = make_prop_repr_plot(
          cluster_id, cluster_df, outdir, styles),
        perc_ded_plot = make_perc_ded_plot(
          cluster_id, cluster_df, df, outdir, styles)
      )
    },
    df = df, outdir = outdir, styles = STYLES
  )
  names(plots_list) <- gms_clusters # Make list indexable by GMS cluster ID
  return(plots_list)
}




### Define patchwork layout functions --------------------------------------

# Comparative parallel coordinates (protein shapes) plot grid
make_shape_grid_plot <- function (plt_list) {
  shape_plts <- future_lapply(plt_list, function (cluster_plots) {
    cluster_plots$shape_plot
  })
  
  shape_grid_plt <- wrap_plots(
    shape_plts, ncol = 4, guides = "collect", axes = "keep", axis_titles = "keep"
  ) +
    plot_annotation(title = "Protein Shape Profiles Across GMS Clusters")
  
  
  pdf(glue("{outdir}/shape_grid_plot.pdf"), height = 12, width = 12)
  print(shape_grid_plt)
  invisible(dev.off())
}


# Comparative GO graphs grid
make_go_graph_grid_plot <- function (plt_list) {
  # TODO in python?
  NULL
}


# Comparative stacked bar (proportion of proteins representing each taxon) grid
make_prop_repr_grid_plot <- function (plt_list) {
  prop_repr_plts <- future_lapply(plt_list, function (cluster_plots) {
    cluster_plots$prop_repr_plot
  })
  prop_repr_grid_plt <- wrap_plots(prop_repr_plts) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Proportions of GMS Cluster Proteins Belonging To Each Taxa"
    )
  pdf(glue("{outdir}/prop_repr_grid_plot.pdf"), height = 12, width = 12)
  print(prop_repr_grid_plt)
  invisible(dev.off())
}


# Comparative "proportion per taxa dedicated" bar plots grid
make_perc_ded_grid_plot <- function (plt_list) {
  perc_ded_plts <- future_lapply(plt_list, function (cluster_plots) {
    cluster_plots$perc_ded_plot
  })
  perc_ded_grid_plt <- wrap_plots(perc_ded_plts) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Proportions of Taxon Proteins Assigned to GMS Clusters"
    )
  pdf(glue("{outdir}/perc_ded_grid_plot.pdf"), height = 12, width = 12)
  print(perc_ded_grid_plt)
  invisible(dev.off())
}


# GMS cluster summary grid
make_gms_cluster_summary_plots <- function (plt_list) {
  # TODO
  NULL
}




### MAIN ------------------------------------------------------------------

main <- function () {
  log_info("Making main plot ...")
  plot_gms_cluster_evo()
  
  log_info("Making grid plots ...")
  individual_plots <- get_plots(df) # keyed by GMS cluster ID, then by plot type
  make_shape_grid_plot(individual_plots)
  #make_go_graph_grid_plot(individual_plots)
  #make_prop_repr_grid_plot(individual_plots)
  #make_perc_ded_grid_plot(individual_plots)
  #make_gms_cluster_summary_plots(individual_plots)
  
  log_info("Done!")
}

main()






