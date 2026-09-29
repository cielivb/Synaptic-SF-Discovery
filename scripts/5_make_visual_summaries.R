# Visualise protein shapes, core GO IDs, and taxa proportions per GMS cluster,
# and visualise evolutionary graph of GMS cluster conservation

outdir <- "results/5_make_visual_summaries"
log_file <- "logs/5_make_visual_summaries.log"
source("scripts/setup.R")
configure_logger(log_file)

library(ggplot2)
library(GGally)


plan(multisession, workers = parallel::detectCores(logical = FALSE))
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv",
  show_col_types = FALSE
)




### Define global plot styling --------------------------------------------

# Define GMS cluster colour map
gms_clusters <- levels(factor(sorted(df$gms_cluster)))
if (0 %in% gms_clusters) { # grey out noise points only if present
  gms_cluster_colour_map <- c( "0" = "grey", # noise GMS "cluster"
    setNames(
      scales::hue_pal()(length(go_module_clusters) - 1),
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
taxa <- levels(factor(sorted(df$taxon)))
taxa_line_style_map <- setNames(
  c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), taxa
)
# TODO: taxa colour map

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
  
  # Build full plot
  full_plt <- ggparcoord(cluster_df,
    columns = starts_with("scaled_"), groupColumn = "taxon",
    order = "anyClass", showPoints = TRUE, alphaLines = 0.3
  ) +
    scale_linetype_manual(values = styles$taxa_line) +
    scale_colour_manual(values = styles$taxa_colour)
    labs(
      title = glue("Protein Shapes (GMS Cluster {cluster_id})"), 
      linetype = "Taxon", color = "Taxon"
    ) +
    theme_ipsum() +
    theme(
      plot.title = element_text(size = 10)
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
  ... # TODO
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
    scale_fill_manual(values = styles$taxa_colors) +
    labs(
      title = glue("Taxon Representation (GMS Cluster {cluster_id})"),
      x = NULL, y = "Proportion", fill = "Taxon"
    ) +
    theme_ipsum() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 12)
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



# Create vertical bar plot showing proportion of proteins within each taxa
# that are assigned to this GMS cluster
make_prop_ded_plot <- function (cluster_id, cluster_df, outdir, styles) {
  
  # Build full plot
  full_plt <- ... # TODO
  
  # Write full plot to file
  svg(glue("{outdir}/gms_{cluster_id}_prop_ded.svg"), width=12, height=8)
  print(full_plt)
  invisible(dev.off())
  
  # Create and return base plot (renamed title) for comparative patchwork grid
  base_plt <- ... # TODO
  return(base_plt)
  
}




### Create "dictionaries" of cluster summary plots, keyed by GMS cluster ID,
# then by plot type
get_plots <- function(df) { 
  gms_clusters <- unique(sorted(df$gms_cluster))
  plots_list <- future_lapply(
    gms_clusters,
    function(cluster_id, df, outdir, styles) {
      cluster_df <- df |> filter(gms_cluster == cluster_id)
      list(
        shape_plot = make_shape_plot(cluster_id, cluster_df, outdir, styles),
        go_graph = make_go_graph_plot(cluster_id, cluster_df, outdir, styles),
        prop_repr_plot = make_prop_repr_plot(cluster_id, cluster_df, outdir, styles),
        prop_ded_plot = make_prop_ded_plot(cluster_id, cluster_df, outdir, styles)
      )
    },
    df = df, outdir = outdir, styles = STYLES
  )
  names(plots_list) <- gms_clusters # Make list indexable by GMS cluster ID
  return(plots_list)
}




### Define patchwork layout functions --------------------------------------

# Comparative parallel coordinates (protein shapes) plot grid

# Comparative GO graphs grid

# Comparative stacked bar (proportion of proteins representing each taxon) grid

# Comparative "proportion per taxa dedicated" bar plots grid

# GMS cluster summary grid




### Create patchwork plots ------------------------------------------------

make_plots <- function (df) {
  individual_plots <- get_plots(df) # keyed by GMS cluster ID, then by plot type
  
  # TODO: get patchwork layouts then stitch em up
}




### MAIN ------------------------------------------------------------------

main <- function () {
  plot_gms_cluster_evo()
  make_plots(df)
}

main()






