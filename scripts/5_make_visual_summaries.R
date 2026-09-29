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

# Define taxon line style map
taxa <- levels(factor(sorted(df$taxon)))
taxa_line_style_map <- setNames(
  c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), taxa
)

STYLES = list(
  gms_cluster = gms_cluster_colour_map,
  go_module = go_module_colour_map,
  taxa = taxa_line_style_map
)



### Define GMS cluster evolutionary conservation plot function (core result)

plot_gms_cluster_evo <- function () {
  NULL # TODO: Call the python script that actually creates the plot lol
}




### Define core visual summary functions ----------------------------------

# Create parallel coordinates plot showing shape features of proteins in a GMS 
# cluster, where lines are styled by taxon
make_shape_plot <- function (cluster_id, cluster_df, outdir, styles) {
  # Create base plot (no titles, legends etc)
  ... # TODO
  # Write full plot (standalone plot with titles etc)
  # Return base plot
}

# Create GO graph showing most common GO IDs in GMS cluster, where GO IDs are
# coloured by their associated GO module
make_go_graph_plot <- function (cluster_id, cluster_df, outdir, styles) {
  ... # TODO
}

# Create horizontal single-column stacked bar plot showing proportion of proteins 
# represented by each taxon present in GMS cluster (e.g., 10% of proteins in 
# GMS cluster are human proteins)
make_prop_repr_plot <- function (cluster_id, cluster_df, outdir, styles) {
  ... # TODO
}

# Create vertical bar plot showing proportion of proteins within each taxa
# that are assigned to this GMS cluster
make_prop_ded_plot <- function (cluster_id, cluster_df, outdir, styles) {
  ... # TODO
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






