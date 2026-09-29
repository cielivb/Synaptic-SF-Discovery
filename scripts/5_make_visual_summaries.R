# Visualise protein shapes, core GO IDs, and taxa proportions per GMS cluster,
# and visualise evolutionary graph of GMS cluster conservation

outdir <- "results/5_make_visual_summaries"
log_file <- "logs/5_make_visual_summaries.log"
source("scripts/setup.R")
configure_logger(log_file)

library(ggplot2)
library(GGally)




### Define GMS cluster evolutionary conservation plot function (core result)

plot_gms_cluster_evo <- function () {
  NULL # TODO: Call the python script that actually creates the plot lol
}



### Define core visual summary functions ----------------------------------

# Create parallel coordinates plot showing shape features of proteins in a GMS 
# cluster, where lines are styled by taxon


# Create GO graph showing most common GO IDs in GMS cluster, where GO IDs are
# coloured by their associated GO module


# Create horizontal single-column stacked bar plot showing proportion of proteins 
# represented by each taxon present in GMS cluster (e.g., 10% of proteins in 
# GMS cluster are human proteins)


# Create vertical bar plot showing proportion of proteins within each taxa
# that are assigned to this GMS cluster




### Create "dictionaries" of cluster summary plots, keyed by GMS cluster ID
get_shape_plots <- function (df) { ... }
get_go_graphs <- function (df) { ... }
get_prop_repr_graphs <- function (df) { ... }
get_prop_ded_graphs <- function (df) { ... }



### Define patchwork layout functions --------------------------------------

# Comparative parallel coordinates (protein shapes) plot grid

# Comparative GO graphs grid

# Comparative stacked bar (proportion of proteins representing each taxon) grid

# Comparative "proportion per taxa dedicated" bar plots grid

# GMS cluster summary grid




### Create patchwork plots ------------------------------------------------

make_plots <- function (df) {
  shape_plots <- get_shape_plots(df)
  go_graphs <- get_go_graphs(df)
  prop_repr_graphs <- get_prop_repr_graphs(df)
  prop_ded_graphs <- get_prop_ded_graphs(df)
  
  # TODO: get patchwork layouts then stitch em up
}


### MAIN ------------------------------------------------------------------

main <- function () {
  plot_gms_cluster_evo()
  df <- read_csv(
    "results/4_cluster_by_gms/clustered_by_gms.csv",
    show_col_types = FALSE
  )
  plan(multisession, workers = parallel::detectCores(logical = FALSE))
  make_plots(df)
}

main()






