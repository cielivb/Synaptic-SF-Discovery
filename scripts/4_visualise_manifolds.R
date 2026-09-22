# Visualise shape, function, and shape+function manifolds both within and across
# taxa, then run HDBSCAN on the manifold with the most structure.

# 7x3 grid of images.


outdir <- "results/4_visualise_manifolds"
log_file <- "logs/4_visualise_manifolds.log"
source("scripts/setup.R")
configure_logger(log_file)

library(dbscan)
library(ggplot2)
library(patchwork)
library(uwot)

set.seed(254)




### Read in dataframe of geometric and functional features ----------------

log_info("Loading data ...")

shape_df <- read_csv(
  "results/2_extract_shape_features/shape_features.csv",
  col_names = TRUE, show_col_types = FALSE
)
go_module_df <- read_csv(
  "results/3_extract_functional_features/functional_features.csv",
  col_names = TRUE, show_col_types = FALSE) |>
  select(taxon, uniprot, starts_with("module_")
)
df <- inner_join(shape_df, go_module_df, by = c("uniprot", "taxon"))




### Define each feature-space's UMAP embedding parameters -----------------

get_shape_umap_embedding <- function (features) {
  umap(
    features, n_neighbors = 15, min_dist = 0.4,
    n_components = 2, metric = "euclidean",
    init = "pca", spread = 1
  )
}

get_func_umap_embedding <- function (features) {
  umap(
    features, n_neighbors = 15, min_dist = 0.4,
    n_components = 2, metric = "euclidean",
    init = "pca", spread = 1
  )
}

get_shapefunc_umap_embedding <- function (features) {
  umap(
    features, n_neighbors = 15, min_dist = 0.4,
    n_components = 2, metric = "euclidean",
    init = "pca", spread = 1
  )
}




### Define manifold visualisation function --------------------------------

plot_manifold <- function (df) {
  plt <- ggplot(df, aes(x = df$UMAP1, y = df$UMAP2)) +
    geom_point(
      alpha = 0.7, size = 2,
      position = position_jitter(width = 0.1, height = 0.1)
    ) + 
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid.major = element_blank(), # Remove major gridlines
      panel.grid.minor = element_blank(), # Remove minor gridlines
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
  return(plt)
}


make_umap_plot <- function (df, taxon, feature_space) {
  log_info(glue("Making {taxon} {feature_space} UMAP manifold plot ..."))
  
  # Use all proteins in df if taxon = "cross-taxa", otherwise subset dataframe
  if (taxon == "cross-taxa") { subset <- df } else {
    subset <- df |> filter(taxon == !!taxon) }
  
  features <- subset |> select(-uniprot, -taxon)
  if (feature_space == "shape") {
    umap_embedding <- get_shape_umap_embedding(features)
  } else if (feature_space == "func") {
    umap_embedding <- get_func_umap_embedding(features)
  } else {
    umap_embedding <- get_shapefunc_umap_embedding(features)
  }
  
  subset$UMAP1 <- umap_embedding[,1]
  subset$UMAP2 <- umap_embedding[,2]
  return(plot_manifold(subset))
}


taxa <- c("human", "mouse", "rat", "fly", "worm", "zebrafish", "cross-taxa")
feature_spaces <- c("shape", "func", "shapefunc")
plots <- list()
for (t in taxa) {
  for (s in feature_spaces) {
    plots[[glue("{t}_{s}")]] <- make_umap_plot(df, t, s)
  }
}


# Build ordered vector for 3x7 / vertical plot
ordered_keys <- unlist(
  lapply(taxa, function(t) paste(t, feature_spaces, sep = "_")))
ordered_plots <- plots[ordered_keys]

# Build vertical plot
plot_cols <- c("", "shape", "func", "shapefunc")
grid_plt <- wrap_plots(ordered_plots, nrow = 7, ncol = 3, guides = "collect")
col_labels <- wrap_plots(
  lapply(plot_cols, function(s) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = s, size = 6) +
      theme_void()
  }),
  ncol = 4, nrow = 1
)
row_labels <- wrap_plots(
  lapply(taxa, function(t) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = t, size = 5) +
      theme_void()
  }),
  nrow = 7, ncol = 1
)
final_grid <- row_labels | grid_plt
final_grid <- col_labels / final_grid
final_grid <- final_grid +
  plot_annotation(
    title = "UMAP Manifolds Across Taxa and Feature Spaces",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )



svg(glue("{outdir}/manifold_grid.svg"), width = 12, height = 20)
print(final_grid)
invisible(dev.off())



log_info("Done!")












