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
taxa <- c("human", "mouse", "rat", "fly", "worm", "zebrafish", "cross-taxa")
feature_spaces <- c("Shape", "GO Module Scores", "Shape\n+\nGO Module Scores")



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




### Define manifold visualisation functions -------------------------------

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




### Define UMAP embedding and clustering functions ------------------------

get_umap_embedding <- function (features) {
  umap(
    features, n_neighbors = 20, min_dist = 0.6,
    n_components = 2, metric = "euclidean",
    init = "pca", spread = 3
  )
}


do_umap <- function (df, taxon, feature_space) {
  # Use all proteins in df if taxon = "cross-taxa", otherwise subset dataframe
  if (taxon == "cross-taxa") { subset <- df } else {
    subset <- df |> filter(taxon == !!taxon) }
  
  features <- subset |> select(-uniprot, -taxon)
  if (feature_space == "Shape") {
    umap_embedding <- get_umap_embedding(
      features |> select(-starts_with("module_")))
  } else if (feature_space == "GO Module Scores") {
    umap_embedding <- get_umap_embedding(
      features |> select(starts_with("module_")))
  } else {
    umap_embedding <- get_umap_embedding(features)
  }
  
  subset$UMAP1 <- umap_embedding[,1]
  subset$UMAP2 <- umap_embedding[,2]
  return(list(plt = plot_manifold(subset), umapped_df = subset))
}




### Visualise manifolds within & across taxa on different feature spaces -- 

log_info("Creating all-manifolds grid plot ...")

# Build list of basic manifold plots and dataframes
plots <- list()
umapped_dfs <- list()
for (t in taxa) {
  for (s in feature_spaces) {
    umap_result <- do_umap(df, t, s)
    plots[[glue("{t}_{s}")]] <- umap_result[["plt"]]
    umapped_dfs[[glue("{t}_{s}")]] <- umap_result[["umapped_df"]]
  }
}

# Build ordered plots vector for 3x7 / vertical grid plot
ordered_keys <- unlist(
  lapply(taxa, function(t) paste(t, feature_spaces, sep = "_")))
ordered_plots <- plots[ordered_keys]

# Construct grid plot
grid_plt <- wrap_plots(ordered_plots, nrow = 7, ncol = 3)
col_labels <- lapply(feature_spaces, function(s) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = s, size = 6) +
    theme_void()
})
row_labels <- lapply(taxa, function(t) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = t, size = 5) +
    theme_void()
})
design <- c(
  area(1, 2, 1, 2), area(1, 3, 1, 3), area(1, 4, 1, 4), # Column labels
  
  # Row labels
  area(2, 1, 2, 1), area(3, 1, 3, 1), area(4, 1, 4, 1), area(5, 1, 5, 1),
  area(6, 1, 6, 1), area(7, 1, 7, 1), area(8, 1, 8, 1),
  
  area(2, 2, 8, 4) # Manifold plots grid
)
all_manifolds_plt <- wrap_plots(
  c(col_labels, row_labels, list(grid_plt)), design = design
) +
  plot_annotation(
    title = "UMAP Manifolds Across Taxa and Feature Spaces",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )

# Save manifolds plot
svg(glue("{outdir}/all_manifolds.svg"), width = 12, height = 20)
print(all_manifolds_plt)
invisible(dev.off())




### Cluster proteins in cross-taxa manifolds & visualise ------------------

# TODO



log_info("Computing shape-function clusters across all proteins ...")

# Create UMAP embedding
set.seed(254)
features <- df |> select(-uniprot, -taxon)
umap_embedding <- get_shapefunc_umap_embedding(features)
df$UMAP1 <- umap_embedding[,1]
df$UMAP2 <- umap_embedding[,2]

# Run HDBSCAN on UMAP embedding to discover initial UMAPped clusters
hdb <- hdbscan(umap_embedding, minPts = 5)
df$umap_cluster <- hdb$cluster # cluster 0 = noise
df$umap_cluster_prob <- hdb$membership_prob

# Visualise clusters on UMAP embedding, where colour represents cluster and 
# point shape represents taxon
cluster_levels <- levels(factor(df$umap_cluster))
cluster_colors <- c("0" = "grey",
                    setNames(
                      scales::hue_pal()(length(cluster_levels) - 1),
                      cluster_levels[cluster_levels != "0"]
                    )
)
umap_plot_path <- glue("{outdir}/manifold_shapefunc_cross_taxa_clustered.svg")
svg(umap_plot_path, width = 8, height = 6)
plt <- ggplot(df, aes(
  x = UMAP1, y = UMAP2, 
  color = factor(umap_cluster), 
  shape = factor(taxon))
) +
  geom_point(
    alpha = 0.7, size = 2,
    position = position_jitter(width = 0.1, height = 0.1)
  ) +
  scale_color_manual(name = "Cluster", values = cluster_colors) +
  scale_shape_discrete(name = "Taxon") +
  guides(
    color = guide_legend(ncol = 2, title = "Cluster"),
    shape = guide_legend(ncol = 1, title = "Taxon")
  ) +
  labs(
    title = stringr::str_wrap(
      "Shape-Function Clusters Across Proteins Present in Chemical Neuron-Neuron Synapses",
      width = 50
    ),
    subtitle = "UMAP of 16-D shape-function features; HDBSCAN clusters (minPts = 5)",
    x = "UMAP1", y = "UMAP2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 12)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.line = element_line(color = "black", linewidth = 0.4)
  )
print(plt)
invisible(dev.off())





log_info("Done!")












