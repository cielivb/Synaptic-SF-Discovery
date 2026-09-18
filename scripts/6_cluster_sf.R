# Are there any distinct clusters when you combine both shape and function into
# a single feature space across taxa instead of treating them separately?


source("scripts/setup.R")
log_file <- "logs/6_cluster_sf.log"
configure_logger(log_file)

library(dbscan)
library(ggplot2)
library(uwot)

dir.create("results/6_cluster_sf", showWarnings = FALSE)


# Read in dataframe of geometric and functional features
df <- read_csv(
  "results/5_get_functional_modules/functional_features.csv",
  col_names = TRUE, show_col_types = FALSE
)

log_info("Computing shape-function clusters across all proteins ...")

# Create UMAP embedding
set.seed(254)
features <- df |> select(-uniprot, -taxon)
umap_embedding <- umap(
  features,
  n_neighbors = 15,
  min_dist = 0.4,
  n_components = 2,
  metric = "euclidean",
  init = "pca",
  spread = 1
)
df$UMAP1 <- umap_embedding[,1]
df$UMAP2 <- umap_embedding[,2]

# Run HDBSCAN on UMAP embedding to discover UMAPped clusters
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
umap_plot_path <- glue("results/6_cluster_sf/umap_cluster.svg")
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
