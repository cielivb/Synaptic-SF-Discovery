# Visualise shape, function, and shape+function manifolds both within and across
# taxa, then run HDBSCAN on the manifold with the most structure.


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

plot_manifold <- function (df, colour_by, colour_map) {
  plt <- ggplot(df, aes(x = UMAP1, y = UMAP2)) + 
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid.major = element_blank(), # Remove major gridlines
      panel.grid.minor = element_blank(), # Remove minor gridlines
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
  if (is.null(colour_by)) {
    plt <- plt + geom_point(
      alpha = 0.7, size = 2,
      position = position_jitter(width = 0.1, height = 0.1)
    )
  } else {
    plt <- plt + geom_point(
        aes(color = factor(.data[[colour_by]])),
        alpha = 0.7, size = 2,
        position = position_jitter(width = 0.1, height = 0.1)
      ) +
      scale_color_manual(name = colour_by, values = colour_map) +
      guides(color = "none") # hide legend
  }
  return(plt)
}

get_3x7_grid_layout <- function () {
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
  return(
    list(col_labels = col_labels, row_labels = row_labels, design = design)
  )
}




### Define UMAP embedding and clustering functions ------------------------

get_umap_embedding <- function (features) {
  set.seed(254)
  umap(
    features, n_neighbors = 20, min_dist = 0.6,
    n_components = 2, metric = "euclidean",
    init = "pca", spread = 3
  )
}


do_umap <- function (df, taxon, feature_space, colour_by=NULL, colour_map=NULL) {
  # Use all proteins in df if taxon = "cross-taxa", otherwise subset dataframe
  if (taxon == "cross-taxa") { subset <- df } else {
    subset <- df |> filter(taxon == !!taxon) }
  
  # Keep row and column order consistent across runs, and drop non-feature cols
  feature_cols <- sort(names(subset)[
    !names(subset) %in% c("uniprot", "taxon", "go_module_cluster", 
                          "go_module_cluster_prob")])
  features <- subset |> arrange(uniprot) |> select(all_of(feature_cols))

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
  return(list(
    plt = plot_manifold(subset, colour_by, colour_map), 
    umapped_df = subset
  ))
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

# Build ordered plots vector for vertical grid plot and construct grid plot
ordered_keys <- unlist(
  lapply(taxa, function(t) paste(t, feature_spaces, sep = "_")))
ordered_plots <- plots[ordered_keys]
grid_plt <- wrap_plots(ordered_plots, nrow = 7, ncol = 3)
layout <- get_3x7_grid_layout()
all_manifolds_plt <- wrap_plots(
  c(layout[["col_labels"]], layout[["row_labels"]], list(grid_plt)), 
  design = layout[["design"]]
) +
  plot_annotation(
    title = "UMAP Manifolds Across Taxa and Feature Spaces",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )

svg(glue("{outdir}/all_manifolds.svg"), width = 12, height = 20)
print(all_manifolds_plt)
invisible(dev.off())




### Cluster cross-taxa GO module scores then recreate grid plot, colouring
# proteins by cross-taxa GO module score cluster membership. --------------

# This will provide evolutionary functional insight and help determine whether 
# it is worth investigating the large cross-taxa shape-function cluster's
# substructure.

log_info("Clustering proteins by cross-taxa GO module scores ...")

# Assign cross-taxa GO module score clusters to each protein
df_to_cluster <- umapped_dfs[["cross-taxa_GO Module Scores"]]
hdb <- hdbscan(df_to_cluster |> select(UMAP1, UMAP2), minPts = 10)
df_to_cluster$go_module_cluster <- hdb$cluster # cluster 0 = noise
df_to_cluster$go_module_cluster_prob <- hdb$membership_prob
df <- inner_join(df, df_to_cluster) |>
  select(-UMAP1, -UMAP2) # drop umap coordinate columns

# Visualise cross-taxa GO module score cluster membership probabilities to
# quickly assess clustering confidence
svg(glue("{outdir}/cross_taxa_go_module_score_cluster_membership_probs.svg"),
    width = 12, height = 12)
plt <- ggplot(df, aes(go_module_cluster_prob)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = glue("HDBSCAN Cross-Taxa GO Module Score Cluster Probabilities"))
print(plt)
invisible(dev.off())


# Create GO Module Score cluster membership colourmap
go_module_clusters <- levels(factor(df$go_module_cluster))
colour_map <- c("0" = "grey", # 0 = noise
  setNames(
    scales::hue_pal()(length(go_module_clusters) - 1),
    go_module_clusters[go_module_clusters != "0"]
  )
)

# Recreate grid plot, colouring proteins by GO Module Score Cluster membership
plots <- list()
umapped_dfs <- list()
for (t in taxa) {
  for (s in feature_spaces) {
    umap_result <- do_umap(
      df, t, s, colour_by = "go_module_cluster", colour_map = colour_map
    )
    plots[[glue("{t}_{s}")]] <- umap_result[["plt"]]
    umapped_dfs[[glue("{t}_{s}")]] <- umap_result[["umapped_df"]]
  }
}
ordered_keys <- unlist(
  lapply(taxa, function(t) paste(t, feature_spaces, sep = "_")))
ordered_plots <- plots[ordered_keys]
grid_plt <- wrap_plots(ordered_plots, nrow = 7, ncol = 3)
layout <- get_3x7_grid_layout()
go_cluster_coloured_manifolds_plt <- wrap_plots(
  c(layout[["col_labels"]], layout[["row_labels"]], list(grid_plt)), 
  design = layout[["design"]]
) +
  plot_annotation(
    title = "UMAP Manifolds Across Taxa and Feature Spaces",
    subtitle = "Proteins Coloured by Cross-Taxa GO Module Score Cluster Membership",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )

svg(glue(
  "{outdir}/all_manifolds_coloured_by_ct_go_cluster_membership.svg"),
  width = 12, height = 20
)
print(go_cluster_coloured_manifolds_plt)
invisible(dev.off())




### Inspect UMAP embedding for large cross-taxa shape+func cluster alone --

# The coloured grid plot shows there is functional substructure within this
# large cluster, but it is a bit difficult to see clearly, likely due to the
# other more extreme clusters in shape+function space. It is worth inspecting
# this large cluster and plotting it side-by-side with its original manifold.

log_info("Inspecting large cross-taxa shape-func cluster ...")

# Use HDBSCAN to discover "coarse" / top-level clusters
df_to_cluster <- do_umap(
  df, "cross-taxa", "Shape\n+\nGO Module Scores"
)[["umapped_df"]]
hdb <- hdbscan(df_to_cluster |> select(UMAP1, UMAP2), minPts = 10)
df_to_cluster$coarse_ct_sf_cluster <- hdb$cluster # cluster 0 = noise
df_to_cluster$coarse_ct_sf_cluster_prob <- hdb$membership_prob

# Visualise cluster membership probabilities to assess coarse clustering confidence
svg(glue("{outdir}/ct_sf_coarse_cluster_probs.svg"), width = 12, height = 12)
plt <- ggplot(df_to_cluster, aes(coarse_ct_sf_cluster_prob)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = glue(
    "HDBSCAN Cross-Taxa Shape + GO Module Score Coarse Cluster Probabilities"
  ))
print(plt)
invisible(dev.off())

# Isolate large coarse cluster - it has the most proteins
large_coarse_cluster_protein_df <- df_to_cluster |>
  add_count(coarse_ct_sf_cluster, name = "coarse_cluster_size") |>
  filter(coarse_cluster_size == max(coarse_cluster_size)) |>
  select(-coarse_cluster_size)
log_info("Protein counts per coarse cluster in cross-taxa shape-func manifold:")
print(large_coarse_cluster_protein_df |> count(coarse_ct_sf_cluster))

# Run UMAP on large coarse cluster alone, including GO Module cluster colourmap
subdf_umap_result <- do_umap(
  large_coarse_cluster_protein_df, "cross-taxa", "Shape\n+\nGO Module Scores",
  colour_by = "go_module_cluster", colour_map = colour_map
)
large_cluster_embedding_plt <- subdf_umap_result[["plt"]]
svg(glue("{outdir}/ct_sf_large_cluster_umap.svg"), width = 8, height = 8)
print(large_cluster_embedding_plt)
invisible(dev.off())
str(umap_result[["umapped_df"]])

# Run HDBSCAN on large cluster UMAP embedding to isolate "fine" clusters
subdf_to_cluster <- subdf_umap_result[["umapped_df"]]
hdb <- hdbscan(subdf_to_cluster |> select(UMAP1, UMAP2), minPts = 10)
subdf_to_cluster$fine_ct_sf_cluster <- hdb$cluster # cluster 0 = noise
subdf_to_cluster$fine_ct_sf_cluster_prob <- hdb$membership_prob

log_info("Protein counts per fine cluster in large coarse cluster:")
print(subdf_to_cluster |> count(fine_ct_sf_cluster))

# Visualise fine clusters on large cluster UMAP embedding
log_info("making colour map ...")
fine_clusters <- levels(factor(subdf_to_cluster$fine_ct_sf_cluster))
cm_fine_clusters <- c(
  "0" = "grey", # 0 = noise
  setNames(
    viridis::viridis(length(fine_clusters) - 1, option = "D"),
    fine_clusters[fine_clusters != "0"]
  )
)
log_info("made colour map")
svg(glue("{outdir}/ct_sf_large_cluster_fine_clusters.svg"), width = 8, height=8)
print(do_umap(
  subdf_to_cluster, "cross-taxa", "Shape\n+\nGO Module Scores",
  colour_by = "fine_ct_sf_cluster", colour_map = cm_fine_clusters
))
invisible(dev.off())

# Visualise cluster membership probabilities to assess fine clustering confidence
svg(glue("{outdir}/ct_sf_fine_cluster_probs.svg"), width = 12, height = 12)
plt <- ggplot(subdf_to_cluster, aes(fine_ct_sf_cluster_prob)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = glue(
    "HDBSCAN Cross-Taxa Shape + GO Module Score Fine Cluster Probabilities"
  ))
print(plt)
invisible(dev.off())


# Display large cluster embedding plot next to full cross-taxa shape-func plot
# with an arrow pointing to represented cluster


log_info("Done!")











