# Visualise shape, function, and shape+function manifolds both within and across
# taxa, then run HDBSCAN on the manifold with the most structure.


outdir <- "results/4_cluster_by_gms"
log_file <- "logs/4_cluster_by_gms.log"
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
  "results/3_extract_go_module_scores/functional_features.csv",
  col_names = TRUE, show_col_types = FALSE) |>
  select(taxon, uniprot, starts_with("module_")
)
df <- inner_join(shape_df, go_module_df, by = c("uniprot", "taxon"))




### Define manifold visualisation functions -------------------------------

plot_manifold <- function (df, taxon, colour_by, colour_map) {
  
  # Define consistent shapes for cross-taxa plots
  taxon_shapes <- c(
    "human" = 16, # circle
    "mouse" = 15, # square
    "rat" = 17, # triangle
    "fly" = 18, # diamond
    "worm" = 8, # star
    "zebrafish" = 4 # cross
  )
  
  plt <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
  
  if (is.null(colour_by)) { # No colour mapping
    if (taxon == "cross-taxa") { # Draw taxa-coded point shapes
      plt <- plt +
        geom_point(
          aes(shape = factor(taxon)),
          alpha = 0.7, size = 2,
          position = position_jitter(width = 0.1, height = 0.1)
        ) +
        scale_shape_manual(values = taxon_shapes) +
        guides(shape = "none")
    } else { # Draw basic point shape
      plt <- plt +
        geom_point(
          alpha = 0.7, size = 2,
          position = position_jitter(width = 0.1, height = 0.1)
        )
    }
  } else { # Colour mapping
    if (taxon == "cross-taxa") { # Draw taxa-coded and custom coloured points
      plt <- plt +
        geom_point(
          aes(
            color = factor(.data[[colour_by]]),
            shape = factor(taxon)
          ),
          alpha = 0.7, size = 2,
          position = position_jitter(width = 0.1, height = 0.1)
        ) +
        scale_color_manual(name = colour_by, values = colour_map) +
        scale_shape_manual(values = taxon_shapes) +
        guides(color = "none", shape = "none")
    } else { # Draw basic point shapes with custom colour map
      plt <- plt +
        geom_point(
          aes(color = factor(.data[[colour_by]])),
          alpha = 0.7, size = 2,
          position = position_jitter(width = 0.1, height = 0.1)
        ) +
        scale_color_manual(name = colour_by, values = colour_map) +
        guides(color = "none")
    }
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

get_2x7_grid_layout <- function () {
  col_labels <- lapply(c("Shape", "GO Module Scores"), function(s) {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = s, size = 6) +
      theme_void()
  })
  row_labels <- lapply(taxa, function(t) {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = t, size = 5) +
      theme_void()
  })
  design <- c(
    # Column labels (row 1, columns 2 and 3)
    area(1, 2, 1, 2), area(1, 3, 1, 3),
    
    # Row labels (rows 2–8, column 1)
    area(2, 1, 2, 1), area(3, 1, 3, 1),
    area(4, 1, 4, 1), area(5, 1, 5, 1),
    area(6, 1, 6, 1), area(7, 1, 7, 1),
    area(8, 1, 8, 1),
    
    # Plot grid (rows 2–8, columns 2–3)
    # (area(top_row, left_col, bottom_row, right_col))
    area(2, 2, 8, 3)
  )
  return(
    list(col_labels = col_labels, row_labels = row_labels, design = design)
  )
}




### Define UMAP embedding and clustering functions ------------------------

get_umap_embedding <- function (features) {
  set.seed(254)
  umap(
    features, n_neighbors = 15, min_dist = 0.2,
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
    !names(subset) %in% c(
      "uniprot", "taxon", "gms_cluster", "gms_cluster_prob",
      "fine_ct_sf_cluster", "fine_ct_sf_cluster_prob", "coarse_ct_sf_cluster",
      "coarse_ct_sf_cluster_prob")])
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
    plt = plot_manifold(subset, taxon, colour_by, colour_map), 
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
df_to_cluster$gms_cluster <- hdb$cluster # cluster 0 = noise
df_to_cluster$gms_cluster_prob <- hdb$membership_prob
df <- inner_join(df, df_to_cluster) |>
  select(-UMAP1, -UMAP2) # drop umap coordinate columns

# Visualise cross-taxa GO module score cluster membership probabilities to
# quickly assess clustering confidence
svg(glue("{outdir}/cross_taxa_go_module_score_cluster_membership_probs.svg"),
    width = 12, height = 12)
plt <- ggplot(df, aes(gms_cluster_prob)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = glue("HDBSCAN Cross-Taxa GO Module Score Cluster Probabilities"))
print(plt)
invisible(dev.off())

# Write clusters to file
write_csv(df_to_cluster, glue("{outdir}/clustered_by_gms.csv"), col_names = TRUE)


# Create GO Module Score cluster membership colourmap
gms_clusters <- levels(factor(sort(df$gms_cluster)))
if (0 %in% gms_clusters) { # grey out noise points only if present
  colour_map <- c(
    "0" = "grey",
    setNames(
      scales::hue_pal()(length(gms_clusters) - 1),
      gms_clusters[gms_clusters != 0]
    )
  )
} else {
  colour_map <- setNames(
    scales::hue_pal()(length(gms_clusters)),
    gms_clusters
  )
}

# Recreate grid plot, colouring proteins by GO Module Score Cluster membership
plots <- list()
umapped_dfs <- list()
for (t in taxa) {
  for (s in c("Shape", "GO Module Scores")) {
    umap_result <- do_umap(
      df, t, s, colour_by = "gms_cluster", colour_map = colour_map
    )
    plots[[glue("{t}_{s}")]] <- umap_result[["plt"]]
    umapped_dfs[[glue("{t}_{s}")]] <- umap_result[["umapped_df"]]
  }
}
ordered_keys <- unlist(
  lapply(taxa, function(t) paste(t, c("Shape", "GO Module Scores"), sep = "_")))
ordered_plots <- plots[ordered_keys]
grid_plt <- wrap_plots(ordered_plots, nrow = 7, ncol = 2)
layout <- get_2x7_grid_layout()
go_cluster_coloured_manifolds_plt <- wrap_plots(
  c(layout[["col_labels"]], layout[["row_labels"]], list(grid_plt)), 
  design = layout[["design"]]
) +
  plot_annotation(
    title = "UMAP Manifolds Across Taxa and Feature Spaces",
    subtitle = "Proteins Coloured by Cross-Taxa GO Module Score Cluster Membership",
    theme = theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

svg(glue(
  "{outdir}/all_manifolds_coloured_by_ct_go_cluster_membership.svg"),
  width = 12, height = 20
)
print(go_cluster_coloured_manifolds_plt)
invisible(dev.off())





log_info("Done!")


