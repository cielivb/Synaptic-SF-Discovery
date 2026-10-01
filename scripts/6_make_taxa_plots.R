# Make plots describing relationships between taxa and GMS clusters, including:
# - Tile plot illustrating GMS cluster presence/absence per taxa
# - Stacked bar plot illustrating proportion of proteins in each taxa assigned 
#   to each GMS cluster, with taxa side-by-side, and bars stacked per taxa
# - One bar plot per GMS cluster illustrating proportion of proteins in each taxa
#   assigned to that cluster
# - One stacked single-column bar per GMS cluster illustrating proportion of 
#   proteins in that cluster belonging to each taxa (colour by taxa)
# - Grid plots for the above two plot types for ease of cluster comparison

outdir <- "results/6_make_taxa_plots"
log_file <- "logs/6_make_taxa_plots.log"
source("scripts/setup.R")
configure_logger(log_file)

library(ggplot2)

# Read in dataframe
df <- read_csv(
  "results/4_cluster_by_gms/clustered_by_gms.csv", show_col_types = FALSE
)




### Create tile plot illustrating GMS cluster presence/absence per taxa ---

log_info("Creating GMS Cluster P/A by Taxon Tile Plot ...")

# taxon | cluster | present (0/1)
pa_df <- df |>
  select(taxon, gms_cluster) |>
  mutate(present = 1) |>
  distinct() |>
  mutate(
    cluster = paste0("Cluster ", gms_cluster),
    present = factor(present),
    taxon = stringr::str_to_title(taxon)
  ) |>
  select(taxon, cluster, present)

# Sort clusters by cluster ID
pa_df$cluster <- factor(
  pa_df$cluster, levels = paste0("Cluster ", sort(unique(df$gms_cluster)))
)

# Sort taxa semi-phylogenetically
taxa_order <- c("Human", "Rat", "Mouse", "Zebrafish", "Fly", "Worm")
pa_df$taxon <- factor(pa_df$taxon, levels = taxa_order)


# Create tile plot
svg(glue("{outdir}/cluster_pa_by_taxa_tile_plot.svg"), height = 6, width = 20)
plt <- ggplot(pa_df, aes(x = cluster, y = taxon, fill = present)) +
  geom_tile(color = "grey80") +
  scale_fill_manual(values = c("0" = "white", "1" = "#2ECC71")) +
  labs(
    title = "GMS Cluster Presence By Taxa", x = "", y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b=12)),
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 55, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
  )
print(plt)
invisible(dev.off())




# TODO
# - Stacked bar plot illustrating proportion of proteins in each taxa assigned 
#   to each GMS cluster, with taxa side-by-side, and bars stacked per taxa
# - One bar plot per GMS cluster illustrating proportion of proteins in each taxa
#   assigned to that cluster
# - One stacked single-column bar per GMS cluster illustrating proportion of 
#   proteins in that cluster belonging to each taxa (colour by taxa)
# - Grid plots for the above two plot types for ease of cluster comparison




log_info("Done!")
