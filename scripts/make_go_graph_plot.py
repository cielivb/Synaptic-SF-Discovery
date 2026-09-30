""" Create a GO ID graph depicting enriched GO IDs a given GMS cluster """

import argparse
import json
import networkx as nx
from networkx.drawing.nx_agraph import to_agraph
import pandas as pd
from pathlib import Path
import textwrap


# Read in cluster ID
parser = argparse.ArgumentParser()
parser.add_argument("cluster_id", type = int)
args = parser.parse_args()
CLUSTER_ID = args.cluster_id




### Read in and wrangle initial/pre-join dataframes -----------------------

# Read in cluster dataframe
cluster_df_path = Path(
  f"results/6_make_visual_summaries/cluster_dfs/{CLUSTER_ID}.csv"
)
cluster_df = pd.read_csv(str(cluster_df_path))

# Read in GO ID metadata dataframe (maps GO IDs to GO names/descriptions)
go_metadata_path = Path(f"results/1_filter_proteins/go_metadata.csv")
go_metadata_df = pd.read_csv(str(go_metadata_path))

# Read in GO module members df (maps GO IDs to their GO modules/GO ID clusters)
go_module_df_path = Path(
  "results/3_extract_functional_features/go_module_members.csv"
)
go_module_df = pd.read_csv(str(go_module_df_path))

# Read in GO ID presence/absence per protein dataframe
go_presence_df_path = Path(
  "results/3_extract_functional_features/functional_features.csv"
)
go_presence_df = pd.read_csv(str(go_presence_df_path))

# Read in GO IDs in cluster with effect size (odds ratio) >= 2.5
enriched_go_ids_path = Path(
  "results/5_find_enriched_go_ids/go_enrichment_results.csv"
)
enriched_go_ids_df = pd.read_csv(str(enriched_go_ids_path))
enriched_go_ids_df = enriched_go_ids_df[
    (enriched_go_ids_df["gms_cluster"] == CLUSTER_ID) &
    (enriched_go_ids_df["odds_ratio"] >= 2.5)
][["go_id", "odds_ratio"]]




### Join dataframes -------------------------------------------------------

# Get dataframe containing GO ID, label, namespace, and GO module for enriched
# GO IDs only
go_metadata_df = (
    go_metadata_df.merge(enriched_go_ids_df, on = "go_id", how = "inner")
    .merge(go_module_df, on = "go_id", how = "inner")
    .rename(columns={"lbl": "label"})
)

# Get dataframe containing only rows present in cluster dataframe, and only GO
# ID columns where GO ID is present in at least one of these rows
cluster_df = cluster_df.merge(go_presence_df, on = "uniprot", how = "inner")
go_cols = [c for c in cluster_df.columns if c.startswith("GO:")]
keep_cols = [c for c in go_cols if cluster_df[c].any()] # >=1 non-zero val
cluster_go_df = cluster_df.drop(
  columns=[c for c in cluster_df.columns if c not in keep_cols]
)
cluster_go_df = cluster_go_df.astype(int)

# Compute GO ID frequencies (percentage of proteins with GO ID) in cluster
go_freq_df = cluster_go_df.mean().mul(100).reset_index()
go_freq_df.columns = ["go_id", "percentage"]

# Keep only those GO IDs present in at least 5% of cluster proteins
go_freq_df = go_freq_df[go_freq_df["percentage"] >= 5]

# Merge GO ID percentage dataframe with GO metadata dataframe
cluster_go_df = go_freq_df.merge(go_metadata_df, on = "go_id", how = "inner")
cluster_go_df.to_csv(f"results/6_make_visual_summaries/cluster_dfs/merged{CLUSTER_ID}.csv")



### Create ancestor -> descendant DAG -------------------------------------

# Read in JSON mapping GO IDs to lists of ancestors
json_path = Path("results/1_filter_proteins/synaptic_go_ancestors.json")
with open(json_path) as file:
    ancestor_map = json.load(file)

G = nx.DiGraph()

# Add nodes plus node metadata
for _, row in cluster_go_df.iterrows():
    G.add_node(
        row["go_id"], 
        label = row["label"], namespace = row["namespace"],
        module = row["module"], percentage = row["percentage"],
        odds_ratio = row["odds_ratio"]
    )

# Add DAG edges connecting nodes to closest ancestors
for child, ancestors in ancestor_map.items():
    if child not in G.nodes: continue
    valid_ancestors = set([a for a in ancestors if a in G.nodes]) 
    keep_ancestors = valid_ancestors.copy()
  
    # Iteratively remove ancestors as closer ancestors are discovered
    # That is, if c (child) -> b, and b -> a, drop c -> a
    for b in valid_ancestors:
        if not ancestor_map.get(b, None): continue
        # We know c -> b by definition of b being an ancestor of c
        # Now remove all c -> a where b -> a exists
        for a in valid_ancestors:
            if a in ancestor_map[b]:
                if a in keep_ancestors:
                    keep_ancestors.remove(a)
          
    # Add DAG edges between remaining ancestors in keep_ancestors and child
    for ancestor in keep_ancestors:
        G.add_edge(ancestor, child)



  
### Make a cool plot! -----------------------------------------------------

A = to_agraph(G)

go_module_cmap_path = Path(
  "results/6_make_visual_summaries/go_module_colours.csv"
)
go_module_cmap_df = pd.read_csv(str(go_module_cmap_path))
go_module_cmap_df["module"] = (
    go_module_cmap_df["module"]
    .str.replace("module_", "", regex=False)
    .astype(int)
)

module_to_colour = dict(zip(go_module_cmap_df.module, go_module_cmap_df.colour))

namespace_to_linestyle = {
    "biological_process": "solid",
    "molecular_function": "dashed",
    "cellular_component": "dotted"
}


# Set node attributes
for meta_node in G.nodes:
    node = A.get_node(meta_node)
    data = G.nodes[meta_node]

    # Node box outline properties
    node.attr["shape"] = "box"
    node.attr["color"] = module_to_colour[data["module"]]  # outline colour
    node.attr["style"] = namespace_to_linestyle[data["namespace"]]
    node.attr["style"] += ",rounded"
    node.attr["penwidth"] = "3"  # box outline width
    node.attr["margin"] = "0.1,0.05"
    node.attr["fixedsize"] = "false"

    # Text properties
    node.attr["fontsize"] = str(8 * 1.1 ** data["percentage"])
    #node.attr["fontsize"] = str(8 + data["percentage"])
    node.attr["label"] = "\n".join(textwrap.wrap(data["label"], width=25))

# Render plot
A.layout("fdp")
A.draw(f"results/6_make_visual_summaries/gms_{CLUSTER_ID}_go_dag.svg")



