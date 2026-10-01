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

# Read in cluster dataframe and keep only GO IDs that have high odds ratio
# (gms_cluster, go_id, prop_in_cluster, prop_in_bg, p_value,
#  odds_ratio, p_adj, enriched, high_odds)
cluster_df_path = Path(
  f"results/5_describe_gms_cluster_func/cluster_dfs/{CLUSTER_ID}.csv"
)
cluster_df = pd.read_csv(str(cluster_df_path))
cluster_df = cluster_df[cluster_df["high_odds"]][["go_id", "odds_ratio"]]

# Read in GO ID metadata dataframe (go_id, lbl, namespace)
go_metadata_path = Path(f"results/1_filter_proteins/go_metadata.csv")
go_metadata_df = pd.read_csv(str(go_metadata_path))
go_metadata_df = go_metadata_df.rename(columns = {"lbl": "label"})

# Read in GO module members df (module, go_id)
go_module_df_path = Path(
  "results/3_extract_go_module_scores/go_module_members.csv"
)
go_module_df = pd.read_csv(str(go_module_df_path))

# Read in GO ID presence/absence per protein dataframe
# (uniprot, taxon, GO:0005515, ..., GO:0071881, module_1, ..., module_10)
go_presence_df_path = Path(
  "results/3_extract_go_module_scores/go_module_scores.csv"
)
go_presence_df = pd.read_csv(str(go_presence_df_path))

# Read in GMS cluster assignment dataframe and keep only proteins in this cluster
# (-,uniprot,--,gms_cluster,)
gms_cluster_assign_path = Path("results/4_cluster_by_gms/clustered_by_gms.csv")
gms_cluster_assign_df = pd.read_csv(str(gms_cluster_assign_path))
gms_cluster_proteins_df = gms_cluster_assign_df[
  gms_cluster_assign_df[gms_cluster] == CLUSTER_ID
][["uniprot", "gms_cluster"]]



### Join dataframes -------------------------------------------------------

# (go_id, module, label, namespace, odds_ratio)
merged = cluster_df
  .merge(go_module_df, on = "go_id", how = "inner")
  .merge(go_metadata_df, on = "go_id", how = "inner")

# GO P/A for GMS cluster proteins only
go_cols = [c for c in go_presence_df.columns if c.startswith("GO:")]
go_pa = gms_cluster_proteins_df
  .merge(go_presence_df, on = "uniprot", how = "inner")[go_cols]

# Compute GO ID frequencies (percentage of proteins with GO ID) in cluster
go_freq_df = go_pa.mean().mul(100).reset_index()
go_freq_df.columns = ["go_id", "percentage"]

# Merge GO ID percentage dataframe with merged GO metadata dataframe
cluster_go_df = go_freq_df.merge(merged, on = "go_id", how = "inner")

# Keep only those GO IDs present in at least 5% of cluster proteins
cluster_go_df = cluster_go_df[cluster_go_df["percentage"] >= 5]

cluster_go_df.to_csv(
  f"results/5_describe_gms_cluster_func/cluster_dfs/merged_{CLUSTER_ID}.csv"
)



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
A.graph_attr.update({
    "nodesep": "0.5",
    "ranksep": "0.5",
    "overlap": "false",
    "splines": "true"
})


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
A.layout("twopi")
A.draw(f"results/6_make_visual_summaries/gms_{CLUSTER_ID}_go_dag.svg")



