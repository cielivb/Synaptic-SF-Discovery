from pathlib import Path
import networkx as nx
import json
import pandas as pd


def extract_go(uri):
    """ Extract GO ID out of GO URI """
    return uri.split("/")[-1].replace("GO_", "GO:")


# Load synaptic GO IDs
synaptic_go_path = Path("results/1_filter_by_func/goslim_synapse_go_ids.csv")
synaptic_go_ids = pd.read_csv(str(synaptic_go_path))["go_id"].tolist()

# Load GO-basic JSON into directed graph
go_graph = nx.DiGraph()
with open("data/go-basic.json") as file:
    go = json.load(file)
for edge in go["graphs"][0]["edges"]: # Add child -> parent edges
    child, parent = edge["sub"], edge["obj"]
    go_graph.add_edge(extract_go(child), extract_go(parent))
for node in go["graphs"][0]["nodes"]: # Add node metadata
    node_id = extract_go(node["id"])
    go_graph.add_node(node_id, **node)


# Map synaptic GO IDs to their ancestor IDs
ancestor_map = {}
for go_id in synaptic_go_ids:
    if go_id in go_graph:
        ancestors = sorted(list(nx.ancestors(go_graph, go_id)))
        ancestor_map[go_id] = ancestors
    else:
        ancestor_map[go_id] = []

# Write JSON mapping
outpath = Path("results/1_filter_by_func/synaptic_go_ancestors.json")
with open(outpath, "w") as outfile:
    json.dump(ancestor_map, outfile, indent=2)


# Write synaptic GO ID + ancestor metadata csv
target_go_ids, rows = set(), []
for synaptic_go_id, ancestors in ancestor_map.items():
  target_go_ids.add(synaptic_go_id)
  target_go_ids.update(ancestors)

for go_id in target_go_ids:
    if go_id in go_graph.nodes:
        attrs = go_graph.nodes[go_id]
        row = {"go_id": go_id}
        row.update(attrs)
        rows.append(row)
metadata_df = pd.DataFrame(rows)
metadata_df.to_csv("results/1_filter_by_func/go_metadata.csv", index=False)
