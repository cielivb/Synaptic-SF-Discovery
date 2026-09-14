import pandas as pd
from pathlib import Path
import networkx as nx
import json


def extract_go(uri):
    """ Extract GO ID out of GO URI """
    return uri.split("/")[-1].replace("GO_", "GO:")


# Load synaptic GO IDs
synaptic_go_path = "results/1_filter_by_func/goslim_synapse_go_ids.csv"
synaptic_go_ids = pd.read_csv(synaptic_go_path)

# Load GO-basic JSON into directed graph
go_graph = nx.DiGraph()
with open("data/go-basic.json") as f:
    go = json.load(f)
for edge in go["graphs"][0]["edges"]: # Add child -> parent edges
    child, parent = edge["sub"], edge["obj"]
    go_graph.add_edge(extract_go(child), extract_go(parent))
for node in go["graphs"][0]["nodes"]: # Add node metadata
    G.add_node(**node)


# Map synaptic GO IDs to their ancestor IDs
synaptic_graph = go_graph.subgraph(list(synaptic_go_ids[["go_id"]])).copy()

