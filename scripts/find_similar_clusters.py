""" Find Similar Clusters

"Cluster" similar clusters (of protein geometries) by finding connected
components of a graph where nodes are cluster labels and edges represent
being similar (i.e., not being statistically different; alpha >= 0.05)

"""

import argparse
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, required=True)
parser.add_argument("--title", type=str, required=True)
args = parser.parse_args()
print(f"Using file {args.file}")

# Load pairwise permanova results and filter for "similar"
# pairs (pval_adj > alpha)
df = pd.read_csv(args.file)
alpha = 0.05
similar = df[df["pval_adj"] > alpha]

# Iteratively build networkx graph then extract connected
# components into list of python sets
G = nx.Graph()
for _, row in similar.iterrows():
    G.add_edge(row["cluster_1"], row["cluster_2"])
components = list(nx.connected_components(G))

# Visualise connected components
pos = nx.spring_layout(G, seed=254)
plt.figure(figsize=(8, 8))
nx.draw(
    G,
    pos,
    with_labels=True,
    node_color="skyblue",
    edge_color="gray",
    node_size=1200,
    font_size=12,
    linewidths=1.5
)
plt.title(args.title)
outpath = Path(args.file).parent / f"{args.title.replace(" ", "_")}.svg"
plt.savefig(str(outpath), format="svg")
