""" Create a GO ID graph depicting most common GO IDs a given GMS cluster """

import argparse
import json
import networkx as nx
import pandas as pd
from pathlib import Path


# Read in cluster ID
parser = argparse.ArgumentParser()
parser.add_argument("cluster_id", type = int)
args = parser.parse_args()
CLUSTER_ID = args.cluster_id




### Read in and wrangle initial/pre-join dataframes -----------------------

# Read in cluster dataframe
cluster_df_path = Path(f"results/5_make_visual_summaries/{CLUSTER_ID}.csv")
CLUSTER_DF = pd.read_csv(str(cluster_df_path))


# Read in GO ID metadata dataframe (maps GO IDs to GO names/descriptions)
go_metadata_path = Path(f"results/1_filter_proteins/go_metadata.csv")
go_metadata_df = pd.read_csv(str(go_metadata_path))

def extract_namespace(meta_str):
  """ The namespace is stuck inside a string representation of a python dict """
    try:
        meta = ast.literal_eval(meta_str)
        bpv = meta.get("basicPropertyValues", [])
        if bpv:
            return bpv[0].get("val")
        return None
    except Exception:
        return None

go_metadata_df["namespace"] = go_metadata_df["meta"].apply(extract_namespace)
go_metadata_df = go_metadata_df[["go_id", "lbl", "namespace"]]
go_metadata_df = go_metadata_df.rename({"label": "lbl"})

# Read in GO module members df(maps GO IDs to their GO modules/GO ID clusters)
go_module_df_path = Path(
  "results/3_extract_functional_features/go_module_members.csv"
)
go_module_df = pd.read_csv(str(go_module_df_path))

# Read in GO ID presence/absence per protein dataframe
go_presence_df_path = Path(
  "results/3_extract_functional_features/functional_features.csv"
)
go_presence_df = pd.read_csv(str(go_presence_df_path))




### Join dataframes -------------------------------------------------------

# Get dataframe containing only proteins present in cluster dataframe, and GO
# IDs present in at least one of these proteins






