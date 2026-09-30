""" Create a GO ID graph depicting most common GO IDs a given GMS cluster """

import argparse
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
cluster_df = pd.read_csv(str(cluster_df_path))


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




### Join dataframes -------------------------------------------------------

# Get dataframe containing GO ID, label, namespace, and GO module
go_metadata_df = go_metadata_df.merge(go_module_df, on = "go_id", how = "inner")

# Get dataframe containing only rows present in cluster dataframe, and only GO
# ID columns where GO ID is present in at least one of these rows
cluster_df = cluster_df.merge(go_presence_df, on = "uniprot", how = "inner")
go_cols = [c for c in cluster_df.columns if c.startswith("GO:")]
keep_cols = [c for c in go_cols if cluster_df[c].any()] # >=1 non-zero val
cluster_go_df = cluster_df.drop(
  columns=[c for c in cluster_df.columns if c not in keep_cols]
)







