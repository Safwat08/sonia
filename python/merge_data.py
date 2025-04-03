# Imports
import pandas as pd
import scanpy as sc
import anndata as ad
from functions import save_ann
import random
import numpy as np
import argparse
import sys
import os

# Set seed for reproducibility
seed = 42
random.seed(seed)
np.random.seed(seed)

if __name__ == "__main__":
    if len(sys.argv) > 3:
        # Parse command line arguments
        merge_object_name = sys.argv[1]          # Output filename
        merge_label_name = sys.argv[2]           # Label for merged objects (e.g. "batch")
        paths_of_objects_to_merge = sys.argv[3:] # Paths to directories containing 'adata.h5ad'

        print("Merge object name:", merge_object_name)
        print("Merge label name:", merge_label_name)
        print("Paths of objects to merge:", paths_of_objects_to_merge)
    else:
        print("Usage: python script.py <output_name> <merge_label> <paths...>")
        sys.exit(1)

merge_dict = {}

# Load each AnnData object from its respective path
for path in paths_of_objects_to_merge:
    adata_path = os.path.join(path, "adata.h5ad")

    try:
        adata = sc.read_h5ad(adata_path)
        print(f"Successfully loaded data from {adata_path}")
        adata_name = path.replace("/", "_")   
        # Make the first column (gene_ids) the index of adata.var
        adata.var.set_index(adata.var['gene_ids'], inplace=True)

        adata.X = adata.layers['rawcounts'].copy()
        
        merge_dict[adata_name] = adata  # Use path as key for labeling
    except Exception as e:
        print(f"Error loading data from {adata_path}: {e}")
        sys.exit(1)

# Merge all loaded AnnData objects
try:
    merged_adata = ad.concat(
            merge_dict,
            label=merge_label_name,
            index_unique="_"
        )
    print("Merge successful!")
except Exception as e:
    print(f"Error during merging: {e}")
    sys.exit(1)

adata1 = list(merge_dict.values())[0]

# Use first adata to get gene names
merged_adata.var['gene_names'] = adata1.var['gene_names'] # Get the gene names from another anndata

# Make gene_ids column
merged_adata.var['gene_ids'] = merged_adata.var.index # Make a new column for ensemble ids

# Set rawcounts layer
merged_adata.layers['rawcounts'] = merged_adata.X.copy()

# Normalizing to 10,000 total counts
sc.pp.normalize_total(merged_adata)

# Logarithmize the data
sc.pp.log1p(merged_adata)

# Save log1p layer
merged_adata.layers['log1p'] = merged_adata.X.copy()

# Optional: save merged object
save_ann(merged_adata, 'output', merge_object_name, full = True)