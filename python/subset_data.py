# Imports
import pandas as pd
import scanpy as sc
from functions import dim_reduc, save_ann
import random
import numpy as np
import sys
import os

# Set seed
seed = 42
random.seed(seed)
np.random.seed(seed)

# Argument parsing
if __name__ == "__main__":
    if len(sys.argv) > 2:
        # First argument: whether to keep clusters ("true" or "false")
        keep_arg = sys.argv[1].lower()
        if keep_arg not in ['true', 'false']:
            print("First argument must be 'true' or 'false'.")
            sys.exit(1)

        keep = keep_arg == "true"
        print("Keep clusters:", keep)

        # Second argument: path to AnnData folder (assumes adata.h5ad inside it)
        path = sys.argv[2]
        adata_path = os.path.join(path, "adata.h5ad")
        print("Loading from path:", adata_path)

        # Remaining arguments: clusters to filter
        strings = sys.argv[3:]
        print("Filtering clusters:", strings)

    else:
        print("Usage: python script.py <true|false> <path_to_folder> <cluster1> <cluster2> ...")
        sys.exit(1)

    # Load the .h5ad file
    try:
        adata = sc.read_h5ad(adata_path)
        print(f"Successfully loaded data from {adata_path}")
    except Exception as e:
        print(f"Error loading data from {adata_path}: {e}")
        sys.exit(1)

    # Ensure column exists
    cluster_col = "leiden_res0_25"
    if cluster_col not in adata.obs.columns:
        print(f"Error: The column '{cluster_col}' does not exist in adata.obs.")
        sys.exit(1)

    # Subset based on keep/remove logic
    if keep:
        subset_adata = adata[adata.obs[cluster_col].isin(strings), :]
    else:
        subset_adata = adata[~adata.obs[cluster_col].isin(strings), :]

    if keep:
        # Check that only the specified clusters are present
        unexpected_clusters = set(subset_adata.obs[cluster_col].unique()) - set(strings)
        if unexpected_clusters:
            print(f"Error: The following unexpected clusters are present in '{cluster_col}': {unexpected_clusters}")
            sys.exit(1)
        else:
            print("Only the specified clusters are present in the subset.")
    else:
        # Check that none of the specified clusters are present
        still_present = set(strings).intersection(subset_adata.obs[cluster_col].unique())
        if still_present:
            print(f"Error: The following entries are still present in '{cluster_col}': {still_present}")
            sys.exit(1)
        else:
            print("None of the provided clusters are present anymore.")

# Rawcounts
subset_adata.X = subset_adata.layers['rawcounts'].copy()

# Normalizing to 10,000 total counts
sc.pp.normalize_total(subset_adata)

# Logarithmize the data
sc.pp.log1p(subset_adata)

subset_adata.layers['log1p'] = subset_adata.X.copy()

# Save the subset data to a new file (optional)

adata_name = path.replace('output/', '')

save_path = adata_name + "_subset"

save_ann(subset_adata, 'output', save_path)

print(f"Subset data saved to {save_path}")

run = False

if run:
    
    # Load the .h5 file using the read_10x_h5 function from Scanpy
    adata = sc.read_h5ad(path)
    
    subset_adata = adata[~adata.obs['leiden_res0_25'].isin(strings), :]
    
    # Load deviant genes
    deviant_genes = set(pd.read_csv('output/subset_raw/deviant_genes.csv', index_col=0).iloc[:,0])
    
    # Create a boolean column for whether a gene is highly deviant or not
    adata.var['highly_deviant_genes'] = adata.var.index.isin(deviant_genes)
    
    adata.X = adata.layers['log1p'].copy()
    	
    # Conduct scanpy pipeline: log1p-norm > PCA > kNN > umap, option Impute
    adata = dim_reduc(adata, impute = True)
    
    #
    adata.write('output/subset_analyzed/adata.h5ad')
    
    # Subset adata
    subset_adata = adata[~adata.obs['leiden_res0_25'].isin(["7", "5", "13", "11", "8", "14"]), :]
    
    # Rawcounts
    subset_adata.X = subset_adata.layers['rawcounts'].copy()
    
    # Normalizing to 10,000 total counts
    sc.pp.normalize_total(subset_adata)
    
    # Logarithmize the data
    sc.pp.log1p(subset_adata)
    
    subset_adata.layers['log1p'] = subset_adata.X.copy()
    
    subset_adata.write('output/subset2_raw/adata.h5ad')