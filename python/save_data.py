# Import
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

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Print all arguments except the script name
        print("All arguments:", sys.argv[1:])
        
        # Store the arguments as a list
        adata_dir = sys.argv[1:]
        print("adata_dir:", adata_dir)
    else:
        print("No arguments passed!")
        adata_dir = []  # Define an empty list if no arguments are passed


adata_dict = {}

# Load the .h5ad file using the read_h5ad function from Scanpy
for path_name in adata_dir:
    print(f"Working on path {path_name}.")
    
    # Remove the "output/" prefix from the path name
    adata_name = path_name.replace("output/", "")
    
    # Check if the directory exists
    if os.path.isdir(path_name):
        file_path = os.path.join(path_name, "adata.h5ad")
        
        # Check if the AnnData file exists
        if os.path.exists(file_path):
            adata = sc.read_h5ad(file_path)  # Load the AnnData object
            
            # output_file_name
            output_file_name = adata_name

            # Save anndata
            save_ann(adata, "output", output_file_name, full = True)