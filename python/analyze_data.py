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
        genes_path = os.path.join(path_name, "deviant_genes.csv")
        
        # Check if the AnnData file exists
        if os.path.exists(file_path):
            adata = sc.read_h5ad(file_path)  # Load the AnnData object
            
            # Check if the deviant genes file exists
            if os.path.exists(genes_path):
                deviant_genes = set(pd.read_csv(genes_path, index_col=0).iloc[:, 0])
                adata.var['highly_deviant_genes'] = adata.var.index.isin(deviant_genes)
                adata.X = adata.layers['log1p'].copy()
	
                # Conduct scanpy pipeline: log1p-norm > PCA > kNN > umap, option Impute
                adata = dim_reduc(adata, impute = True, integrate = False, scvi = False)

                adata.var.index = adata.var['gene_names'].copy()

                # output_file_name
                output_file_name = adata_name + "_analyzed"

                # Save anndata
                save_ann(adata, "output", output_file_name, full = True)

                if 'batch' in adata.obs.columns:
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color = "batch",
                        frameon=False, 
                        cmap = 'turbo'
                    )
                    
                if 'sample' in adata.obs.columns:
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color = "sample",
                        frameon=False, 
                        cmap = 'turbo'
                    )


                # Contamination
                sc.pl.embedding(
                    adata,
                    basis="force",
                    layer="log1p",
                    color = 'leiden_res0_25',
                    legend_loc="on data",
                    frameon=False, 
                    cmap = 'turbo'
                )

                print(f"adata_name: {adata_name}")  # To check the value of adata_name
                
                if adata_name == 'fm_ilc_neg_postqc/':
                    # Define the gene lists
                    islet_genes = ['TOP2A','REG1A','CPA1','GP2', 'RBJL', 'CELA3A', 'RBFOX1', 'SUSD4', 'EGF', 'SOX6', 'ONECUT1', 'LGR4', 'SCTR', 'SLC4A4', 'CFTR', 'DCDC2', 'ANXA4', 'FGFR2', 'GLIS3', 'RFX6', 'GCG', 'CHGA', 'INS', 'SLC30A8', 'ROBO2', 'SST', 'ST18']
                    factors = ['VEGFA', 'VEGFB', 'VEGFC', 'AREG', 'COL18A1', 'DPP4', 'TIMP1', 'PDGF', 'ANG', 'IGFBP2', 'IGFPB3', 'SERPINE1', 'THBS1', 'PIGF']
                    
                    # Filter genes that are present in the adata object
                    islet_genes_present = [gene for gene in islet_genes if gene in adata.var_names]
                    factors_present = [factor for factor in factors if factor in adata.var_names]
                    
                    # Contamination: Plot only the genes that are present in adata
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color=islet_genes_present,
                        frameon=False, 
                        cmap='turbo'
                    )
                    
                    # Contamination: Plot only the factors that are present in adata
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color=factors_present,
                        frameon=False, 
                        cmap='turbo'
                    )
                    
                else:                 
                    vascular_hierarchy_genes = ['Chd5', 'Pecam1', 'Vwf', 'Mgp', 'Gja4','Efnb2','Efhd1','Cd36','Esm1', 'Ephb4','Nrp2','Top2a','Plvap', 'Ttn', 'Sele', 'Mki67']
                    contamination_genes = ['Ins2','Nkx2-2','Gcg','Chga', 'Ghr', 'Celf2', 'Zeb2']
                
                    # vascular hierarchy genes
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color = vascular_hierarchy_genes,
                        frameon=False, 
                        cmap = 'turbo'
                    )
                        
                    # Contamination
                    sc.pl.embedding(
                        adata,
                        basis="force",
                        layer="log1p",
                        color = contamination_genes,
                        frameon=False, 
                        cmap = 'turbo'
                    )
                    
                print(f"Processed {genes_path} and updated {file_path}.")
            else:
                print(f"Error: {genes_path} not found.")
        else:
            print(f"Error: {file_path} not found.")
    else:
        print(f"Error: Directory {path_name} does not exist.")

