import scanpy as sc
import numpy as np

# Load the .h5 file using the read_10x_h5 function from Scanpy
adata_mixed = sc.read_10x_h5('rawdata/fm_ilc_neg_mixed/filtered_feature_bc_matrix.h5')

adata_human = sc.read_10x_h5('rawdata/fm_ilc_neg_human/filtered_feature_bc_matrix.h5')

# Create boolean masks for human and mouse genes
human_mask = adata.var['genome'] == 'GRCh38'
mouse_mask = adata.var['genome'] == 'GRCm39'

# Sum counts for human and mouse genes in each cell
human_counts = np.sum(adata.X[:, human_mask], axis=1)  # Summing across features
mouse_counts = np.sum(adata.X[:, mouse_mask], axis=1)

# Total counts in each cell
total_counts = human_counts + mouse_counts

# Calculate proportions
human_proportion = human_counts / total_counts
mouse_proportion = mouse_counts / total_counts

# Add proportions to `adata.obs` for each cell
adata.obs['human_proportion'] = human_proportion
adata.obs['mouse_proportion'] = mouse_proportion


# Identify indices of cells with human_proportion > 0.8
high_human_proportion_cells = adata.obs['human_proportion'] > 0.8

# Subset the AnnData object to include only these cells
adata_human_subset = adata[high_human_proportion_cells, :]