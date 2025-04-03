# Import
from functions import preprocess_adata, save_ann
import os
import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np

# Set seed
seed = 42
random.seed(seed)
np.random.seed(seed)

# Check consistent seed
print(random.randint(1, 100))
print(np.random.rand(3))           

# Path to the rawdata directory
rawdata_directory = "rawdata"
output_directory = 'output'
figures_directory = 'figures'

# Process negative fraction seperately

print(f"Loading adata for negative fraction seperately")

# Load the .h5 file using the read_10x_h5 function from Scanpy
adata_mixed = sc.read_10x_h5('rawdata/fm_ilc_neg_mixed/filtered_feature_bc_matrix.h5')
adata_human = sc.read_10x_h5('rawdata/fm_ilc_neg_human/filtered_feature_bc_matrix.h5')

# Create boolean masks for human and mouse genes
human_mask = adata_mixed.var['genome'] == 'GRCh38'
mouse_mask = adata_mixed.var['genome'] == 'GRCm39'

# Sum counts for human and mouse genes in each cell
human_counts = np.sum(adata_mixed.X[:, human_mask], axis=1)  # Summing across features
mouse_counts = np.sum(adata_mixed.X[:, mouse_mask], axis=1)

# Total counts in each cell
total_counts = human_counts + mouse_counts

# Calculate proportions
human_proportion = human_counts / total_counts
mouse_proportion = mouse_counts / total_counts

# Add proportions to `adata.obs` for each cell
adata_mixed.obs['human_proportion'] = human_proportion
adata_mixed.obs['mouse_proportion'] = mouse_proportion

# Identify indices of cells with human_proportion > 0.8
high_human_proportion_cells = adata_mixed.obs[adata_mixed.obs['human_proportion'] > 0.8].index

# Identify indices of cells with human_proportion > 0.8
high_mouse_proportion_cells = adata_mixed.obs[adata_mixed.obs['mouse_proportion'] > 0.8].index

# Ambiguous indices
ambiguous_cells = adata_mixed.obs[
    (adata_mixed.obs['mouse_proportion'] < 0.8) & (adata_mixed.obs['human_proportion'] < 0.8)
].index

# Find matching indices for high human proportion cells
matching_human_indices = high_human_proportion_cells.intersection(adata_human.obs.index)

# Find matching indices for high mouse proportion cells
matching_mouse_indices = high_mouse_proportion_cells.intersection(adata_human.obs.index)

# Find matching indices for ambiguous cells
matching_ambiguous_indices = ambiguous_cells.intersection(adata_human.obs.index)

# Initialize the filtering_results column with "not_matched" as default
adata_human.obs['filtering_results'] = 'not_matched'

# Label the high human proportion cells
adata_human.obs.loc[matching_human_indices, 'filtering_results'] = 'high_human'

# Label the high mouse proportion cells
adata_human.obs.loc[matching_mouse_indices, 'filtering_results'] = 'high_mouse'

# Label the ambiguous cells
adata_human.obs.loc[matching_ambiguous_indices, 'filtering_results'] = 'ambiguous'

# Preprocess data
adata_human = preprocess_adata(adata_human, dataset='fm_ilc_neg', species= 'Human')

# Dictionary to store the names of objects
adata_dict = {}

print(f"Loading adata from rawdata")

# List dir except neg fraction
data_dirs = [d for d in os.listdir(rawdata_directory) if d not in ["fm_ilc_neg_mixed", "fm_ilc_neg_human"]]

# Iterate over each folder
for folder_name in data_dirs:
    print(f"Working on folder {folder_name}.")
    folder_path = os.path.join(rawdata_directory, folder_name)
    if os.path.isdir(folder_path):
        file_path = os.path.join(folder_path, "filtered_feature_bc_matrix.h5")
        if os.path.exists(file_path):
            # Generate a valid variable name
            adata_name = folder_name.replace('-', '_').replace(' ', '_')
            
            # Process the data
            adata = preprocess_adata(file_path, dataset=adata_name)
            
            # Store the processed object in the dictionary
            adata_dict[adata_name] = adata
            
            print(f"Variable {adata_name} created and stored in the dictionary.")

adata_dict['fm_ilc_neg'] = adata_human

print(f"Saving pre-qc adata")

# Make data frame of outlier 
qc_results = pd.DataFrame(columns=["Dataset", "Total", "Remaining", "Filtered_Mitochondrial_Only", "Filtered_OtherMetrics_Only", "Filtered_Both"])

# Save preqc adata
for adata_name, adata in adata_dict.items():
    print(f"Working on adata {adata_name}.")
    
    # Make output file path
    output_file_name = adata_name + '_preqc'

    # Save the AnnData object using save_ann
    try:
        save_ann(adata, output_directory, output_file_name, full = True)
        print(f"Saved {adata_name} to {output_directory}.")
    except Exception as e:
        print(f"Error saving {adata_name}: {e}")

    # Get values
    total = adata.shape[0]
    remaining = ((adata.obs["mt_outlier"] == False) & (adata.obs["mt_outlier"] == False)).sum()
    mt_filtered = ((adata.obs["mt_outlier"] == True) & (adata.obs["outlier"] == False)).sum()
    other_filtered = ((adata.obs["outlier"] == True) & (adata.obs["mt_outlier"] == False)).sum()
    both_filtered = ((adata.obs["mt_outlier"] == True) & (adata.obs["mt_outlier"] == True)).sum()

    # Make row
    adata_row = [adata_name, total, remaining, mt_filtered, other_filtered, both_filtered]

    # Add row
    qc_results.loc[len(qc_results)] = adata_row

    # Save plots
    outlier_scatter_file_name = "_" + output_file_name + "_qcplot1_outlier.png"
    mt_outlier_scatter_file_name = "_" + output_file_name + "qcplot2_mt_outlier.png"
    outlier_violin_file_name = "_" + output_file_name + "qcplot3_outlier.png"
    mt_outlier_violin_file_name = "_" + output_file_name + "qcplot4_mt_outlier.png"


    # Change boolean to categorical for easier visualization
    adata.obs["outlier"] = pd.Categorical(adata.obs["outlier"])
    adata.obs["mt_outlier"] = pd.Categorical(adata.obs["mt_outlier"])
    
    # QC Plot 1: Scatter of ngenes by count & total counts, colored by outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="outlier",
        save=outlier_scatter_file_name
    )

    # QC Plot 2: Scatter of ngenes by count & total counts, colored by mt_outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="mt_outlier",
        save=mt_outlier_scatter_file_name
    )

    # QC Plot 3: Violin plot of n genes by count, total counts and pct_counts_mt, colored by outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'outlier',
        save = outlier_violin_file_name
    )

    # QC Plot 4: Violin plot of n genes by count, total counts and pct_counts_mt, colored by mt_outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'mt_outlier',
        save = mt_outlier_violin_file_name
    )

    if adata_name == 'fm_ilc_neg':
        # QC Plot 1: Scatter of ngenes by count & total counts, colored by outlier
        sc.pl.scatter(
            adata_human,
            x="n_genes_by_counts", 
            y="total_counts",     
            color="filtering_results",
            save = adata_name + "_qcplot5_species.png"
        )
        
        # QC Plot 3: Violin plot of n genes by count, total counts and pct_counts_mt, colored by outlier
        sc.pl.violin(
            adata_human,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True, groupby= 'filtering_results',
            save = adata_name + "_qcplot6_species.png"
        )

    # Change back to bool
    adata.obs["outlier"] = adata.obs["outlier"].astype(bool)
    adata.obs["mt_outlier"] = adata.obs["mt_outlier"].astype(bool)


# Save qc_results
qc_results.to_csv('output/qc_results.csv', index=False)

print(f"QC Results:")

print(qc_results)

# Loop through each row to create and save a pie chart
for index, row in qc_results.iterrows():
    # Extract the dataset name from the first column
    dataset_name = row["Dataset"]
    
    # Extract the last 3 columns for the row
    data = row.iloc[-4:]  # Get last 3 columns
    labels = data.index   # Use column names as labels
    
    # Plot the pie chart
    plt.figure()
    plt.pie(data, labels=labels, autopct="%1.1f%%", startangle=90)
    plt.title(f"Pie Chart for {dataset_name}")
    
    # Save the pie chart to the "figures" directory
    file_name = f"pie_{dataset_name}_qcplot5.png"
    plt.savefig(os.path.join(figures_directory, file_name))
    plt.show() # Show plot
    plt.close()  # Close the figure to free memory

    print(f"Saved pie chart for {dataset_name}")

print(f"Saving post-qc adata")

# Dictionary to store the names of objects
adata_filtered_dict = {}

# Save postqc adata
for adata_name, adata in adata_dict.items():
    print(f"Working on adata {adata_name}.")
    
    # Filter based on outlier
    adata_filtered = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    if adata_name == 'fm_ilc_neg':
        adata_filtered = adata_filtered[adata_filtered.obs['filtering_results'] == 'high_human', :]

    # Store the processed object in the dictionary
    adata_filtered_dict[adata_name] = adata_filtered
    
    # Make output file path
    output_file_name = adata_name + '_postqc'

    # Normalizing to 10,000 total counts
    sc.pp.normalize_total(adata_filtered)
    
    # Logarithmize the data
    sc.pp.log1p(adata_filtered)
    
    adata_filtered.layers['log1p'] = adata_filtered.X.copy()

    # Save the AnnData object using save_ann
    try:
        save_ann(adata_filtered, output_directory, output_file_name)
        print(f"Saved {adata_name} to {output_directory}.")
    except Exception as e:
        print(f"Error saving {adata_name}: {e}")

# Merge BM
adata_bm_merge = ad.concat(
    [adata_filtered_dict["bm_ilc"], adata_filtered_dict["bm_only"]],
    label="sample",
    index_unique="-",
    keys=["bm_ilc", "bm_only"],  # Custom labels for batches
)

# Merge FM
adata_fm_merge = ad.concat(
    [adata_filtered_dict["fm_ilc_pos"], adata_filtered_dict["fm_only_pos"]],
    label="sample",
    index_unique="-",
    keys=["fm_ilc", "fm_only"],  # Custom labels for batches
)

# Create merged dict
merge_dict = {
    "bm": adata_bm_merge,
    "fm": adata_fm_merge,
}

# Get one of the objects for modifying var names
adata1 = list(adata_filtered_dict.values())[0]

# Save all merges
for adata_name, adata_merge in merge_dict.items():
    print(f"Working on merged adata for {adata_name}.")

    # Use first adata to get gene names
    adata_merge.var['gene_names'] = adata1.var['gene_names'] # Get the gene names from another anndata

    # Make gene_ids column
    adata_merge.var['gene_ids'] = adata_merge.var.index # Make a new column for ensemble ids

    # Set rawcounts layer
    adata_merge.layers['rawcounts'] = adata_merge.X.copy()
    
    # Normalizing to 10,000 total counts
    sc.pp.normalize_total(adata_merge)
    
    # Logarithmize the data
    sc.pp.log1p(adata_merge)

    # Save log1p layer
    adata_merge.layers['log1p'] = adata_merge.X.copy()

    # Make output file path
    output_file_name = adata_name + '_merge'

    # save ann
    save_ann(adata_merge, output_directory, output_file_name)