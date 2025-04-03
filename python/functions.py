def change_cluster_name(adata, 
                        old: str, 
                        new: str, 
                        orig_meta_name: str, 
                        meta_name: str = "cluster"):
    """
    Change the name of a cluster in an AnnData object.

    Parameters:
    - adata: AnnData object where the cluster names are to be updated.
    - old: Character string specifying the old cluster ID.
    - new: Character string specifying the new cluster name.
    - orig_meta_name: The name of the existing column in adata.obs containing the original cluster names.
    - meta_name: The name of the new column to be added to adata.obs with updated cluster names.
    - inc_meta: If True, a different approach is taken to include the cluster name in the metadata.

    Returns:
    - Updated AnnData object.
    """
    import warnings
    # Make sure to call copy as a method with parentheses
    meta = adata.obs[orig_meta_name].copy()
    
    # Iterate over the old cluster names to update the relevant entries directly in meta
    for i in range(len(old)):
        # Get the index locations where the value matches old[i]
        old_idx = meta[meta == old[i]].index
        
        if not old_idx.empty:
            # Replace the values at the identified indices
            meta.iloc[old_idx] = new
        else:
            # Raise a warning if the old cluster name is not found
            warnings.warn(f"Cluster '{old[i]}' is not available in the current metadata '{orig_meta_name}' and will be ignored.", UserWarning)
    
    # Assign the updated metadata to the new column
    adata.obs[meta_name] = meta
        
def is_outlier(adata, 
               metric: str, 
               nmads: int):
    """
    Identify outliers in the AnnData object based on the specified metric and 
    the number of median absolute deviations (MADs). 
    Anna Schaar. Theiss Lab

    Parameters:
    -----------
    adata : AnnData
        The AnnData object containing the observation data.
    metric : str
        The name of the column in `adata.obs` on which to detect outliers.
    nmads : int
        The number of median absolute deviations (MADs) to use as the threshold 
        for determining outliers.

    Returns:
    --------
    outlier : pandas.Series
        A boolean Series where `True` indicates that the corresponding 
        observation is an outlier.
    """
    import numpy as np 
    from scipy.stats import median_abs_deviation # For QC
    
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def preprocess_adata(path_or_adata: str,
                    dataset: str,
                    species: str = "Mouse"):
    """
    Load an h5 file using the Scanpy library and modify the AnnData object.
    
    Parameters:
    -----------
    path_h5 : str
        The file path to the .h5 file to be loaded.
    
    Returns:
    --------
    adata : AnnData
        The processed AnnData object where:
        - A new column 'gene_names' is created from the index of `adata.var`.
        - The first column (gene_ids) of `adata.var` is set as the new index.
        - percent.mt and other qc metrics are added and object is filtered
    """
    # process_ann_data(path_h5)
    import scanpy as sc  # Import scanpy at the top?
    import pandas as pd

    # Load the .h5 file or use the passed AnnData object
    if isinstance(path_or_adata, str):
        adata = sc.read_10x_h5(path_or_adata)
    elif isinstance(path_or_adata, sc.AnnData):
        adata = path_or_adata
    else:
        raise ValueError("path_or_adata must be either a file path (str) or an AnnData object.")
        
    # Insert a column 'gene_names' with the values from the index of adata.var
    adata.var['gene_names'] = adata.var.index
    
    # Make the first column (gene_ids) the index of adata.var
    adata.var.set_index(adata.var['gene_ids'], inplace=True)

    # Handle species-specific mitochondrial gene detection
    if species.lower() == "human":
        mito_prefix = "MT-"
    elif species.lower() == "mouse":
        mito_prefix = "mt-"
    else:
        raise ValueError("Species not recognized. Use 'Human' or 'Mouse'.")

    # Identify mitochondrial genes based on species-specific prefix
    adata.var["mt"] = adata.var["gene_names"].str.startswith(mito_prefix)
    
    # Calculate qc_metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, log1p=True, percent_top=[20]
    )
    
    # Filter by 5 Median Absolute Deviations (MADs)
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    
    # Filter pct_counts_mt at 3 MADs and and > 8
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 5) | (
    adata.obs["pct_counts_mt"] > 8)

    # Makea a column with the name of the dataset
    adata.obs['dataset'] = dataset

    # Make a new layer 'rawcounts' to keep the raw data
    adata.layers['rawcounts'] = adata.X.copy()
    
    # Return the modified AnnData object
    return adata

def save_ann(adata, path: str, name: str, full = True):
    """
    Save an AnnData object to disk in a specified format.
    
    This function creates a folder with the specified name at the provided path.
    Inside the folder, it saves:
    - `adata.obs` as 'barcodes.csv'
    - `adata.var` as 'features.csv'
    - Each layer in `adata.layers` as '{layer_name}_matrix.mtx'
    - Each layer in `adata.obsm` and 'adata.varm' as '{layer_name}_obsm.csv'/'_varm.csv'
    - The entire AnnData object as 'adata.h5ad'
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object to be saved.
    path : str
        The directory path where the folder will be created.
    name : str
        The name of the folder to be created.
    
    Returns:
    --------
    None
    """
    import os
    import pandas as pd
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix
    
    # Create the directory if it doesn't exist
    folder_path = os.path.join(path, name)
    os.makedirs(folder_path, exist_ok=True)

    # Save the entire AnnData object as adata.h5ad
    adata_file = os.path.join(folder_path, 'adata.h5ad')
    adata.write(adata_file)

    if full:
        # Save adata.obs as barcodes.csv
        barcodes_file = os.path.join(folder_path, 'barcodes.csv')
        adata.obs.to_csv(barcodes_file, index=True, header=True)
        
        # Save adata.var as features.csv
        features_file = os.path.join(folder_path, 'features.csv')
        adata.var.to_csv(features_file, index=True, header=True)
        
        # Save each layer in adata.layers as a separate .mtx file
        for layer_name in adata.layers.keys():
            layer_data = adata.layers[layer_name]
            
            # Convert to CSR matrix if it's not already in sparse format
            if not isinstance(layer_data, csr_matrix):
                layer_data = csr_matrix(layer_data)
            
            # Define the output file path for each layer
            layer_matrix_file = os.path.join(folder_path, f"{layer_name}_matrix.mtx")
            
            # Write the CSR matrix to a Matrix Market file
            mmwrite(layer_matrix_file, layer_data)
            print(f"Saved {layer_name} as a Matrix Market file at {layer_matrix_file}")
    
        # Save each component in adata.obsm as CSV
        for obsm_name in adata.obsm.keys():
            obsm_data = adata.obsm[obsm_name]  # Corrected to adata.obsm
            
            obsm_csv_file = os.path.join(folder_path, f"{obsm_name}_obsm.csv")
            
            pd.DataFrame(obsm_data).to_csv(obsm_csv_file, index=True)
            
            print(f"Saved {obsm_name} (obsm) as CSV at {obsm_csv_file}")
        
        # Save each component in adata.varm as CSV
        for varm_name in adata.varm.keys():
            varm_data = adata.varm[varm_name]
            
            varm_csv_file = os.path.join(folder_path, f"{varm_name}_varm.csv")
            
            pd.DataFrame(varm_data).to_csv(varm_csv_file, index=True)
            
            print(f"Saved {varm_name} (varm) as CSV at {varm_csv_file}")

def dim_reduc(adata, impute: bool, integrate: bool, scvi: bool, batch=None):
    """
    Perform a series of preprocessing and analysis steps on the input AnnData object using Scanpy.
    
    This function:
    1. Performs PCA on highly variable genes.
    2. Optional: Conducts imputation based on diffusion map
    2. Creates a PCA scatter plot.
    3. Computes nearest neighbors, UMAP, and Leiden clustering with multiple resolutions.
    4. Visualizes plots colored by Leiden clustering at various resolutions.
    
    Parameters:
    -----------
    adata : AnnData
        The input AnnData object to be analyzed.
    impute : bool
        Whether to impute using MAGIC or not
    
    Returns:
    --------
    None
    """

    if integrate and batch is None:
        raise ValueError("The 'batch' parameter must be provided when 'integrate' is True.")

    # Necessary imports for the function
    import scanpy as sc
    import scanpy.external as sce
    import palantir
    import networkx as nx
    import numpy as np
    from fa2_modified import ForceAtlas2

    # Perform PCA on the highly variable genes
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_deviant_genes")

    if integrate:
        sce.pp.harmony_integrate(adata, batch)

    if impute:
        if integrate:
            palantir.utils.run_diffusion_maps(adata, n_components=5, pca_key = 'X_pca_harmony')
        else:
            palantir.utils.run_diffusion_maps(adata, n_components=5, pca_key = 'X_pca')
            
        palantir.utils.determine_multiscale_space(adata)
        palantir.utils.run_magic_imputation(adata)
        adata.X = adata.layers['MAGIC_imputed_data'].copy()
        
        # Create ForceAtlas2 object with desired parameters
        forceatlas2 = ForceAtlas2(
                                  # Behavior alternatives
                                  outboundAttractionDistribution=True,  # Dissuade hubs
                                  linLogMode=False,  # NOT IMPLEMENTED
                                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                                  edgeWeightInfluence=1.0,
        
                                  # Performance
                                  jitterTolerance=1.0,  # Tolerance
                                  barnesHutOptimize=True,
                                  barnesHutTheta=1.2,
                                  multiThreaded=False,  # NOT IMPLEMENTED
        
                                  # Tuning
                                  scalingRatio=2.0,
                                  strongGravityMode=False,
                                  gravity=1.0,
        
                                  # Log
                                  verbose=True)
        # Perform PCA on the highly variable genes
        sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_deviant_genes")

        # Extract the KNN graph from Scanpy
        knn_graph = adata.obsp['DM_Kernel']  # Sparse adjacency matrix of the KNN graph
        
        # Convert the sparse matrix to a NetworkX graph
        G = nx.from_scipy_sparse_array(knn_graph)
        
        # Use NetworkX to create a force-directed layout (spring layout)
        pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=1000)
        
        # Convert positions into numpy array for storing in adata
        coordinates = np.array(list(pos.values()))
        
        # Store the coordinates in adata.obsm under the key 'X_force'
        adata.obsm['X_force'] = coordinates

    if scvi:
        adata.var['highly_variable_genes'] = adata.var['highly_deviant_genes']
        scvi.model.SCVI.setup_anndata(adata, layer="rawcounts")
        model = scvi.model.SCVI(adata)
        
    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors = 15)
    sc.tl.umap(adata)
    
    # Perform Leiden clustering with different resolutions
    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

    # Perform phate
    sc.external.tl.phate(adata)

    return(adata)