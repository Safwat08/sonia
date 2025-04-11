import decoupler as dc
import scanpy as sc
from python.functions import save_ann
import pandas as pd

adata = sc.read_h5ad('output/fm_bm_fresh_merge_analyzed/adata.h5ad')

# Step 1: Set 'Arterial' → 'Capillary 3' for specific datasets
adata.obs.loc[
    (adata.obs['vascular_subtype'] == 'Arterial') &
    (adata.obs['dataset'].isin(['fm_only_pos', 'fm_ilc_pos'])),
    'vascular_subtype'
] = 'Capillary 3'

# Step 2: Change 'Arterial' and 'Artery' → 'Arteriole'
adata.obs['vascular_subtype'] = adata.obs['vascular_subtype'].replace({
    'Arterial': 'Arteriole',
    'Artery': 'Arteriole'
})

# Step 3: Change 'Venous' → 'Venule'
adata.obs['vascular_subtype'] = adata.obs['vascular_subtype'].replace({
    'Venous': 'Venule'
})


# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='dataset',
    groups_col='vascular_subtype',
    layer='MAGIC_imputed_data',
    mode='mean',
    min_cells=10,
    min_counts=100
)
pdata
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# OR convert all non-numeric columns to string (safely)
for col in pdata.obs.columns:
    if not pd.api.types.is_numeric_dtype(pdata.obs[col]):
        pdata.obs[col] = pdata.obs[col].astype(str)

# Save metadata
pdata.obs[['vascular_subtype', 'dataset']].to_csv("output/pseudobulk_metadata.csv")

# Save PCA embeddings
pd.DataFrame(
    pdata.obsm['X_pca'],
    index=pdata.obs_names,
    columns=[f'PC{i+1}' for i in range(pdata.obsm['X_pca'].shape[1])]
).to_csv("output/pseudobulk_pca_coordinates.csv")

sc.pl.pca(pdata, color=['dataset', 'vascular_subtype'], ncols=1, size=300, components = '1,2')
sc.pl.pca_variance_ratio(pdata)