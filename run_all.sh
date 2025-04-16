#!/bin/bash

set -e  # Stop on error

# ======================================
# Step 1. Preprocessing data
# ======================================
conda activate py_env
echo "Step 1: Preprocessing data..."
python "python/load_data.py"
conda deactivate

# ======================================
# Step 2. Find deviant genes (pre-analysis)
# ======================================

echo "Step 2: Finding deviant genes (initial)..."
Rscript "R/run_outputDeviantGenes.R" "output/bm_merge/"
Rscript "R/run_outputDeviantGenes.R" "output/fm_merge/"
Rscript "R/run_outputDeviantGenes.R" "output/fm_ilc_neg_postqc/"

# ======================================
# Step 3. Analyze data
# ======================================
echo "Step 3: Analyzing data..."
conda activate py_env
python "python/analyze_data.py" \
  "output/bm_merge" \
  "output/fm_merge" \
  "output/fm_ilc_neg_postqc"

# ======================================
# Step 4. Subset data
# ======================================
echo "Step 4: Subsetting data..."

# Subset fm_merge
python "python/subset_data.py" "true" \
  "output/fm_merge_analyzed" \
  "1" "2" "0" "9" "8" "14"

# Subset bm_merge
python "python/subset_data.py" "true" \
  "output/bm_merge_analyzed" \
  "6" "2" "7" "3" "11" "5" "1" "0"

conda deactivate

# ======================================
# Step 5. Re-run find deviant genes (post-subset)
# ======================================
echo "Step 5: Finding deviant genes (post-subset)..."
Rscript "R/run_outputDeviantGenes.R" "output/fm_merge_analyzed_subset/"
Rscript "R/run_outputDeviantGenes.R" "output/bm_merge_analyzed_subset/"

# ======================================
# Step 6. Re-analyze data (post-subset)
# ======================================
echo "Step 6: Re-analyzing data (post-subset)..."
conda activate py_env
python "python/analyze_data.py" \
  "output/fm_merge_analyzed_subset" \
  "output/bm_merge_analyzed_subset" \
conda deactivate


##!!!! Note you have to go in manually to change PC_varm.csv to X_pca_varm.csv!!!
find output/ -type f -name "PCs_varm.csv" -execdir mv {} X_pca_varm.csv \; 


# ======================================
# Step 7. Create seurat objects to visualize in R
# ======================================
echo "Step 7: Creating seurat object..."
Rscript R/run_createSeurat.R "output/fm_merge_analyzed_subset_analyzed"
Rscript R/run_createSeurat.R "output/bm_merge_analyzed_subset_analyzed"
Rscript R/run_createSeurat.R "output/fm_ilc_neg_postqc_analyzed"

# ======================================
# Step 8. Run module score for EC subtypes
# ======================================
echo "Step 8: Running EC subtype scoring..."

Rscript R/run_ec_subtypeScore.R "output/fm_merge_analyzed_subset_analyzed"
Rscript R/run_ec_subtypeScore.R "output/bm_merge_analyzed_subset_analyzed"

# ======================================
# Done
# ======================================
echo "All scripts completed successfully"