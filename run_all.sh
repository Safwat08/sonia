#!/bin/bash

set -e  # Stop on error

# ======================================
# Enable 'conda/mamba activate' in scripts
# ======================================
source ~/miniconda3/etc/profile.d/conda.sh

# ======================================
# Step 1. Preprocessing data
# ======================================
echo "Step 1: Preprocessing data..."
mamba activate py_env
python "python/load_data.py"
mamba deactivate

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
mamba activate py_env
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

# Subset fm_ilc_neg_postqc
python "python/subset_data.py" "true" \
  "output/fm_ilc_neg_postqc_analyzed" \
  "0" "1" "4" "5"
mamba deactivate

# ======================================
# Step 5. Re-run find deviant genes (post-subset)
# ======================================
echo "Step 5: Finding deviant genes (post-subset)..."
Rscript "R/run_outputDeviantGenes.R" "output/fm_merge_analyzed_subset/"
Rscript "R/run_outputDeviantGenes.R" "output/bm_merge_analyzed_subset/"
Rscript "R/run_outputDeviantGenes.R" "output/fm_ilc_neg_postqc_analyzed_subset/"

# ======================================
# Step 6. Re-analyze data (post-subset)
# ======================================
echo "Step 6: Re-analyzing data (post-subset)..."
mamba activate py_env
python "python/analyze_data.py" \
  "output/fm_merge_analyzed_subset" \
  "output/bm_merge_analyzed_subset" \
  "output/fm_ilc_neg_postqc_analyzed_subset"
mamba deactivate

# ======================================
# Step 7. Create seurat objects to visualize in R
# ======================================
echo "Step 7: Creating seurat object..."

Rscript R/run_createSeurat.R "output/fm_merge_analyzed_subset"
Rscript R/run_createSeurat.R "output/bm_merge_analyzed_subset"
Rscript R/run_createSeurat.R "output/fm_ilc_neg_postqc_analyzed_subset"

# ======================================
# Step 8. Run module score for EC subtypes
# ======================================
echo "Step 8: Running EC subtype scoring..."

Rscript R/run_ec_subtypeScore.R "output/fm_merge_analyzed_subset"
Rscript R/run_ec_subtypeScore.R "output/bm_merge_analyzed_subset"

# ======================================
# Done
# ======================================
echo "All scripts completed successfully"