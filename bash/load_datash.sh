#!/bin/bash

set -e  # Stop on error

# ======================================
# Step 1. Preprocessing data
# ======================================
source /Users/vasconcelos_lab/miniforge3/etc/profile.d/conda.sh
conda activate py_env
echo "Step 1: Preprocessing data..."
python "python/load_data.py"
conda deactivate
