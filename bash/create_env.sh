#!/bin/bash

echo "Setting up python env, source conda.sh prior to running"
# Enable 'conda activate' inside this script
# source /Users/vasconcelos_lab/miniforge3/etc/profile.d/conda.sh

# Recreate py_env
conda create -n py_env -f environment_py.yml -y
conda activate py_env
conda deactivate