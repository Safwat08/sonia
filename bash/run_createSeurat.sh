#!/bin/bash

# createSeurat

# Define the list properly as an array
list=(
	output/bm_fresh_postqc_analyzed
)
#list=(
 #   output/bm_fresh_postqc_analyzed_subset_analyzed
  #  output/bm_merge_analyzed_subset_analyzed
   # output/fm_bm_merge_analyzed
   # output/fm_merge_analyzed_subset_analyzed
#)

# Iterate over the list
for item in "${list[@]}"; do
    Rscript R/run_createSeurat.R "$item"
done
