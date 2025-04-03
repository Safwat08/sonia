#!/bin/bash

# Find deviant genes for the merge dataset

# Define a listt
list='output/fm_bm_fresh_merge/'
#list="output/bm_merge_analyzed_subset/ output/fm_merge_analyzed_subset/ output/bm_fresh_postqc_analyzed_subset/"
#list="output/fm_merge_analyzed/_subset/"
#list="output/bm_fresh_postqc/"
#list="output/fm_merge/ output/bm_merge/ output/fm_ilc_neg_postqc/ output/bm_fresh_postqc/"
#list="output/fm_merge_analyzed_subset/ output/bm_merge_analyzed_subset/ output/all_merge_analyzed_subset/"
#list="output/fm_merge_analyzed_subset/ output/bm_merge_analyzed_subset/ output/all_merge_analyzed_subset/"
#list="output/bm_merge_analyzed_subset_analyzed_subset/"

# Iterate over the list
for item in $list; do
    # Use the variable in the Rscript command
    Rscript R/run_outputDeviantGenes.R "$item"
done

