SONIA (Single nuclei Organotypic eNdothelial Integrated Analysis) 

This repository contains the analysis code for the following publication - TBD

Please create the following directories before analyzing:
rawdata/
output/
figures/

cellranger outputs must be placed under "rawdata" in their own directory
	fMV-Only GFP + Fraction - fm_only_pos
	fMV+ILC GFP + Fraction - fm_ilc_pos
	bMV-Only GFP + Fraction - bm_only
	bMV+ILC GFP + Fraction - bm_ilc
	fMV+ILC GFP - Fraction (Aligned to dual human-mouse genome) - fm_ilc_neg_mixed
	fMV+ILC GFP - Fraction (Aligned to only human genome) - fm_ilc_neg_human

Recreate objects run - run_all.sh 
	run_all.sh is a bash script that calls python and Rscripts to produce analyzed objects from cellranger(v8.0.1) 
	** Important between Step 6 and 7, please change name of PCs_varm.csv to X_pca_varm.csv
	** Please note that fresh rat fat microvessel data was analyzed using a prior pipeline as detailed in the publication. Please make the object following the instructions and then place in the output folder prior to running analysis_workflow.rmd

Recreate visualization and cellchat analysis - analysis_workflow.rmd
	analysis_workflow.rmd uses analyzed seurat objects to recreate visualizations and run cellchat analysis

Please consider citing the above publication if using any of the available codes for your own analysis :)

Disclaimer- Generative AI assited to troubleshoot and clean the data
