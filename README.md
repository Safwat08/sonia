# SONIA  
**Single nuclei Organotypic eNdothelial Integrated Analysis**

SONIA is reproducible analysis pipeline for studying endothelial heterogeneity using snRNA-seq data.
**→ Publication Link: _TBD_**

---

## 📁 Directory Setup

Before running the pipeline, please create the following folders:

rawdata/
output/
figures/

## 📦 Input Data (CellRanger outputs)

Place the **CellRanger (v8.0.1)** outputs into the `rawdata/` directory under the following folder names:

| Sample Description                                             | Folder Name            |
|----------------------------------------------------------------|------------------------|
| fMV-Only GFP+ Fraction                                         | `fm_only_pos`          |
| fMV+ILC GFP+ Fraction                                          | `fm_ilc_pos`           |
| bMV-Only GFP+ Fraction                                         | `bm_only`              |
| bMV+ILC GFP+ Fraction                                          | `bm_ilc`               |
| fMV+ILC GFP– Fraction (dual human-mouse aligned)               | `fm_ilc_neg_mixed`     |
| fMV+ILC GFP– Fraction (human-only aligned)                     | `fm_ilc_neg_human`     |

---

## 🚀 Reproduce Analysis Objects

To generate the processed Seurat objects from CellRanger outputs, run:

```bash
bash run_all.sh
```
	run_all.sh is a bash script that sequentially calls Python and R scripts to process and analyze the data.

Important Notes:

** Between Step 6 and 7, please rename PCs_varm.csv to X_pca_varm.csv

** The fresh rat fat microvessel data was analyzed using a prior pipeline. Please generate this object manually (as described in the publication) and place it in the output/ folder before running analysis_workflow.rmd.


Recreate visualization and cellchat analysis - analysis_workflow.rmd
	analysis_workflow.rmd uses analyzed seurat objects to recreate visualizations and run cellchat analysis

Please consider citing the above publication if using any of the available codes for your own analysis :)

Disclaimer- Generative AI assited to troubleshoot and clean the data
