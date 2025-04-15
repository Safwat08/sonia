# SONIA  
**Single nuclei Organotypic eNdothelial Integrated Analysis**

SONIA is reproducible analysis pipeline for studying endothelial heterogeneity using snRNA-seq data.
**→ Publication Link: _TBD_**

---
## Directory Setup
---

Before running the pipeline, please create the following folders:

```bash
mkdir -p rawdata output figures
```

---
## Input Raw Data (CellRanger outputs)
---

Place the **CellRanger (v8.0.1)** outputs into the `rawdata/` directory under the following folder names:

| Sample Description                                             | Folder Name            |
|----------------------------------------------------------------|------------------------|
| fMV-Only GFP+ Fraction                                         | `fm_only_pos`          |
| fMV+ILC GFP+ Fraction                                          | `fm_ilc_pos`           |
| bMV-Only GFP+ Fraction                                         | `bm_only`              |
| bMV+ILC GFP+ Fraction                                          | `bm_ilc`               |
| fMV+ILC GFP– Fraction (dual human-mouse aligned)               | `fm_ilc_neg_mixed`     |
| fMV+ILC GFP– Fraction (human-only aligned)                     | `fm_ilc_neg_human`     |

**The fresh rat fat microvessel data was analyzed using a prior pipeline. Please generate this object (as described in the publication) and place it in the 'output/' folder before running analysis_workflow.rmd.**

---
## Reproduce Analysis Objects
---

Step 1. To generate the processed Seurat objects from CellRanger outputs using Scanpy, run:

```bash
bash run_all.sh
```
Note
* run_all.sh recreates environments for R and Python using conda with the install.R and environment.yml files respectively

Step 2. To recreate visualization and cellchat analysis use - analysis_workflow.rmd

```bash
conda activate r_env
Rscript rmarkdown::render('analysis_workflow.Rmd')
```

If you use any part of this pipeline in your own work, please consider citing the corresponding publication.

Note: Some steps in this pipeline were iteratively refined using generative AI assistance (e.g., ChatGPT) to troubleshoot issues and improve reproducibility.