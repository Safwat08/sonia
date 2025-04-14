# install.R

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ------------------------------
# Install CRAN packages with versions
# ------------------------------
install.packages("Seurat", version = "5.2.0")
install.packages("dplyr", version = "1.1.4")
install.packages("tidyverse", version = "2.0.0")
install.packages("ggplot2", version = "3.5.1")
install.packages("stringr", version = "1.5.1")
install.packages("cowplot", version = "1.1.3")
install.packages("patchwork", version = "1.3.0")
install.packages("ggvenn", version = "0.1.10")
install.packages("ggrepel", version = "0.9.6")
install.packages("Matrix", version = "1.7.1")
install.packages("scales", version = "1.3.0")

# ------------------------------
# Install Bioconductor packages
# ------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("progeny", version = "3.18")
BiocManager::install("decoupleR", version = "3.18")

# ------------------------------
# Install GitHub packages (locked to latest release tags)
# ------------------------------
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("lazappi/SCpubr@v2.0.2")
devtools::install_github("sqjin/CellChat@v2.1.2")