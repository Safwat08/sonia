# install.R

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ------------------------------
# Install CRAN packages with versions
# ------------------------------
packages <- list(
  Seurat     = "5.2.0",
  dplyr      = "1.1.4",
  tidyverse  = "2.0.0",
  ggplot2    = "3.5.1",
  stringr    = "1.5.1",
  cowplot    = "1.1.3",
  patchwork  = "1.3.0",
  ggvenn     = "0.1.10",
  ggrepel    = "0.9.6",
  Matrix     = "1.7.1",
  scales     = "1.3.0"
)

for (pkg in names(packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, version = packages[[pkg]])
  }
}

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

devtools::install_github("enblacar/SCpubr", ref = "v2.0.0-dev-stable")
devtools::install_github("jinworks/CellChat", ref = "v2.1.2")
