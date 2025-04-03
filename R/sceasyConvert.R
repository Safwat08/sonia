# Load necessary libraries
library(sceasy)
library(reticulate)

# Use a specific conda environment
use_condaenv('py_env')

# Specify the directory containing the folders
base_dir <- "output/"  # Change this to your actual base directory

# Get a list of subdirectories (folders) in the base directory
folders <- list.dirs(base_dir, recursive = FALSE)  # Only list top-level directories

# Loop over each folder and apply the sceasy::convertFormat function
for (i in 1:length(folders)) {
  
  h5ad_file <- file.path(folders[i], "adata_palantir.h5ad")  # Define path to .h5ad file
  rds_file <- file.path(folders[i], "seuratobj_palantir.rds") # Define output path
  
  # Check if the file exists before trying to convert
  if (file.exists(h5ad_file)) {
    message(paste("Converting:", h5ad_file, "to", rds_file))
    
    # Convert the file
    sceasy::convertFormat(h5ad_file, from="anndata", to="seurat", outFile=rds_file)
  } else {
    message(paste("File not found:", h5ad_file))
  }
}