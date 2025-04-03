# Load necessary libraries
library(Seurat)
library(Matrix)  # For handling sparse matrices
library(readr)   # For reading CSV files
library(stringr)

createSeurat <- function(output_dir) {

  # Load the metadata and count matrices
  gene_meta <- read.csv(file.path(output_dir, "features.csv"), header = T)
  cell_meta <- read.csv(file.path(output_dir, "barcodes.csv"), header = T, row.names = 1)

  # Load all matrix files
  matrix_files <- list.files(output_dir, pattern = '.mtx')

  # Get matrix names from file name
  matrix_names <- sapply(matrix_files, function(x) {
    str_split_fixed(x, pattern = '_matrix.mtx', n = 2)[1]
  })

  # Initiate an empty list to store all matrix
  matrix_list <- list()

  # Loop over matrix files and store in matrix list
  for (i in 1:length(matrix_files)) {

    # Load matrix and transpose to fit seurat object
    matrix <- t(as(as.matrix(readMM(file.path(output_dir, matrix_files[i]))), "dgCMatrix"))

    # Gene names
    rownames(matrix) <- gene_meta['gene_names']$gene_names # Assuming the first column of features.csv contains gene names

    # Cell barcodes
    colnames(matrix) <- rownames(cell_meta)

    # Put matrix in list
    matrix_list[[i]] <- matrix

  }

  # Name matrix list by matrix names
  names(matrix_list) <- matrix_names

  # Create Seurat object using rawcounts
  seurat_obj <- CreateSeuratObject(counts = matrix_list[['rawcounts']], data = matrix_list[['rawcounts']])

  # If more than 1 matrix present
  if (length(matrix_list) > 1) {

    matrix_list[['rawcounts']] <- NULL

    # Assuming first matrix is rawcounts
    for (i in 1:length(matrix_list)) {

      # Add other matrix to assays
      seurat_obj[[matrix_names[i]]] <- CreateAssay5Object(counts = matrix_list[[i]], data = matrix_list[[i]])

    }

  }

  # Add cell metadata from barcodes.csv
  seurat_obj <- AddMetaData(seurat_obj, metadata = cell_meta)

  # Add gene metadata to the RNA assay
  # List all assays in the Seurat object
  #assay_names <- names(seurat_obj@assays)

  #for (assay in assay_names) {

  # gene_meta <- gene_meta[, !colnames(gene_meta) %in% "gene_meta.1"]

  # Check if the assay contains gene data (not all assays may be relevant)
  #seurat_obj@assays[[assay]]@meta.data <- gene_meta
  #}

  # Function to load embeddings
  load_embeddings <- function(file_path) {
    # Load the matrix
    data <- read_csv(file_path)

    # Set the first column as rownames and remove it from the data
    rownames(data) <- data[[1]]  # Set the first column as rownames (cell names)
    data <- data[,-1]  # Remove the first column which is now the rownames

    # Convert to a numeric matrix directly
    numeric_matrix <- as.matrix(data)

    colnames(numeric_matrix) <- 1:ncol(numeric_matrix)

    return(numeric_matrix)
  }


  # Load matrix files
  obsm_files <- list.files(output_dir, pattern = '_obsm.csv')

  # IF obsm files present
  if(length(obsm_files) > 0) {

    obsm_names <- sapply(obsm_files, function(x) {
      str_split_fixed(x, pattern = '_obsm.csv', n = 2)[1]
    })

    obsm_list <- list()

    for (i in 1:length(obsm_files)) {

      obsm <- load_embeddings(file.path(output_dir, obsm_files[i]))

      rownames(obsm) <- rownames(cell_meta)

      obsm_list[[i]] <- obsm

    }

    names(obsm_list) <- obsm_names

    for (i in 1:length(obsm_list)) {

      KEY = paste0(toupper(obsm_names[[i]]), "_")

      seurat_obj[[obsm_names[i]]] <- CreateDimReducObject(embeddings = obsm_list[[i]], key = KEY)

    }
  }


  # Load matrix files
  varm_files <- list.files(output_dir, pattern = '_varm.csv')

  if(length(varm_files) > 0) {

    varm_names <- sapply(varm_files, function(x) {
      str_split_fixed(x, pattern = '_varm.csv', n = 2)[1]
    })

    varm_list <- list()

    for (i in 1:length(varm_files)) {

      varm <- load_embeddings(file.path(output_dir, varm_files[i]))

      KEY = paste0(gsub("_", "", toupper(varm_names[[1]])), "_")

      rownames(varm) <- gene_meta$gene_names
      colnames(varm) <- paste0(KEY, colnames(varm))

      varm_list[[i]] <- varm

    }

    names(varm_list) <- varm_names

    for (i in 1:length(varm_list)) {

      seurat_obj[[varm_names[i]]]@feature.loadings <- varm_list[[i]]
    }
  }

  # View Seurat object
  seurat_obj

  save_file <- file.path(output_dir, "seurat_obj.rds")

  saveRDS(seurat_obj, file = save_file)

}

