library(Seurat)
library(Matrix)

# Read converted anndata object (with force)
object <- readRDS('output/subset_analyzed/seuratobj_palantir.rds')

# Load impuated matrix
imputed_matrix <- t(readMM('output/subset_analyzed/imputed_matrix.mtx'))

# Load log1p matrix
log1p_matrix <- t(readMM('output/subset_analyzed/lop1p_matrix.mtx'))

# Load genes and barcodes
genes <- read.csv('output/subset_analyzed/features.csv')

barcodes <- read.csv('output/subset_analyzed/barcodes.csv')

rownames(imputed_matrix) <- genes[,1]

colnames(imputed_matrix) <- barcodes[,1]

rownames(log1p_matrix) <- genes[,1]

colnames(log1p_matrix) <- barcodes[,1]

# Step 3: Add the imputed matrix as a new assay in the Seurat object
object[["MAGIC"]] <- CreateAssay5Object(data = imputed_matrix)
object[["LOG"]] <- CreateAssay5Object(data = log1p_matrix)

# Step
object@assays[["MAGIC"]]@layers[["data"]]@Dimnames[[1]] <- genes[,1]

object@assays[["LOG"]]@layers[["data"]]@Dimnames[[1]] <- genes[,1]

object@assays[["MAGIC"]]@layers[["data"]]@Dimnames[[2]] <- barcodes[,1]

object@assays[["LOG"]]@layers[["data"]]@Dimnames[[2]] <- barcodes[,1]


saveRDS(object, file = 'output/subset_analyzed/seuratobj_palantir.rds')
