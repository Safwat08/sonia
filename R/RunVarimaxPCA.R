library(pracma)
library(stats)

## --- RunVarimaxPCA -------------------------------------------------------- ##

RunVarimaxPCA <- function(object,
                          features = 'default',
                          reduction = "pca",
                          ncomp = 'default',
                          norm = 'RNA',
                          rotmat = T,
                          scaleby_sdev = T,
                          reduction.name = 'vpca') {

  #' Runs varimax rotation on PC
  #'
  #' Conducts varimax on the PCA loadings and then adds the new rotated loadings
  #' to a different dimension object in th seurat object
  #'
  #'@param object Seurat object with PCA already done
  #'@param reduction The name of the reduction the pca results is stored under
  #'@param norm Which assay to get scaled data from
  #'@param ncomp Number of principal components to use for rotation
  #'@param rotmat Whether to use rotmat from varimax to get new scores (Default)
  #'or mutliply by the inverse of the rotated loadings
  #'@param scaleby_sdev Whether to scale the pca loadings by the standard deviation
  #'
  #'@importFrom stats varimax
  #'@importFrom pracma pinv
  #'
  #'
  #'@return is a seurat object with varimax pca in a new reduction
  #'
  #'@export

  # Get pca results
  # prcomp_results$x or irlba_results$u %*% diag(irlba_results$d)is stored in cell.embeddings
  # prcomp_results$rotation or irlba_results$v is stored under feature.loadings

  pca_cell_embeddings <- object@reductions[[reduction]]@cell.embeddings
  pca_feature_loadings  <- object@reductions[[reduction]]@feature.loadings
  sdev <- object@reductions[[reduction]]@stdev


  if ('default' %in% features) {
    # Get data for variable features
    features <- VariableFeatures(object)
    data <- FetchData(object, vars = features, slot = 'data') #order is off in feature loadings
    data <- t(data[,rownames(pca_feature_loadings)])
  } else {
    features <- intersect(features, rownames(object))
    data <- FetchData(object, vars = features, slot = 'data') #order is off in feature loadings
  }

  if (ncomp == 'default') {
    ncomp <- ncol(pca_feature_loadings)
  }

  pca_cell_embeddings <- pca_cell_embeddings[,1:ncomp]

  pca_feature_loadings <- pca_feature_loadings[,1:ncomp]

  # Modified from amoeba, 2013, stackoverflow
  # Scale by standard deviation or not
  if (scaleby_sdev == T) {
    sdev <- apply(pca_cell_embeddings, 2, sd)
    rawLoadings <- pca_feature_loadings[,1:ncomp] %*% diag(sdev, ncomp, ncomp)
  } else {
    rawLoadings <- pca_feature_loadings
  }

  if (rotmat == T) {
    # Scale
    rotatedLoadings <- varimax(rawLoadings)$loadings
    scores <- scale(pca_cell_embeddings) %*% varimax(rawLoadings)$rotmat
  } else {
    rotatedLoadings <- varimax(rawLoadings)$loadings
    invLoadings     <-  pinv(rotatedLoadings)
    scores          <- scale(data) %*% invLoadings
  }

  #Convert to matrix (adapted from stackoverflow jay.sf, Dec, 2018)
  rotatedLoadings_matrix <- as.matrix(data.frame(matrix(as.numeric(rotatedLoadings),attributes(rotatedLoadings)$dim,dimnames=attributes(rotatedLoadings)$dimnames)))

  vpc_names <- paste0("VPC_",1:ncomp)

  colnames(rotatedLoadings_matrix) <- vpc_names
  colnames(scores) <- vpc_names

  vpca <- CreateDimReducObject(
    embeddings = scores,
    loadings = rotatedLoadings_matrix,
    assay = norm,
    key = 'VPC_',
  )

  object@reductions[[reduction.name]] <- vpca

  return(object)

}
