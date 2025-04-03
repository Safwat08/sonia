library(celda)
library(Seurat)

#decontObject
decontObject <- function(object) {

  #'@title decontObject
  #'
  #'@description This function decontaminates counts data using decontx
  #'
  #'@param object Seurat S4.
  #'
  #'
  #'@importFrom Seurat CreateAssay5Object GetAssayData
  #'@importFrom celda decontX
  #'
  #'@return is a seurat S4 object, with decontaminated counts and contamination metadata
  #'
  #'@export
  #'
  #'@examples
  #' # To decontaminate a seurat object
  #' object <- decontObject(object)

  # Decontamination Analysis

  # Get assay data
  counts <- GetAssayData(object, assay = 'RNA', slot = 'counts')

  # Get decontaminated counts
  decont_results <- decontX(counts)
  decont_counts <- decont_results$decontXcounts

  # Make decontaminated assay
  decontx_assay <- CreateAssay5Object(counts = decont_counts)

  # Assay
  object[["decontx"]] <- decontx_assay
  object@meta.data[['contamination']] <- decont_results[["contamination"]]

  return(object)

}
