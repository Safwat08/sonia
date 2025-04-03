library(scry)

outpuDeviantGenes <- function(data,
                              n = 4000) {
  #'@title outputDeviantGenes
  #'
  #'@description This function finds the top 4000 most deviant genes for feature selection
  #'
  #'@param data SparseMatrix, with genes as rownames and columns as cell names
  #'@param n Numeric, top n genes to use. Default is 4000
  #'
  #'@importFrom scry devianceFeatureSelection
  #
  #'@return is vector of genes
  #'
  #'@export
  #'
  #'@examples
  #' # Use top 3000 deviant genes
  #' deviant_genes <- outputDeviantGenes(object, n = 3000)


  deviant_genes <- devianceFeatureSelection(data)
  deviant_genes <- sort(deviant_genes, decreasing = T)
  top_4000 <- names(deviant_genes[1:4000])
  return(top_4000)

}
