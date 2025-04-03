library(Seurat)
library(dplyr)

# Define a function that updates the markers data with mean and median expression for each cluster and gene
getMeanExp <- function(object, markers, assay = 'magic', layer = 'data', cluster_col = 'annotation_1') {

  # Get the expression data for the specified assay and layer
  data <- GetAssayData(object, assay = assay, layer = layer)

  # Loop over each unique cluster in the metadata
  for (cluster in unique(object@meta.data[[cluster_col]])) {

    # Get cells for the current cluster
    cells_cluster <- object@meta.data %>%
      filter(!!sym(cluster_col) == cluster) %>%
      rownames(.)

    # Get cells for the rest of the data (non-cluster cells)
    cells_rest <- object@meta.data %>%
      filter(!!sym(cluster_col) != cluster) %>%
      rownames(.)

    # Subset the expression data for cluster and rest of the cells
    data_cluster <- data[, cells_cluster, drop = FALSE]
    data_rest <- data[, cells_rest, drop = FALSE]

    # Loop over each gene in the markers table
    for (gene in markers$gene) {
      # Make sure the gene exists in the data
      if (gene %in% rownames(data)) {

        # Calculate mean and median expression for the current gene in the cluster
        mean_gene_cluster <- mean(data_cluster[gene, , drop = FALSE], na.rm = TRUE)
        median_gene_cluster <- median(data_cluster[gene, , drop = FALSE], na.rm = TRUE)

        # Calculate mean and median expression for the current gene in the rest of the cells
        mean_gene_rest <- mean(data_rest[gene, , drop = FALSE], na.rm = TRUE)
        median_gene_rest <- median(data_rest[gene, , drop = FALSE], na.rm = TRUE)

        # Update the markers data frame with the calculated values for the current cluster and gene
        markers[markers$gene == gene & markers$cluster == cluster, 'mean_exp_cluster'] <- mean_gene_cluster
        markers[markers$gene == gene & markers$cluster == cluster, 'median_exp_cluster'] <- median_gene_cluster
        markers[markers$gene == gene & markers$cluster == cluster, 'mean_exp_rest'] <- mean_gene_rest
        markers[markers$gene == gene & markers$cluster == cluster, 'median_exp_rest'] <- median_gene_rest
      }
    }
  }

  # Return the updated markers data frame
  return(markers)
}
