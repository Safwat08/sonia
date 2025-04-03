library(Seurat)
library(ggtree)
library(cowplot)
library(tidyr)
library(ggplot2)

betterDotPlot <- function(object, features, metadata, cluster_group = TRUE, cluster_feature = TRUE, scaledots = T) {

  #'@title betterDotPlot
  #'
  #'@description This function creates a dotplot for a set of genes with optional hierarchical clustering.
  #'
  #'@param object Seurat S4 object.
  #'@param features Character vector of features (genes) to plot.
  #'@param metadata Character, metadata to group the analysis.
  #'@param cluster_group Logical, whether to perform hierarchical clustering on the groups.
  #'@param cluster_feature Logical, whether to perform hierarchical clustering on the features (genes).
  #'
  #'@import ggplot2
  #'@importFrom Seurat SetIdent DotPlot
  #'@importFrom celda decontX
  #'@importFrom ggtree ggtree
  #'
  #'@return A list containing the dotplot and ordered features.
  #'
  #'@examples
  #'plot <- betterDotPlot(object, features, "metadata_column")

  object <- SetIdent(object, value = metadata)

  # Generate the initial dot plot
  p <- DotPlot(object, features = features)

  # Process data for clustering
  if (scaledots) {
    data <- p$data %>%
      filter(!is.na(avg.exp.scaled)) %>%
      mutate(weighted.avg.exp.scaled = avg.exp.scaled * pct.exp / 100) %>%
      select(-pct.exp, -avg.exp, -avg.exp.scaled) %>%
      pivot_wider(names_from = features.plot, values_from = weighted.avg.exp.scaled) %>%
      as.data.frame()

  } else {
    data <- p$data %>%
      filter(!is.na(avg.exp.scaled)) %>%
      mutate(weighted.avg.exp.scaled = avg.exp * pct.exp / 100) %>%
      select(-pct.exp, -avg.exp, -avg.exp.scaled) %>%
      pivot_wider(names_from = features.plot, values_from = weighted.avg.exp.scaled) %>%
      as.data.frame()
  }


  if (cluster_group | cluster_feature) {

    # Set row names and convert to matrix for clustering
    rownames(data) <- as.vector(data[,1])
    data <- data[,-1]
    mat <- as.matrix(data)

    # Perform hierarchical clustering
    clust <- hclust(dist(mat))  # Cluster for groups (samples)
    clust_features <- hclust(dist(t(mat)))  # Cluster for features (genes)

    # Order the groups based on clustering
    my_levels <- sub("[.]", " ", clust$labels[clust$order])

    # Create dendrogram plots using ggtree
    ggtree_plot <- ggtree(clust)
    ggtree_features <- ggtree(clust_features)

  }


  # Re-level the object based on clustering if required
  if (cluster_group) {
    object@active.ident <- factor(x = object@active.ident, levels = my_levels)
  }

  # Order features based on clustering if required
  if (cluster_feature) {
    features <- clust_features$labels[clust_features$order]
  }

  # Generate the final dot plot with ordered features



  dPlot <- DotPlot(object, features = features, dot.scale = 5, scale = scaledots)

  if (scaledots) {
    dPlot <- dPlot + scale_color_gradient2(low = 'blue3', high = 'red3', mid = 'white', midpoint = 0,
                                   guide = guide_colorbar(
                                     frame.colour = "black",
                                     ticks.colour = "black",
                                     title = 'Avg. Exp'
                                   ))

  } else {
    dPlot <- dPlot + scale_color_viridis_c(option = "turbo",
                                           guide = guide_colorbar(
                                             frame.colour = "black",
                                             ticks.colour = "black",
                                             title = 'Avg. Exp'
                                           ))

  }
    dPlot <- dPlot +
    geom_point(aes(size = pct.exp),  # Ensure size is linked to pct.exp as in DotPlot
               shape = 21,           # Shape 21 allows both fill and stroke
               color = "black",      # Black stroke (outline)
               stroke = 0.5) +       # Stroke proportional to size
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8, face = "bold", colour = "black", vjust = 1, hjust = 1, angle = 90),
      axis.text.y = element_text(size = 8, face = "bold", colour = "black"),
      legend.text = element_text(size = 8, face = "bold", colour = "black"),
      legend.title = element_text(size = 8, face = "bold", colour = "black"),
      legend.key.width = unit(0.5, 'cm'),
      legend.key.height = unit(0.3, 'cm'),
      legend.position = 'right'
    )

  # Combine dendrogram and dot plot if clustering groups
  if (cluster_group) {
    plot <- plot_grid(ggtree_plot, dPlot, nrow = 1, rel_widths = c(0.1, 1), align = 'h')
  } else {
    plot <- dPlot
  }

  # Return the combined plot and ordered features
  return(list(plot = plot, features = features))
}

#markerDotPlot
markerDotPlot <- function(object,
                          metadata,
                          topn = 5,
                          cluster = T) {
  #'@title markerDotPlot
  #'
  #'@description This function creates a dotplot with top markers for each clusters
  #'
  #'@param object Seurat S4.
  #'@param metadata Character, metadata to conduct analysis on
  #'@param topn Numeric, number of genes to use in dotplot
  #'@param cluster Logical, whetehr to conduct hierarchical clustering or not
  #'
  #'@importFrom Seurat SetIdent FindAllMarkers
  #'@importFrom celda decontX
  #'
  #'@return is a ggplot2 object
  #'
  #'@export
  #'
  #'@examples
  #' # To decontaminate a seurat object
  #' object <- decontObject(object)

  # Ident
  object <- SetIdent(object, value = metadata)

  # Marker analysis
  markers <- FindAllMarkers(object)

  genes <- markers %>%
    filter(pct.1 > 0.1 &
             avg_log2FC > 0.25 &
             p_val_adj < 0.05) %>%
    mutate(feature = rownames(.),
           pct.diff = pct.1 - pct.2,
           norm.pct.diff = pct.1/pct.1 - pct.2/pct.1,
           neglogpval = -log(p_val_adj),
           neglogpval = ifelse(is.infinite(neglogpval),
                               max(-log(.$p_val_adj)[!is.infinite(-log(.$p_val_adj))]),
                               neglogpval),
           log2FC.neglogpval = avg_log2FC*neglogpval) %>%
    group_by(cluster) %>%      # Group by column A
    slice_max(order_by = log2FC.neglogpval, n = topn) %>%  # Replace 'some_column' with the column you want to rank by
    ungroup() %>%
    pull(gene)

  p <- betterDotPlot(object,
                     features = unique(genes),
                     metadata = metadata,
                     cluster_group = cluster,
                     cluster_feature = F)

  return(p)

}
