stackVlnPlot <- function(object,
                         features,
                         cluster_name = 'default',
                         psize = 0,
                         n_break = 2,
                         color_pal = NULL,
                         ...) {
  #' Makes a stacked violin plot
  #'
  #' @param object Seurat object
  #' @param features Character vector of features to plot (can be named by cluster)
  #' @param cluster_name Metadata column to use for clustering (default = 'default', uses active ident)
  #' @param psize Dot size on violins (default = 0)
  #' @param n_break Number of y-axis breaks (default = 2)
  #' @param color_pal Optional named vector of colors for fill (names must match cluster levels)
  #' @return A stacked ggplot object
  #' @export

  library(Seurat)
  library(ggplot2)
  library(patchwork)

  # Set identity if a different cluster column is specified
  if (cluster_name != 'default') {
    object <- SetIdent(object, value = cluster_name)
  }

  # Get cluster levels
  cluster_levels <- levels(Idents(object))

  # Validate color palette
  if (!is.null(color_pal)) {
    if (!all(cluster_levels %in% names(color_pal))) {
      stop("All cluster levels must be named in color_pal.")
    }
  }

  # Generate violin plots
  plot_list <- lapply(seq_along(features), function(i) {
    show_x <- (i == length(features))  # Only show x-axis text on the last plot

    p <- VlnPlot(object, features = features[i], pt.size = psize, ...) +
      scale_y_continuous(n.breaks = n_break) +
      theme_minimal(base_family = "Arial") +
      theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold", color = "black", angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", color = "black", angle = if (show_x) 90 else 0),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks = element_line(color = "black")
      ) +
      ylab(features[i])

    if (!show_x) {
      p <- p + theme(axis.text.x = element_blank())
    }

    if (!is.null(color_pal)) {
      p <- p + scale_fill_manual(values = color_pal)
    }

    return(p)
  })

  return(wrap_plots(plot_list, ncol = 1))
}
