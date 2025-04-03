# Custom FeaturePlot function for multiple features with grid layout control
customFeaturePlot <- function(object,
                              features,
                              reduction = "umap",
                              metadata = NULL,
                              add_border = TRUE,
                              nrow = NULL,
                              ncol = NULL,
                              pt_size = 1) {

  # Only use features present in data
  features <- intersect(features, c(rownames(object), names(object@meta.data)))

  # List to store individual plots
  plot_list <- list()

  for (feature in features) {
    # Fetch reduction coordinates
    reduction_coordinates <- Embeddings(object, reduction = reduction)

    # Fetch expression/meta data
    data_to_plot <- FetchData(object, vars = c(feature, metadata))

    # Combine data and remove NA
    plot_data <- cbind(reduction_coordinates, data_to_plot)
    plot_data <- plot_data[!is.na(plot_data[, feature]), ]

    # Sort by expression (lower first)
    plot_data <- plot_data[order(plot_data[, feature], decreasing = FALSE), ]

    cols <- colnames(plot_data)

    if(any(grepl('-', c(cols)))) {

      print(paste0(feature, 'has - in name will be replaced with .'))

      dash_ind <- which(grepl('-', c(cols)) == T)

      new_col <- gsub("-", "_", cols[dash_ind])

      cols[dash_ind] <- new_col

      colnames(plot_data) <- cols

      feature <- new_col

      print(paste0('New feature is ', feature))

    }

    # Create the feature plot
    if (add_border) {
      feature_plot <- ggplot(plot_data, aes_string(x = colnames(reduction_coordinates)[1],
                                                   y = colnames(reduction_coordinates)[2])) +
        geom_point(color = "black", size = pt_size*1.5, alpha = 1) +
        geom_point(aes_string(color = feature), size = pt_size)
    } else {
      feature_plot <- ggplot(plot_data, aes_string(x = colnames(reduction_coordinates)[1],
                                                   y = colnames(reduction_coordinates)[2])) +
        geom_point(aes_string(color = feature), size = pt_size)
    }

    feature_plot <- feature_plot +
      scale_color_viridis_c(option = "turbo") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            plot.title = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(face = "bold", size = 12, family = "Arial"),
            axis.ticks = element_blank(),
            axis.line = element_line(),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12)) +
      labs(color = feature, x = "DIM 1", y = "DIM 2")

    plot_list[[feature]] <- feature_plot
  }

  # Use specified layout or guess one
  if (!is.null(nrow) || !is.null(ncol)) {
    grid_plot <- wrap_plots(plot_list, nrow = nrow, ncol = ncol)
  } else {
    n_features <- length(plot_list)
    if (n_features == 1) {
      grid_plot <- plot_list[[1]]
    } else if (n_features <= 4) {
      grid_plot <- wrap_plots(plot_list, ncol = 2, nrow = 2)
    } else if (n_features <= 9) {
      grid_plot <- wrap_plots(plot_list, ncol = 3, nrow = 3)
    } else {
      grid_plot <- wrap_plots(plot_list, ncol = 4, nrow = ceiling(n_features / 4))
    }
  }

  return(grid_plot)
}
