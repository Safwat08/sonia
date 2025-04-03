



#Custom GeyserPlot function for multiple features with position_jitter and seed





customGeyserPlot <- function(object,
                             plot_feature = TRUE,
                             features = NULL,
                             metadata = NULL,
                             add_border = TRUE,
                             add_box = TRUE,
                             order_feature = TRUE) {

  if (plot_feature) {

      # Fetch expression data for the given features
    data_to_plot <- FetchData(object, vars = c(features, metadata))

    # Melt the data into long format (for ggplot)
    plot_data_long <- reshape2::melt(data_to_plot, variable.name = "Feature", value.name = "Expression")

    # Order features if specified
    if (order_feature) {
      # Calculate mean expression for each feature
      mean_expression <- aggregate(Expression ~ Feature, data = plot_data_long, FUN = mean)
      ordered_features <- mean_expression$Feature[order(mean_expression$Expression, decreasing = TRUE)]
      plot_data_long$Feature <- factor(plot_data_long$Feature, levels = ordered_features)
    } else {
      # Convert 'Feature' to a factor for proper ordering
      plot_data_long$Feature <- as.factor(plot_data_long$Feature)
    }

    if (!is.null(metadata)) {
      if (add_border) {

        geyser_plot <- ggplot(plot_data_long, aes(x = Feature, y = Expression, color = !!sym(metadata))) +
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2.5, stroke = 0.5, shape = 21, color = "black", aes(fill = !!sym(metadata))) +  # First layer with black outline
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2)  # Second layer with main colors +

      } else {

        geyser_plot <- ggplot(plot_data_long, aes(x = Feature, y = Expression, color = !!sym(metadata))) +  # First layer with black outline
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2)  # Second layer with main colors
      }
    } else {
      # Add jittered points with or without border
      if (add_border) {
        geyser_plot <- ggplot(plot_data_long, aes(x = Feature, y = Expression)) +
          geom_point(position = position_jitter(width = 0.2, height = 0, seed = 0),
                     shape = 21, color = "black", fill = "white", size = 2.5) +
          geom_point(aes(color = Expression),
                     position = position_jitter(width = 0.2, height = 0, seed = 0), size = 2)
      } else {
        geyser_plot <- ggplot(plot_data_long, aes(x = Feature, y = Expression)) +
          geom_point(aes(color = Expression),
                     position = position_jitter(width = 0.2, height = 0, seed = 0), size = 2)
      }

    }

    # Add color scale based on usage of metadata
    if (is.null(metadata)) {
      geyser_plot <- geyser_plot + scale_color_viridis_c(option = "turbo", name = "Expression")
    } else {
      geyser_plot <- geyser_plot +
        scale_color_manual(values = do_ColorPalette('steelblue', n = length(unique(plot_data_long[, 1]))))
    }

  } else {

    # Fetch expression data for the given features
    data_to_plot <- FetchData(object, vars = c(features, metadata))

    # Melt the data into long format (for ggplot)
    plot_data_long <- reshape2::melt(data_to_plot, variable.name = "Feature", value.name = "Expression")

    # Order features if specified
    if (order_feature) {
      # Calculate mean expression for each feature
      mean_expression <- aggregate(Expression ~ Feature, data = plot_data_long, FUN = mean)
      ordered_features <- mean_expression$Feature[order(mean_expression$Expression, decreasing = TRUE)]
      plot_data_long$Feature <- factor(plot_data_long$Feature, levels = ordered_features)
    } else {
      # Convert 'Feature' to a factor for proper ordering
      plot_data_long$Feature <- as.factor(plot_data_long$Feature)
    }

    if (!is.null(feature)) {
      if (add_border) {

        geyser_plot <- ggplot(plot_data_long, aes(x = !!sym(metadata), y = Expression, color = Feature)) +
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2.5, stroke = 0.5, shape = 21, color = "black", aes(fill = Feature)) +  # First layer with black outline
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2)  # Second layer with main colors +

      } else {

        geyser_plot <- ggplot(plot_data_long, aes(x = !!sym(metadata), y = Expression, color = Feature)) + # First layer with black outline
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 42),
                     size = 2)  # Second layer with main colors
      }
    } else {
      # Add jittered points with or without border
      if (add_border) {
        geyser_plot <- ggplot(plot_data_long, aes(x = !!sym(metadata), y = Expression)) +
          geom_point(position = position_jitter(width = 0.2, height = 0, seed = 0),
                     shape = 21, color = "black", fill = "white", size = 2.5) +
          geom_point(aes(color = Expression),
                     position = position_jitter(width = 0.2, height = 0, seed = 0), size = 2)
      } else {
        geyser_plot <- ggplot(plot_data_long, aes(x = !!sym(metadata), y = Expression)) +
          geom_point(aes(color = Expression),
                     position = position_jitter(width = 0.2, height = 0, seed = 0), size = 2)
      }

    }
  }
    # Add color scale based on usage of metadata
    if (is.null(metadata)) {
      geyser_plot <- geyser_plot + scale_color_viridis_c(option = "turbo", name = "Score")
    } else {
      geyser_plot <- geyser_plot +
        scale_color_manual(values = do_ColorPalette('steelblue', n = length(unique(plot_data_long[, 1]))))
    }

  # Customize theme
  geyser_plot <- geyser_plot +
    theme_minimal() +
    theme(panel.grid = element_blank(),  # Remove grid lines
          axis.line = element_line(color = "black"),  # Add x and y axis lines
          axis.ticks = element_line(color = "black"),  # Add ticks for axes
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = 'bold', family = 'Arial'),  # Angle the feature names for better readability
          axis.text.y = element_text(size = 12, face = 'bold', family = 'Arial'),  # Angle the feature names for better readability
          axis.title.y = element_text(size = 12, face = 'bold', family = 'Arial'),
          legend.title = element_text(size = 12, face = 'bold', family = 'Arial'),
          legend.text = element_text(size = 12, face = 'bold', family = 'Arial')) +
    labs(x = "Features", y = "Score", color = if (is.null(metadata)) "Score" else metadata)

  return(geyser_plot)

}

