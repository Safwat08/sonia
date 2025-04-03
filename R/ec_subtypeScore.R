# Load necessary libraries
library(Seurat)
library(dplyr)  # For handling sparse matrices
library(ggplot2)   # For reading CSV files
library(stringr)

source('R/addScore.R')
source('R/markerDotPlot.R')

ec_subtypeScore <- function(object, meta, ec_df, name, return = F, plot = T) {

  # Get all the labels
  labels <- unique(ec_df$label)

  score_names <- c()

  for (i in 1:length(labels)) {

    genes <- ec_df %>%
      filter(label == labels[i]) %>%
      pull(genes)

    # Use addScore function to obtain scores using Seurat's method
    object <- addScore(object, genes, score_name = labels[i], score_type = 'seurat')

    # Get the score names
    score_names[i] <- paste0(labels[i], "_score_seurat")

  }

  # organ score

  # Split into vascular subtype scores and organ-wide scores
  organ_score_ind = grep("organ", score_names)

  score_subtypeonly = score_names[-organ_score_ind]

  score_organonly = score_names[organ_score_ind]

  # Subtype score
  dotPlot_list <- betterDotPlot(object, features = score_subtypeonly, metadata = meta, cluster_group = T, cluster_feature = T)


  dotPlot <- dotPlot_list[[1]]

  features <-  dotPlot_list[[2]]

  tile_df <- tibble(
    label = str_to_lower(features),               # First column: feature labels
    x = 1:length(features),    # Second column: sequence from 1 to length of features
    y = 1                    # Third column: all values set to 1
  )

  tile_df <- tile_df %>%
    mutate(fill = case_when(
      str_detect(label, "large") & str_detect(label, "artery") ~ "Large Artery",
      str_detect(label, "artery") & !str_detect(label, "large") ~ "Artery",
      str_detect(label, "arterial") | str_detect(label, "arteriole") ~ "Arteriole",
      str_detect(label, "vein") ~ "Vein",
      str_detect(label, "venous") | str_detect(label, "venule") ~ "Venule",
      str_detect(label, "proliferating") ~ "Proliferating",
      str_detect(label, "capillary") & !str_detect(label, "arterial") &
        !str_detect(label, "venous") & !str_detect(label, "arteriole") &
        !str_detect(label, "venule") ~ "Capillary",
      str_detect(label, "lymphatic") | str_detect(label, "lypmhatic") ~ "Lymphatic",
      str_detect(label, "angiogenic") ~ "Angiogenic",
      str_detect(label, "choroidplexus") | str_detect(label, "glomeruli") |
        str_detect(label, "colon_...7") ~ "Capillary",
      str_detect(label, "interferon") ~ "Immunogenic",
    ))

  # Custom color palette for each label
  color_pallete <- c("Large Artery" = "#AD002A",
                     "Artery" = '#EE0000',
                     "Arteriole" = '#EFC000',
                     "Vein" = '#631879',
                     "Venule" = '#AF46B4',
                     "Capillary" = '#008280',
                     "Proliferating" = '#BA6338',
                     "Lymphatic" = '#F7F7F7',
                     'Angiogenic' = '#5050FF',
                     "Immunogenic" = 'yellow'
  )

  tile_df$fill <- factor(tile_df$fill, levels = names(color_pallete))

  # Plot the geom_tile with custom colors, including duplicates
  # Create the tile plot
  tilePlot <- ggplot(tile_df, aes(x = x, y = y, fill = fill)) +
    geom_tile(color = "black", width = 0.9, height = 0.9, size = 0.9) +
    scale_fill_manual(values = color_pallete) +
    theme_minimal() +
    labs(x = "", y = "", fill = "Vascular Subtype") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank()) +
    coord_fixed(ratio = 1)


  combinedPlot <- plot_grid(dotPlot, NULL, tilePlot, ncol = 1, rel_heights = c(1, -0.2, 0.7))

  fileName_obj <- paste0("figures/score_dotplot_vascular_subtype_", name)

  fileName_obj_ext <- paste0(fileName_obj, ".png")

  ggsave(combinedPlot,
         filename = fileName_obj_ext,
         width = 17,
         height = 7)


}
