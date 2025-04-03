## --- VlnDotPlot ----------------------------------------------------------- ##

VlnDotPlot <- function(object,
                       feature,
                       group,
                       fill_name = 'Average Score',
                       y_name = 'Absolute Score',
                       order = F) {

  data <- FetchData(object, vars = c(feature, group)) %>%
    rename_all(~ c('feature', 'group'))

  data_means <- data %>%
    group_by(group) %>%
    dplyr::summarise(mean_feature = mean(feature))

  data <- data %>%
    left_join(data_means, by = join_by(group == group))

  if (order == T) {

    group_order <- data_means %>%
      arrange(mean_feature) %>%
      pull(group)


  } else {

    group_order <- data_means %>%
      pull(group)

  }

  p <- ggplot(data, aes(x = factor(group, levels = group_order),
                        y = feature,
                        fill = mean_feature)) +
    geom_violin()  +
    scale_fill_viridis(na.value = 'black',
                       option = 'D',
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black",
                                              title = 'Avg. Score')) +
    stat_summary(aes(x = group,
                     y = feature), fun="mean",
                 geom="point", color="red", size = 1) +

    theme_classic() +
    theme( plot.title = element_blank(),
           axis.title.y = element_blank(),
           axis.title.x = element_text(size=8, face = "bold", colour = "black"),
           axis.text.x = element_text(size =8, angle = 90, face = "bold", colour = "black"),
           axis.text.y = element_text(size =8, face = "bold", colour = "black"),
           legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
           legend.text = element_text(size=8, face = "bold", colour = "black"),
           legend.title.align = 0.5,
           legend.key.width = unit(0.5, 'cm'),
           legend.key.height = unit(0.3, 'cm'))+
    labs(fill = fill_name,
         y = y_name) +
    coord_flip()

  return(p)


}
