colorObject <- function(object, group) {

  #Get number of clusters
  clusters <- unique(object@meta.data[[group]])
  nclusters <- length(clusters)


  color_pallete = do_ColorPalette(colors.use = "steelblue", n = nclusters)

  return(color_pallete)

}

