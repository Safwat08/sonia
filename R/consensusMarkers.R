pairwiseMarkers <- function(object,
                            metadata,
                            population1,
                            population2 = 'all',
                            pairwise,
                            pairwise_consensus,
                            pval_cutoff = 0.05,
                            pct_cutoff = 0.2,
                            log2FC_cutoff = 1) {

  if (!is.factor(object@meta.data[[metadata]])) {
    # Convert to factor if not already
    object@meta.data[[metadata]] <- as.factor(object@meta.data[[metadata]])
  }

  all_classes <- levels(object@meta.data[[metadata]])

  object <- SetIdent(object, value = metadata)

  if (length(population2) < 2 && any(population2 == 'all')) {
    population2 <- setdiff(all_classes, population1)
  }

  pop1_intersect <- intersect(population1, all_classes)
  pop1_setdiff <- setdiff(population1, all_classes)

  if (length(pop1_setdiff) > 0) {
    warning("Following population 1 not detected and will not be used in DE:", pop1_setdiff)
  }

  if (length(pop1_intersect) < 1) {
    stop("None of the population 1 classes exist in metadata:", metadata)
  }

  pop2_intersect <- intersect(population2, all_classes)
  pop2_setdiff <- setdiff(population2, all_classes)

  if (length(pop2_setdiff) > 0) {
    warning("Following population 2 not detected and will not be used in DE:", pop2_setdiff)
  }

  if (length(pop2_intersect) < 1) {
    stop("None of the population 2 classes exist in metadata:", metadata)
  }

  if (length(pop2_intersect) < 2 && pairwise) {
    warning("Less than two population2 classes detected, simply conducting FindMarkers between", population1, "and", pop2_intersect)
    pairwise <- FALSE
  }

  markers <- c()

  if (pairwise) {
    for (i in seq_along(pop2_intersect)) {
      message("Conducting pairwise DE between:", population1, "and", pop2_intersect[i])

      object_markers <- FindMarkers(object,
                                    ident.1 = pop1_intersect,
                                    ident.2 = pop2_intersect[i])

      object_markers <- object_markers %>%
        filter(p_val_adj < pval_cutoff &
                 pct.1 > pct_cutoff &
                 avg_log2FC > log2FC_cutoff) %>%
        mutate(feature = rownames(.))

      new_markers <- object_markers$feature

      if (length(new_markers) < 1) {
        warning("No markers found between", population1, "and", pop2_intersect[i])
      }

      if (i == 1) {
        markers <- new_markers
      } else {
        if (pairwise_consensus) {
          markers <- intersect(markers, new_markers)
        } else {
          markers <- unique(c(markers, new_markers))
        }
      }
    }
  } else {
    message("Conducting pairwise DE between:", population1, "and", pop2_intersect)

    object_markers <- FindMarkers(object,
                                  ident.1 = pop1_intersect,
                                  ident.2 = pop2_intersect)

    object_markers <- object_markers %>%
      filter(p_val_adj < pval_cutoff &
               pct.1 > pct_cutoff &
               avg_log2FC > log2FC_cutoff) %>%
      mutate(feature = rownames(.))

    markers <- object_markers$feature

    if (length(markers) < 1) {
      warning("No markers found between", population1, "and", pop2_intersect)
    }
  }

  return(markers)
}



consensusMarkers <- function(objects,
                             method = c("all",
                                        "wilcox_raw",
                                        "wilcox_norm",
                                        "xgboost",
                                        "LDA"),
                             control_object = "None",
                             features = 'highest_variant',
                             metadata,
                             population1,
                             population2 = 'all',
                             pairwise = T,
                             pairwise_consensus = T,
                             object_consensus = T,
                             pval_cutoff = 0.05,
                             pct_cutoff = 0.5,
                             log2FC_cutoff = 1) {

  #'@title pairwiseMarkers
  #'
  #'@description This function calculates pairwise DE across all population2 classes using FindMarkers
  #'
  #'@param object Seurat S4 object
  #'@param method Character, options are "wilcox_raw" wilxocon rank sum test in a
  #'@param metadata Character, metadata to conduct FindMarkers on
  #'@param population1 Character, population 1
  #'@param population2 Character, population 2. IF 'all' then all remaining population will be used. Default is 'all.
  #'@param consensus Logical, whether to output the consensus across all pairwise DE or include all. Default is T
  #'
  #'@import dplyr
  #'
  #'@importFrom Seurat FindMarkers
  #'@importFrom magrittr %>%
  #
  #'@return is a vector of markers
  #'
  #'@export
  #'
  #'@examples
  #' # To get all consensus pairwise markers between Endothelial and remaining populations
  #' markers <- pairwiseMarkers(object = pancreas, metadata = 'seurat_clusters', population1 = 'Endothelial', population2 = 'all', consensus = T)


  # Helper functions:

  pairwiseMarkers <- function(object,
                              metadata,
                              population1,
                              population2,
                              pairwise,
                              pairwise_consensus,
                              pval_cutoff,
                              pct_cutoff,
                              log2FC_cutoff) {

    # Get all classes in metadata
    all_classes <- levels(object@meta.data[[metadata]])

    # Set active ident to metadata first
    object <- SetIdent(object, value = metadata)

    # If population2 is "all" get remaining populations
    if(length(population2) < 2 & any(population2 == 'all')) {
      population2 <- setdiff(all_classes,
                             population1)
    }

    # Check if classes exist
    pop1_intersect <- intersect(population1, all_classes)
    pop1_setdiff <- setdiff(population1, all_classes)

    if(length(pop1_setdiff) > 0) {
      warning(paste("Following population 1 not detected and will not be used in DE:"), pop1_setdiff)
    }

    if(length(pop1_intersect) < 1) {
      stop(paste("None of the population 1 classes exist in metadata:", metadata))

    }

    pop2_intersect <- intersect(population2, all_classes)
    pop2_setdiff <- setdiff(population2, all_classes)

    if(length(pop2_setdiff) > 0) {
      warning(paste("Following population 2 not detected and will not be used in DE:"), pop2_setdiff)
    }

    if(length(pop2_intersect) < 1) {
      stop(paste("None of the population 2 classes exist in metadata:", metadata))

    }

    if(length(pop2_intersect) < 2 & pairwise) {
      warning("Less than two population2 classes detected, simply conducting FindMarkers between", " ", population1, " and", pop2_intersect)
      pairwise <- F
    }

    markers <- c()

    if(pairwise) {

      for(i in 1:length(pop2_intersect)) {

        message("Conducting pairwise DE between:", " ", pop1_intersect, " and ", pop2_intersect[i])

        object_markers <- FindMarkers(object,
                                      ident.1 = pop1_intersect,
                                      ident.2 = pop2_intersect[i])

        object_markers <- object_markers %>%
          filter(p_val_adj < pval_cutoff &
                   pct.1 > pct_cutoff &
                   avg_log2FC > log2FC_cutoff) %>%
          mutate(feature = rownames(.))

        new_markers <- object_markers$feature

        if(length(new_markers) < 1) {
          warning("No markers found between ", pop1_intersect, " and ", pop2_intersect[i])
        }

        if(i == 1) {
          markers <- new_markers
        } else {
          if(pairwise_consensus) {
            markers <- intersect(markers, new_markers)
          } else {
            markers <- unique(c(markers, new_markers))
          }
        }
      }
    } else {
      message("Conducting pairwise DE between:", " ", pop1_intersect, " and ", pop2_intersect)

      object_markers <- FindMarkers(object,
                                    ident.1 = pop1_intersect,
                                    ident.2 = pop2_intersect)

      object_markers <- object_markers %>%
        filter(p_val_adj < 0.05 &
                 pct.1 > 0.1 &
                 avg_log2FC > 0.5) %>%
        mutate(feature = rownames(.))

      markers <- object_markers$feature

      # Warn if no markers
      if(length(markers) < 1) {
        warning("No markers found between ", pop1_intersect, " and ", pop2_intersect)
      }
    }

    return(markers)
  }

  for(j in 1:length(objects)) {

    object <- objects[[j]]

    object_markers <- pairwiseMarkers(object,
                                      metadata,
                                      population1,
                                      population2,
                                      pairwise,
                                      pairwise_consensus)

    if(j == 1) {
      markers <- object_markers
    } else {
      if(object_consensus) {
        markers <- intersect(markers, object_markers)
      } else {
        markers <- unique(c(markers, object_markers))
      }
    }

  }

  if(length(markers) < 1) {
    stop("No consensus markers found across objects.")
  }

  return(sort(markers))

}



#venn <- ggvenn(
  #venn_list,
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #stroke_size = 0.5,
  #set_name_size = 1,
 # show_percentage = F, text_size = 4,
#)

#ggsave(venn, filename = '../figures/venn_all_endothelial_genes.svg', width = 3, height = 3)

#endothelial_genes <- unique(as.vector(unlist(venn_list)))

#saveRDS(endothelial_genes, file = 'tests/testthat/objects/endothelial_genes.rds')
