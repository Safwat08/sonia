library(Seurat)
library(AUCell)
library(UCell)

addScore <- function(input,
                     score_features,
                     score_name = 'feature',
                     score_type = 'all',
                     force_rankings = FALSE,
                     store_rankings = TRUE,
                     ...) {

  # --- Documentation ------------------------------------------------------- ##
  #' @title addScore
  #'
  #' @description Adds a collective score for a specified gene set to a Seurat object or dgCMatrix.
  #'
  #' @param input Seurat object or dgCMatrix with features as rows and samples as columns.
  #' @param score_features Character vector specifying genes to score.
  #' @param score_name Character specifying the score name to be stored in metadata.
  #' @param score_type Character specifying which type of scoring to use. Options are "all", "seurat", "aucell", or "ucell".
  #' @param store_rankings Logical indicating whether to store rankings in the Seurat object (default: TRUE).
  #' @param force_rankings Logical indicating whether to force recalculation of rankings (default: FALSE).
  #' @param ... Additional parameters for individual scoring methods.
  #'
  #' @importFrom UCell StoreRankings_UCell ScoreSignatures_UCell
  #' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
  #' @importFrom Seurat CreateSeuratObject AddModuleScore FetchData GetAssayData
  #'
  #' @return A Seurat object with scores stored in metadata, or a list with scores and rankings if input is a dgCMatrix.
  #'
  #' @export
  #'
  #' @example
  # --- End of Documentation ------------------------------------------------ ##

  # Check input type
  if (inherits(input, 'Seurat')) {
    message("Seurat object detected")
    input_type <- "seurat"
    exp_matrix <- GetAssayData(input, slot = "data") # Retrieve the RNA assay data matrix
    cells <- colnames(exp_matrix) # Extract cell names
    # Check for existing UCell and AUCell rankings in the Seurat object
    misc <- names(input@assays[["RNA"]]@misc)
    rank_ucell_present <- "ranks_ucell" %in% misc
    rank_aucell_present <- "ranks_aucell" %in% misc
  } else if (inherits(input, 'dgCMatrix')) {
    message("dgCMatrix detected")
    input_type <- "matrix"
    exp_matrix <- input
    cells <- colnames(exp_matrix)
    rank_ucell_present <- FALSE
    rank_aucell_present <- FALSE
  } else {
    stop("Input must be a Seurat object or dgCMatrix")
  }

  # Determine which scores to calculate
  calc_seurat <- calc_ucell <- calc_aucell <- FALSE
  if ("all" %in% score_type) {
    calc_seurat <- calc_ucell <- calc_aucell <- TRUE
  }
  if ('seurat' %in% score_type) {
    calc_seurat <- TRUE
  }
  if ('ucell' %in% score_type) {
    calc_ucell <- TRUE
  }
  if ('aucell' %in% score_type) {
    calc_aucell <- TRUE
  }

  # Initialize an empty dataframe to store scores
  output_meta <- data.frame(row.names = cells)

  # Seurat Module Score Calculation
  if (calc_seurat) {
    message("Calculating Seurat module scores")
    seurat_object <- CreateSeuratObject(exp_matrix)
    seurat_object <- AddModuleScore(seurat_object, features = list(score_features), name = 'score', slot = 'counts')
    scores <- FetchData(seurat_object, vars = "score1")
    colnames(scores) <- paste0(score_name, "_score")
    output_meta <- cbind(output_meta, scores)
  }

  # UCell Score Calculation
  if (calc_ucell) {
    if (rank_ucell_present && !force_rankings) {
      ranks_ucell <- input@assays[["RNA"]]@misc$ranks_ucell
      message("Using pre-existing UCell rankings")
    } else {
      ranks_ucell <- StoreRankings_UCell(exp_matrix)
      message("Generated UCell rankings")
    }
    scores <- ScoreSignatures_UCell(matrix = exp_matrix, features = list(score_features), precalc.ranks = ranks_ucell)
    colnames(scores) <- paste0(score_name, "_score_ucell")
    output_meta <- cbind(output_meta, scores)
  }

  # AUCell Score Calculation
  if (calc_aucell) {
    if (rank_aucell_present && !force_rankings) {
      ranks_aucell <- input@assays[["RNA"]]@misc$ranks_aucell
      message("Using pre-existing AUCell rankings")
    } else {
      ranks_aucell <- AUCell_buildRankings(exp_matrix, plotStats = FALSE)
      message("Generated AUCell rankings")
    }
    list_score_features <- list("score" = score_features)
    auc_scores <- AUCell_calcAUC(geneSets = list_score_features, rankings = ranks_aucell)
    auc_data <- t(as.data.frame(auc_scores@assays@data$AUC))
    colnames(auc_data) <- paste0(score_name, "_score_aucell")
    output_meta <- cbind(output_meta, auc_data)
  }

  # Store the scores in the appropriate format
  if (input_type == "seurat") {
    message("Storing scores in Seurat metadata")
    input <- AddMetaData(input, output_meta)
    if (store_rankings) {
      if (calc_ucell) {
        message("Storing UCell rankings in Seurat object")
        input@assays[["RNA"]]@misc$ranks_ucell <- ranks_ucell
      }
      if (calc_aucell) {
        message("Storing AUCell rankings in Seurat object")
        input@assays[["RNA"]]@misc$ranks_aucell <- ranks_aucell
      }
    }
    return(input)
  } else if (input_type == "matrix") {
    message("Returning scores and rankings as a list")
    output_list <- list("Scores" = output_meta)
    if (store_rankings) {
      if (calc_ucell) {
        output_list$ranks_ucell <- ranks_ucell
      }
      if (calc_aucell) {
        output_list$ranks_aucell <- ranks_aucell
      }
    }
    return(output_list)
  }
}
