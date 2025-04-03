## --- changeClusterName ---------------------------------------------------- ##

changeClusterName <- function(object,
                              old_name,
                              new_name,
                              replace = FALSE,
                              orig_col,
                              new_col = "cluster") {

  # --- Documentation ----------------------------------------------------- ##

  #' changeClusterName
  #'
  #'  This function changes the name of a cluster in a Seurat object.
  #'
  #' @param object Seurat S4 object, a Seurat object containing metadata
  #' @param old_name Character, specifying the old name of the cluster
  #' @param new_name Character, specifying the new name of the cluster
  #' @param replace Logical, whether to replace the original metadata column or create a new one
  #' @param orig_col Character, the name of the original metadata column to modify
  #' @param new_col Character, the name of the new metadata column to create if replace is FALSE
  #' @return Seurat object with updated cluster names in metadata
  #' @export

  # --- End of Documentation ---------------------------------------------- ##

  meta <- object@meta.data

  # Ensure the original column exists in metadata
  if (!orig_col %in% colnames(meta)) {
    stop(paste("Error: The column", orig_col, "does not exist in the Seurat object metadata."))
  }

  if (replace) {
    # Replace directly in the original column
    meta[[orig_col]] <- ifelse(meta[[orig_col]] == old_name, new_name, meta[[orig_col]])
  } else {
    # Create a new column with updated cluster names
    meta[[new_col]] <- ifelse(meta[[orig_col]] == old_name, new_name, meta[[orig_col]])
  }

  # Update the Seurat object with modified metadata
  object@meta.data <- meta

  # Return the modified object
  return(object)
}
