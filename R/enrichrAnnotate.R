library(enrichR)
library(stringr)

enrichAnnotateR <- function(input,
                            meta_name = 'default',
                            cluster_name = 'default',
                            compare_to = 'all',
                            database = 'annot',
                            custom_dbs = NULL,
                            pct_cutoff = 0.1,
                            log2fc_cutoff = 1,
                            pval_cutoff = 0.05,
                            ngenes = 50) {

  #' @title enrichAnnotateR
  #'
  #' @description This function annotates a clusters based on gene set enrichment analysis
  #'
  #' @param input Seurat object or vector of genes
  #' @param database Character, specifies whether to compare to Annotation dbs: 'annot', function dbs: 'function', 'all' dbs on enrichR., 'custom' for custom dbs, enter custom_dbs parameter
  #' @param meta_name Character, specifies the metadata of cluster to annotate if input is seurat object
  #' @param cluster_name Character, specifies the cluster to annotate if input is seurat object
  #' @param pct_cutoff,log2fc_cutoff,pval_cutoff,ngenes Numerics, cutoffs for percentage expression in cluster, average log2FC, adjusted p value and top 'n' genes to use for enrichr.
  #'
  #' @import enrichR
  #' @importFrom enrichR enrichr
  #' @importFrom Seurat FindMarkers
  #' @importFrom magrittr %>%
  #' @importFrom dplyr mutate filter arrange
  #' @importFrom stringr str_split_fixed
  #'
  #' @details enrichAnnotateR identifies the top markers using Seurat's FindMarkers
  #' then conducts GSEA using enrichR with the following databases: "CellMarker_2024",
  #' "Azimuth_2023", "Tabula_Sapiens", "Descartes_Cell_Types_and_Tissue_2021", "Tabula_Muris".
  #'
  #' @return is a dataframe of the most likely annotations and corresponding enrichr GSEA metrics
  #'
  #' @export

  if (class(input) == 'Seurat') {

    message("Seurat class detected, taking input seurat object")
    object <- input

    # free space
    rm(input)

    # Check if metadata exists in object
    if (!(meta_name %in% colnames(object@meta.data))) {
      stop("Metadata name not found in object")
    }

    # Get the list of available clusters from the object metadata
    available_clusters <- unique(object@meta.data[[meta_name]])

    # Identify indices of missing cluster_name entries
    missing_idx1 <- which(!(cluster_name %in% available_clusters))

    # Identify indices of missing compare_to entries
    missing_idx2 <- which(!(compare_to %in% available_clusters))

    # Check if cluster_name clusters exist in the object
    if (length(missing_idx1) > 0) {
      warning("Clusters ", paste(cluster_name[missing_idx1], collapse = ", "), " not present")
    }

    if (length(missing_idx1) == length(cluster_name)) {
      stop("None of the input clusters found in metadata")
    }

    # Check if compare_to clusters exist in the object (if compare_to is not 'all')
    if (!all(compare_to == 'all')) {
      if (length(missing_idx2) > 0) {
        warning("Clusters ", paste(compare_to[missing_idx2], collapse = ", "), " not present")
      }

      if (length(missing_idx2) == length(compare_to)) {
        stop("None of the compare clusters found in metadata")
      }
    }

    # Set ident to metadata
    object <- SetIdent(object, value = meta_name)

    # Find Markers
    if (!all(compare_to == 'all')) {
      # Compare to rest of clusters by default
      Markers <- FindMarkers(object, ident.1 = cluster_name)
    } else {
      # Else compare to custom input
      Markers <- FindMarkers(object, ident.1 = cluster_name, ident.2 = compare_to)
    }

    # Filter Markers
    Markers_filt <- Markers %>%
      filter(pct.1 > pct_cutoff &
               avg_log2FC > log2fc_cutoff &
               p_val_adj < pval_cutoff) %>%
      slice_max(order_by = avg_log2FC, n = ngenes)

    # Need to ensure Marker_filt is not errored

    # Get gene names
    genes <- rownames(Markers_filt)

  } else if (class(input) == 'character') {
    # Take input as genes directly
    message("Character class detected, taking input as genes directly")
    genes <- input
  }

  # Check number of genes
  if (length(genes) > 1) {
    message(paste(length(genes), "genes being used"))
  } else {
    stop("No genes detected")
  }

  if (database == 'annot') {

    # Input data sets
    dbs_pathway <- c("CellMarker_2024",
                     "Azimuth_2023",
                     "Tabula_Sapiens",
                     "Descartes_Cell_Types_and_Tissue_2021",
                     "Tabula_Muris")

  } else if (database == 'function') {

    dbs_pathway <- c("WikiPathway_2023_Human",
                     "Reactome_2022",
                     "MSigDB_Hallmark_2020",
                     "KEGG_2021_Human",
                     "GO_Biological_Process_2023")

  } else if (database == 'all') {

    db <- listEnrichrDbs()

    dbs_pathway <- db$libraryName

  } else if(database == 'custom') {

    db <- listEnrichrDbs()

    dbs <- db$libraryName

    dbs_pathway <- intersect(custom_dbs, dbs)

  }

  # Get enrichr output
  enrichr_output <- enrichr(genes, dbs_pathway)

  # Find sets without any terms
  enrichr_output_nogenes <- as.vector(unlist(lapply(enrichr_output, function(x) nrow(x) < 1)))

  # If there are sets with no terms, remove them
  if (any(enrichr_output_nogenes == T)) {
    enrichr_output <- enrichr_output[-c(which(enrichr_output_nogenes == T))]
  }

  # Write in the database name for the different enrichr_output list elements
  for (i in 1:length(enrichr_output)) {
    enrichr_output[[i]]$Database <- names(enrichr_output[i])
  }

  # Filter important terms
  enrichr_output <- Reduce(rbind, enrichr_output) %>% # bind all database outputs
    mutate(intersection_size = as.numeric(str_split_fixed(Overlap, "/", 2)[,1])) %>% #include an intersection size
    mutate(term_size = as.numeric(str_split_fixed(Overlap, "/", 2)[,2])) %>% #include a term size
    mutate(gene_ratio = intersection_size / term_size) %>% #include a gene ratio
    filter(intersection_size > 1) %>% # filter terms with only 1 gene
    filter(Adjusted.P.value < 0.05) %>% #filter terms wih an adjusted p value > 0.05
    arrange(Adjusted.P.value) # arrange by adjusted pvalue

  return(enrichr_output)

}
