# Libraries
library(homologene)
library(stringdist)
library(stringr)
library(dplyr)
library(tidyr)

hgnc <- readRDS("rutils/HGNC_Apr2023.rds")


getCorrectInput <- function(input) {

  #' @title getCorrectInput
  #'
  #' @description Checks whether an input is a feature vector, convert to feature vector
  #'
  #' @param input input data
  #'
  #' @return is a character vector
  #'
  #' @export

  # Check type of input
  if (class(input) == 'Seurat') {
    message("Class input 'Seurat' detected")
    features <- rownames(input)
  } else if (class(input) == 'dgCMatrix') {
    message("Class input 'dgCMatrix' detected")
    features <- rownames(input)
  } else if (class(input) == 'character') {
    message("Class input 'character' detected")
    features <- input
  }

  # Check validity of features
  if (!any(is.character(features))) {
    stop("Non-character elements detected in features")
  }

  return(features)

}

# Similarity function based on Levenshtein distance
similarScore <- function(x,y) {

  #' @title similarScore
  #'
  #' @description Scores similarity of two strings
  #'
  #' @param x,y Character, specifying strings to compare
  #'
  #' @importFrom stringdist stringdist
  #'
  #' @export

  score <- (1 - stringdist(x, y, method = "lv")/ max(nchar(x), nchar(y)))

  return(score)

}


convertGenes <- function(input,
                         convert_to = "human",
                         return_full = F) {

  # --- Documentation ------------------------------------------------------- ##

  #' @title convertGenes
  #'
  #' @description This is a function to convert mouse to human genes or vice verse using the homologene package
  #'
  #' @param input Seurat object, dgCMatrix, genes
  #' @param conver_to Character, specifying which genome to convert to. Either of "human/hsapiens", i.e convert mouse to human genes, or "mouse/mmusculus", i.e. convert human to mouse genes.
  #' @param return_full Logical, if TRUE then full conersion table outputed
  #'
  #' @return A vector or seurat object containing converted genes
  #'
  #' @export
  #'
  #' @importFrom homologene mouse2human human2mouse
  #' @importFrom stringr str_to_title

  # --- End of Documentation ------------------------------------------------ ##

  features <- getCorrectInput(input)

  if (convert_to == "human" | convert_to == "hsapiens") {

    genes_table <- mouse2human(features) # Convert using homologene package
    genes_table <- genes_table[,1:2] # Get only first two relevant columns
    colnames(genes_table) = c("original", "updated")

    # Remove ghost genes
    genes_table <- genes_table %>%
      filter(original %in% features)

    # Get absent features
    features_not_present <- setdiff(features, genes_table$original)

    # Make new dataframe with absent features
    unchaged_genes_df <- data.frame("original" = features_not_present,
                                    "updated" = toupper(features_not_present))

    # Bind to original genes_table
    genes_table <- rbind(genes_table, unchaged_genes_df)

    genes_table <- genes_table %>%
      mutate(cased = toupper(original), # Make a case changed column for comparisions
             alternate = cased, # Make an alternate column that will be switched based on iteration
             final = updated, # Make final column that will be assessed for duplicates
             duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F), # Get the first set of final duplicates
             duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F), # Get the first set of original duplicates
             similarity = similarScore(final, cased)) # Get the similarity score between final and "original" or cased

  } else if (convert_to == "mouse" | convert_to == "mmusculus") { # Repeat for mouse

    genes_table <- human2mouse(features)
    genes_table <- genes_table[,1:2]
    colnames(genes_table) = c("original", "updated")

    genes_table <- genes_table %>%
      filter(original %in% features)

    features_not_present <- setdiff(features, genes_table$original)

    unchaged_genes_df <- data.frame("original" = features_not_present,
                                    "updated" = str_to_title(features_not_present)) # str_to_title instead of toupper

    genes_table <- rbind(genes_table, unchaged_genes_df)

    genes_table <- genes_table %>%
      mutate(cased = str_to_title(original), # str_to_title instead of toupper
             alternate = cased,
             final = updated,
             duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F),
             duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F),
             similarity = similarScore(final, cased))
  }


  genes_table_orig_dup <- genes_table %>%
    filter(duplicated_original == F) # Get rid of original duplicates, which are handled later on

  finaldups <- any(genes_table_orig_dup$duplicated_final == T) # Check for final duplicates that are not original duplicates

  counter <- 1 # Counter for iteration

  while(finaldups && counter < 15) {

    # Get all duplicated final genes
    final_duplicates <- genes_table %>%
      filter(duplicated_final == T) %>%
      pull(final)

    # Split into two datasets
    # No duplicates in neither updated or cased (original)
    genes_table_no_dup <- genes_table %>%
      filter(!(updated %in% final_duplicates) & !(cased %in% final_duplicates))

    # With duplicates in either updated or cased (original)
    genes_table_dup <- genes_table %>%
      filter((updated %in% final_duplicates) | (cased %in% final_duplicates)) %>%
      mutate(similarity = similarScore(final, cased)) %>% # recalculate similarity score with original
      group_by(final) %>% # group_by final
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # rank based on similarity score, window function works on groups
      ungroup() %>% # ungroup
      mutate(final = ifelse(rank == 1, final, alternate), # Keep final if similarity rank == 1, change to alternate if not
             alternate = ifelse(final == cased, updated, cased)) # Change alternate based on what final is

    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F)) # Get new duplicated finals


    # Check for duplicates in final again
    genes_table_orig_dup <- genes_table %>%
      filter(duplicated_original == F)

    finaldups <- any(genes_table_orig_dup$duplicated_final == T)

    counter <- counter + 1

  }

  if (any(genes_table$duplicated_original == T)) {
    # Now handle original duplicates
    # Get all original duplicate genes
    original_duplicates <- genes_table %>%
      filter(duplicated_original == T) %>%
      pull(original)

    # Split into two datasets based on duplicates
    genes_table_no_dup <- genes_table %>%
      filter(!(original %in% original_duplicates))

    genes_table_dup <- genes_table %>%
      filter(original %in% original_duplicates) %>%
      mutate(similarity = similarScore(final, cased)) %>% # Get similarity score again
      group_by(cased) %>% # Group by ORIGINAL/Cased
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # Rank similiarty
      slice_max(order_by = desc(rank)) %>% # Remove all other duplicates except first rank
      ungroup() # ungroup

    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F))
  }

  if (any(genes_table$duplicated_final == T)) {

    message("scVascularAnalysis: Warning duplicates still exist")

    inds <- which(genes_table$duplicated_final == T)

    dup_still <- genes_table$original[inds]

    if (length(inds) > 50) {

      stop("scVascularAnalysis: More than 50 genes have duplicates")

    } else {

      message("scVascularAnalysis: Following genes have duplicates still\n",
              paste(dup_still, collapse = ","))

      message("scVascularAnalysis: Converted back to originals")

      genes_table$final[inds] <- genes_table$original[inds]

    }

  }

  genes_table <- genes_table %>%
    arrange(match(original, features)) # Match the original order

  # Return
  if (return_full == T) {

    # Return entire genes_table
    return(genes_table)

  } else if (return_full == F & class(input) == 'Seurat') {
    #If S4 input

    # Get genes
    newgenes <- genes_table$final

    # Get data
    data <- GetAssayData(input, slot = 'counts')

    # Get metadata
    meta <- input@meta.data

    # Get gene names
    rownames(data) <- newgenes

    # Create new object
    object <- CreateSeuratObject(data,
                                 meta.data = meta)
    return(object)

  } else if (return_full == F & class(input) == 'dgCMatrix') {
    #If S4 input

    # Get gene names
    rownames(input) <- genes_table$final

    return(input)

  } else if (return_full == F & class(input) == "character") {
    # If input was character

    return(genes_table$final)

  }

}

updateGenes<- function(input,
                       db = hgnc,
                       return_full = F) {

  # --- Documentation ------------------------------------------------------- ##

  #' updateGenes
  #'
  #'This is a function to update all genes using gprofilers gconvert function
  #'
  #'
  #' @param input Either of character vector or seurat object containint Ggnes to convert
  #' @param db  HGNC symbol dataframe.
  #' @param return_full Logical, if TRUE then full conersion table outputed
  #'
  #'
  #' @importFrom stringr str_split_fixed
  #' @importFrom stringdist stringdist
  #' @importFrom Seurat GetAssayData
  #'
  #' @export

  # --- End of Documentation ------------------------------------------------ ##

  features <- getCorrectInput(input)

  # Make geneDB original variable
  geneDB_orig <- db

  # Get previous symbols
  previous_symbols <- str_split_fixed(geneDB_orig$Previous.symbols, pattern = ", ", n = Inf)

  # Get alias symbols
  alias_symbols <- str_split_fixed(geneDB_orig$Alias.symbols, pattern = ", ", n = Inf)

  # Merge symbols
  all_symbols <- cbind(previous_symbols, alias_symbols)

  # Change symbol names
  colnames(all_symbols) <- mapply(paste0, rep("original_symbol_",
                                              ncol(all_symbols)), c(1:ncol(all_symbols)))

  # Pivot longer
  geneDB <- geneDB_orig %>%
    select(Approved.symbol) %>% # Select only the approved column
    rename(updated = Approved.symbol) %>% # rename to "approved"
    cbind(all_symbols) %>% # cbind all the rest of the symbols
    mutate_all(~ ifelse(. == "", NA, .)) %>% # Make all "" into NAs
    pivot_longer(cols = starts_with("original"),
                 values_to = "original",
                 values_drop_na = T) %>% #pivot longer
    select(updated, original) # Select only approved and alternate columns

  # Get approved features
  features_approved <- unique(intersect(features, geneDB$updated))

  # Get all genes that are not approved
  features_not_approved <- unique(setdiff(features, features_approved))

  # Get all genes that are not present
  features_not_present <- unique(setdiff(features_not_approved, geneDB$original))

  # Get all features that can't be changed
  features_unchaged <- unique(c(features_approved, features_not_present))

  # Get all genes that can be updated
  features_need_changing <- unique(setdiff(features_not_approved, features_not_present))

  # Remove ghost genes
  genes_table <- geneDB %>%
    filter(original %in% features_need_changing)

  # Make new dataframe with unchangeable features
  unchaged_genes_df <- data.frame("original" = features_unchaged,
                                  "updated" = features_unchaged)


  # Bind to original genes_table
  genes_table <- rbind(genes_table, unchaged_genes_df)

  genes_table <- genes_table %>%
    mutate(alternate = original, # Make an alternate column that will be switched based on iteration
           final = updated, # Make final column that will be assessed for duplicates
           duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F), # Get the first set of final duplicates
           duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F), # Get the first set of original duplicates
           similarity = similarScore(final, original)) # Get the similarity score between final and "original" or cased


  genes_table_orig_dup <- genes_table %>%
    filter(duplicated_original == F) # Get rid of original duplicates, which are handled later on

  finaldups <- any(genes_table_orig_dup$duplicated_final == T) # Check for final duplicates that are not original duplicates

  counter <- 1 # Counter for iteration

  while(finaldups && counter < 15) {

    # Get all duplicated final genes
    final_duplicates <- genes_table %>%
      filter(duplicated_final == T) %>%
      pull(final)

    # Split into two datasets
    # No final duplicates in neither updated or  original
    genes_table_no_dup <- genes_table %>%
      filter(!(updated %in% final_duplicates) & !(original %in% final_duplicates))

    # With final duplicates in either updated or original
    genes_table_dup <- genes_table %>%
      filter((updated %in% final_duplicates) | (original %in% final_duplicates)) %>%
      mutate(similarity = similarScore(final, original)) %>% # recalculate similarity score with original
      group_by(final) %>% # group_by final
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # rank based on similarity score, window function works on groups
      ungroup() %>% # ungroup
      mutate(final = ifelse(rank == 1, final, alternate), # Keep final if similarity rank == 1, change to alternate if not
             alternate = ifelse(final == original, updated, original)) # Change alternate based on what final is

    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F)) # Get new duplicated finals


    # Check for duplicates in final again
    genes_table_orig_dup <- genes_table %>%
      filter(duplicated_original == F)

    finaldups <- any(genes_table_orig_dup$duplicated_final == T)

    counter <- counter + 1

  }

  if (any(genes_table$duplicated_original == T)) {

    # Now handle original duplicates
    # Get all original duplicate genes
    original_duplicates <- genes_table %>%
      filter(duplicated_original == T) %>%
      pull(original)

    # Split into two datasets based on duplicates
    genes_table_no_dup <- genes_table %>%
      filter(!(original %in% original_duplicates))

    genes_table_dup <- genes_table %>%
      filter(original %in% original_duplicates) %>%
      mutate(similarity = similarScore(final, original)) %>% # Get similarity score again
      group_by(original) %>% # Group by ORIGINAL/Cased
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # Rank similiarty
      slice_max(order_by = desc(rank)) %>% # Remove all other duplicates except first rank
      ungroup() # ungroup

    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F))

  }

  if (any(genes_table$duplicated_final == T)) {

    stop("Error duplicates still exist")

  }

  genes_table <- genes_table %>%
    arrange(match(original, features)) # Match the original order

  # Return
  if (return_full == T) {

    # Return entire genes_table
    return(genes_table)

  } else if (return_full == F & class(input) == 'Seurat') {
    #If S4 input

    # Get genes
    newgenes <- genes_table$final

    # Get data
    data <- GetAssayData(input, slot = 'counts')

    # Get metadata
    meta <- input@meta.data

    # Get gene names
    rownames(data) <- newgenes

    object <- CreateSeuratObject(data,
                                 meta.data = meta)
    return(object)

  } else if (return_full == F & class(input) == 'dgCMatrix') {
    #If S4 input

    # Get gene names
    rownames(input) <- genes_table$final

    return(input)

  } else if (return_full == F & class(input) == 'character') {
    # If input was character

    return(genes_table$final)

  }


}
