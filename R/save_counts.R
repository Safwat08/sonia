## Save as anndata
# Laod Library
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are passed
if (length(args) == 0) {

  stop("Provide at least 1 argument specifying dir of rds object")

  } else if (length(args) > 2) {

    print("Only two arguments accepted, rds object and save dir, rest of arguments will not be used")
}

# Arg1 is location of rds object to convert
arg1 <- as.character(args[1])

# Arg2 is save directory
if (length(args) > 1) {

  arg2 <- as.character(args[2])

} else {

  arg2 <- ""
}

# Print object
print(paste("Loading object", arg1))

# Read object
object <- readRDS(arg1)

# Get metadata and counts
meta <- object@meta.data
counts <- GetAssayData(object, assay = 'RNA', slot = 'counts')

# Get file name
file <- arg1 %>%
  sub(".*/", "", .) %>%
  gsub("\\.rds$", "", .) %>% # Remove .rds extension
  gsub("\\.", "_", .) %>% # replace dots with _
  tolower(.)


# Save files

# Dir name
slash <- substr(arg2, nchar(arg2), nchar(arg2))

# Add "/" to dir if not present
if (slash != "/") {

  save_dir <- paste0(arg2, "/")

} else {

  save_dir <- arg2
}

# Print save dir
print(paste("---> Saving to dir", save_dir))

# Save files
writeMM(counts, file = paste0(save_dir, file, "_counts.mtx"))

print("---> Counts file written")

write.csv(colnames(counts), file = paste0(save_dir, file, "_barcodes.csv"))

print("---> Barcodes file written")

write.csv(rownames(counts), file = paste0(save_dir, file, "_features.csv"))

print("---> Features file written")

write.csv(meta, file = paste0(save_dir, file, "_meta.csv"))

print("---> Meta file written")
