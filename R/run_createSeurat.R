# Run Output Deviant Genes.R

# Take command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Input path not valid")
}

# Load library
library(Matrix)

# run_outputDeviantGenes.R
source('R/createSeurat.R')

# Load data
path <- args[1]

createSeurat(path)
