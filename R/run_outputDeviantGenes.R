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
source('R/outputDeviantGenes.R')

# Load data
path <- args[1]

data <- t(readMM(paste0(path,'log1p_matrix.mtx')))
genes <- read.csv(paste0(path,'features.csv'))[,1]

rownames(data) <- genes

deviant_genes <- outpuDeviantGenes(data)

write.csv(deviant_genes, file = paste0(path,'deviant_genes.csv'))
