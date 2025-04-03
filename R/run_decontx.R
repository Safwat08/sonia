# Load library
library(Matrix)
library(decontX)

# Load data
path <- 'output/subset_raw/'

data <- t(readMM(paste0(path,'matrix.mtx')))
genes <- read.csv(paste0(path,'features.csv'))[,1]

rownames(data) <- genes

decont_results <- decontX(data)

decont_counts <- decont_results$decontXcounts

writeMM(t(decont_counts), file = paste0(path,'decont_matrix.mtx'))
