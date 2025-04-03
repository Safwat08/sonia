# Run Output Deviant Genes.R

# Take command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Input path not valid")
}

# run_outputDeviantGenes
source('R/ec_subtypeScore.R')

# Load data
path <- args[1]

path_obj <- paste0(path, "/seurat_obj.rds")

object <- readRDS(path_obj)

core_name <- sub("^output/([^_]+(?:_[^_]+)?).*", "\\1", path)

object <- SetIdent(object, value = 'leiden_res0_25')

DefaultAssay(object) <- 'MAGIC_imputed_data'

# Load dataframe with signatures
ec_df <- read.csv('rutils/ec_vascular_subtype_signatures_dataframe.csv')[,-1]

ec_subtypeScore(object, 'leiden_res0_25', ec_df, core_name)
