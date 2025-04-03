# Library
library(dplyr)
library(tidyr)
library(ggrepel)
source('R/convertGenes_updateGenes.R')
source('R/addScore.R')
source('R/markerDotPlot.R')

# Load object
object <- readRDS('output/subset_analyzed/seuratobj_palantir.rds')


organs <- c("brain", "testis", "liver", "spleen", "lung", "kidney", "colon", "small_intestine", "heart", "EDL", "soleus", "pancreas", "organ")

is_any_capitalized <- function(vec) {
  any(grepl("^[A-Z]+$", vec))
}

for (i in 1:11) {
  
  markers <- readxl::read_excel('rutils/kalucka_murine_ecatlas_allvascularmarkers.xlsx', sheet = i)
  
  for (j in 1:ncol(markers)) {
    genes <- as.vector(na.omit(markers[,j][[1]]))
    
    if(is_any_capitalized(genes)) {
      
      genes <- convertGenes(genes, convert_to = 'mouse')
    }
    
    name <- gsub(" ", "_", names(markers[,j]))
    
    name <- paste0(organs[i], "_", name)
    
    object <- addScore(object, genes, score_name = name, score_type = 'seurat')
  }
  
}

meta <- colnames(object@meta.data)[18:95]

p <- betterDotPlot(object, features = meta, metadata = 'leiden_res1')

# Integrated Atlas DE Analysis
int_atlas_de <- readRDS('robjects/integrated_atlas_mast_de.rds') %>%
  mutate(cluster = case_when(
    cluster == "Arteriole" ~ "Arteriole-Precapillary",
    cluster == "Artery 1" ~ "Artery-Large",
    cluster == "Artery 2" ~ "Arteriole-Feeding",
    cluster == "Capillary 1" ~ "Capillary-Islet",
    cluster == "Capillary 2" ~ "Capillary-Exocrine",
    cluster == "Venule" ~ "Venule-Postcapillary"))

# Expanded signature genes
sig_genes_exp <- readRDS('robjects/signature_genes_expanded.rds')

# Signature genes for each vascular subpopulation
int_atlas_de_filt <- int_atlas_de %>%
  mutate(feature = gene,
         pct.diff = pct.1 - pct.2,
         norm.pct.diff = pct.1/pct.1 - pct.2/pct.1,
         neglogpval = -log(p_val_adj),
         neglogpval = ifelse(is.infinite(neglogpval), 
                             max(-log(.$p_val_adj)[!is.infinite(-log(.$p_val_adj))]), 
                             neglogpval),
         log2FC.neglogpval = avg_log2FC*neglogpval) %>%
  filter(p_val_adj < 0.05 &
           avg_log2FC > 0.5 & 
            feature %in% sig_genes_exp) %>%
  group_by(feature) %>% 
  slice_max(order_by = avg_log2FC, n = 1) %>% # remove duplicates by keeping the gene for the cluster with the highest log2FC for that gene
  ungroup()

# Get list for each subpopulation 
marker_list <- int_atlas_de_filt %>%
  group_by(cluster) %>%
  dplyr::summarise(Markers = list(feature)) %>%
  pull(Markers)

#marker_list <- int_atlas_de_filt %>%
# group_by(cluster) %>%
#slice_max(order_by = avg_log2FC, n = 20)%>%
#dplyr::summarise(Markers = list(feature)) %>%
#pull(Markers)

marker_list <- marker_list[-3]

marker_list <- marker_list[c(1,2,4,3,5)]

names(marker_list) <- c("ArterioleFeeding", "ArteriolePrecapillary","CapillaryIslet", "CapillaryExocrine", "VenulePostcapillary")

# Signature genes for each vascular subpopulation
cap_markers <- readRDS('robjects/integrated_atlas_capillary_subset_analysis_mast_de.rds') %>%
  mutate(feature = gene,
         pct.diff = pct.1 - pct.2,
         norm.pct.diff = pct.1/pct.1 - pct.2/pct.1,
         neglogpval = -log(p_val_adj),
         neglogpval = ifelse(is.infinite(neglogpval), 
                             max(-log(.$p_val_adj)[!is.infinite(-log(.$p_val_adj))]), 
                             neglogpval),
         log2FC.neglogpval = avg_log2FC*neglogpval) %>%
  filter(p_val_adj < 0.05 &
           avg_log2FC > 0.25 & 
           feature %in% sig_genes_exp) %>% 
  group_by(feature) %>% 
  slice_max(order_by = avg_log2FC, n = 1) %>% # remove duplicates by keeping the gene for the cluster with the highest log2FC for that gene
  ungroup()

marker_list_cap <- cap_markers %>%
  group_by(cluster) %>%
  dplyr::summarise(Markers = list(gene)) %>%
  pull(Markers)

names(marker_list_cap) <- c("CapillaryIslet", "CapillaryExocrine")

marker_list$CapillaryExocrine <- intersect(marker_list$CapillaryExocrine, marker_list_cap$CapillaryExocrine)

marker_list$CapillaryIslet <- intersect(marker_list$CapillaryIslet, marker_list_cap$CapillaryIslet)

# Score names 
score_names <- paste0(names(marker_list))

# Remove COL4A1 and MCAM
marker_list$ArteriolePrecapillary <- setdiff(marker_list$ArteriolePrecapillary, c("MCAM", "COL4A1"))

# Convert to mouse

marker_list <- lapply(marker_list, convertGenes, convert_to = 'mouse')

for (i in 1:length(marker_list)) {
  
  object <- addScore(object, score_features = marker_list[[i]], score_name = score_names[i], score_type = c("seurat", "ucell"))
  
}


meta <- colnames(object@meta.data)[94:ncol(object@meta.data)]

plotlist <- list

for (i in 1:11) {
  
  markers <- readxl::read_excel('rutils/kalucka_murine_ecatlas_allvascularmarkers.xlsx', sheet = i)
  
  for (j in 1:ncol(markers)) {
    name <- gsub(" ", "_", names(markers[,j]))
    
    name <- paste0(organs[i], "_", name)
    
    score_name[j] <- paste0(name, "_score_seurat")
    
    
    plotlist[[i]] <- FeaturePlot(object, features = )
    
  }
  
}


meta <- object_sub@meta.data

features <- colnames(meta[,22:ncol(meta)])

DotPlot(object_sub, features = features)