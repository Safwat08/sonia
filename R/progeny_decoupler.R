library(decoupleR)
library(progeny)
library(tidyverse)
library(dplyr)
library(tidyr)
library(Seurat)
source('R/convertGenes_updateGenes.R')
source('R/addScore.R')

organs <- c("brain", "testis", "liver", "spleen", "lung", "kidney", "colon", "small_intestine", "heart", "EDL", "soleus", "pancreas", "organ")

is_any_capitalized <- function(vec) {
  any(grepl("^[A-Z]+$", vec))
}

all_genes <- c()

ec_df <- tibble(EC = character(),
                Genes = character())

for (i in 1:13) {

  markers <- readxl::read_excel('rutils/kalucka_murine_ecatlas_allvascularmarkers.xlsx', sheet = i)

  for (j in 1:ncol(markers)) {
    genes <- as.vector(na.omit(markers[,j][[1]]))

    if(i == 12 | (i == 13 & j == 11)) {

      genes <- convertGenes(genes, convert_to = 'mouse')
    }

    all_genes <- unique(c(all_genes,genes))

    name <- gsub(" ", "_", names(markers[,j]))

    name <- paste0(organs[i], "_", name)

    df = tibble(rep(name, length(genes)), genes)

    ec_df = rbind(ec_df, df)

    object <- addScore(object, genes, score_name = name, score_type = 'seurat')
  }

}

write.csv(ec_df, file = 'rutils/all_ec_genes_longer.csv')

#
gene_trends = palantir.presults.compute_gene_trends(
  ad,
  expression_key="MAGIC_imputed_data",
)


object <- readRDS('output/subset_analyzed/seuratobj_force_magic.rds')

write.csv(model, file = 'rutils/progeny_model_mouse_500.csv')

mat <- GetAssayData(object, slot = 'data', assay = 'RNA')

mat_sub <- mat[as.vector(na.omit(match(all_genes, rownames(mat)))),]

acts <- run_mlm(mat=mat, net=model, .source='pathway', .target='gene',
                .mor='weight', minsize = 5)

# Extract mlm and store it in pathwaysmlm in data
pathobj <- object
pathobj[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = pathobj) <- "pathwaysmlm"

# Scale the data
pathobj <- ScaleData(pathobj)
pathobj@assays$pathwaysmlm@data <- pathobj@assays$pathwaysmlm@scale.data


# TF

net <- read.csv('rutils/collectri_mouse.csv')

acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)


# Extract mlm and store it in pathwaysmlm in data
tfobj <- object
tfobj[['tfulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = tfobj) <- "tfulm"

# Scale the data
tfobj <- ScaleData(tfobj)
tfobj@assays$tfulm@data <- tfobj@assays$tfulm@scale.data


# Genes of interest
FeaturePlot(object, reduction = 'force', features = c('Igfbp7', 'Igfbp3', 'Insr', "Plvap", "Ramp3", "Plcb1", "Hspg2", "Col4a1", "Col4a2", "Card10", "Robo4", "Meox1","Pcdh17", "Col15a1", "Mcam", "Nedd9", "Plxnd1", "Cd34", "Mcf2l", "Tbc1d1", "Cdc42ep1", "Vwa1", "Rapgef4"),
        order = T, cols = c('lightgrey','red'))

+
  theme_classic() +
  theme( plot.title = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(size =8,face = "bold", colour = "black", vjust = 1, hjust= 1, angle = 90),
         axis.text.y = element_text(size =8,face = "bold", colour = "black"),
         legend.text = element_text(size =8, angle = 0, face = "bold", colour = "black", vjust = 0.5),
         legend.title = element_text(size =8, angle = 0, face = "bold", colour = "black", vjust = 0.5),
         legend.key.width = unit(0.5, 'cm'),
         legend.key.height = unit(0.3, 'cm'),
         legend.position = 'right') +
  guides(color = guide_colorbar(title = "Avg. Exp",
                                order = 1),
         size = guide_legend(title = "% Exp",
                             order = 2))
