library(Seurat)
library(slingshot)
library(SingleCellExperiment)


object <- readRDS("output/subset_analyzed/seuratobj_force_magic.rds")
reduction <- 'force'
metadata <- 'leiden_res1'

dimred <- object@reductions[[reduction]]@cell.embeddings
clustering <- as.vector(object@meta.data[,metadata])

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

set.seed(42)
lineages <- getLineages(dimred, clustering)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)

plot(p, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black')

object <- slingshot(object, clusterLabels = 'cluster', reducedDim = 'force')


library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')


lin1 <- getLineages(dimred, clustering, start.clus = 'Capillary 2')



data <- int_atlas

seu <- int_atlas

sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$cluster,
                 start.clus = 'Venule', stretch = 0)



cell_colors <- cell_pal(seu$cluster, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$cluster, hue_pal())



sce <- as.SingleCellExperiment(int_atlas)


sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = "PCA",
                 allow.breaks = FALSE)
# get the lineages:
lnes <- getLineages(reducedDim(sce,"PCA"), sce$cluster)
lnes@lineages


# get the lineages:
lnes <- getLineages(reducedDim(sce,"PCA"), sce$ident)
lnes@lineages


# get the lineages:
lnes <- getLineages(reducedDim(sce,"PCA"), sce$ident)
lnes@lineages

# Save the objects as separate matrices for input in slingshot
dimred <- data@reductions$umap@cell.embeddings
clustering <- data$nn_pca_res.0.2
counts <- data@assays$LOG@counts[data@assays$LOG@var.features,]

# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))


suppressPackageStartupMessages({
  library(slingshot)
})

# Run default Slingshot lineage identification
set.seed(1)

lineages <- getLineages(data = as.data.frame(t(counts)), clusterLabels = clustering)

lineages


#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black')



#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black')
