# Library
library(dplyr)
library(tidyr)
library(ggrepel)
source('R/convertGenes_updateGenes.R')
source('R/addScore.R')
source('R/RunVarimaxPCA.R')
DefaultAssay(object) <- "MAGIC"

# Load object
object <- readRDS('output/subset_analyzed/seuratobj_palantir.rds')

# Load signatures
signature <- readRDS("robjects/.rds")

# Control
signature_control <- sample(object@assays$RNA@counts@Dimnames[[1]], length(signature))
signature_mouse <- convertGenes(signature, convert_to = 'mouse')

# 
object <- addScore(object, score_features = signature_mouse, score_name = 'signature', score_type = c('all'), force_rankings = F)
object <- addScore(object, score_features = signature_control, score_name = 'control', score_type = c('all'), force_rankings = F)

pca_object <- ScaleData(object, features = all_genes)
pca_object <- RunPCA(pca_object, features = all_genes)
pca_object <- RunVarimaxPCA(pca_object, features = all_genes)

pca_df <- pca_object@reductions[['vpca']]@feature.loadings %>%
  as.data.frame() %>%
  mutate(
    genes = rownames(.),
    label = ifelse(
      (rank(VPC_2) <= 10 | rank(-VPC_4) <= 10), 
      genes, 
      NA)
  )

p <- ggplot(pca_df, aes(x=VPC_2, y=VPC_4, label = label)) +
  geom_point(size = 1, color = 'darkred') + 
  geom_label_repel(fill="white", 
                   color = 'black', 
                   max.overlaps = Inf, 
                   show.legend = F,
                   segment.linetype = 3,
                   segment.size = 1,
                   fontface = "bold",
                   size =2.7) +
  theme_classic() +  
  labs( x = "VPC 1",
        y = "VPC 2") + 
  theme(axis.text=element_text(size=8, face = 'bold'),
        axis.title.x = element_text(size=8, face="bold"),
        axis.title.y = element_text(size=8, face="bold"),
        legend.title = element_text(size=8, face = "bold"),
        legend.text = element_text(size=8),
        legend.position = 'none')

# Get feature loadings 
pca_df <- pca_object@reductions[['vpca']]@feature.loadings %>%
  as.data.frame %>%
  mutate(feature = rownames(.))

# DE 
DefaultAssay(object) <- "RNA"

object <- SetIdent(object, value = 'leiden_res0_25')

Markers6 <- FindMarkers(object, ident.1 = '5') %>%
  mutate(feature = rownames(.),
         pct.diff = pct.1 - pct.2,
         norm.pct.diff = pct.1/pct.1 - pct.2/pct.1,
         neglogpval = -log(p_val_adj),
         neglogpval = ifelse(is.infinite(neglogpval), 
                             max(-log(.$p_val_adj)[!is.infinite(-log(.$p_val_adj))]), 
                             neglogpval),
         log2FC.neglogpval = avg_log2FC*neglogpval) %>%
  mutate(color1 = case_when(
    avg_log2FC > 1 & p_val_adj < 0.05 ~ "A",   # Positive differential expression
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "B",        # Negative differential expression
    (avg_log2FC >= -1 & avg_log2FC <= 1) | p_val_adj > 0.05 ~ "C"  # Neutral or non-significant
  ),
  color2 = ifelse(feature %in% signature_mouse, "A", "C")) %>% 
  arrange(factor(color2, levels = c("C", "B", "A")))
         
# Filter the top 10 points where color1 is "steelblue"
top_labels <- Markers6 %>%
  filter(color1 == "A") %>%
  slice_max(log2FC.neglogpval, n = 10)
         
# DE Volcano
vPlot <- ggplot(data = Markers6,  
                aes(x = avg_log2FC, 
                    y = neglogpval,
                    color = color1,
                    label = feature),
                show.legend = F) +
  geom_point(size = 1*1.5,
             color = 'black') +
  geom_point(size = 1) + 
 scale_color_manual(values = c("B" = "#B44656", "A" = "steelblue", "C"= 'grey')) +
  geom_hline(yintercept=c(-log(0.05)), 
             col="black", 
             linetype="dotted", 
             size = 1) +
  geom_vline(xintercept=c(1,-1), 
             col="black", 
             linetype="dotted", 
             size = 1) + 
  geom_label_repel(data = top_labels,  # Use filtered data for labeling
                   aes(label = feature),
                   fill="white", 
                   color = 'black', 
                   max.overlaps = Inf, 
                   show.legend = F,
                   segment.linetype = 3,
                   segment.size = 1,
                   fontface = "bold",
                   size =2.5) +
  theme_classic() +  
  labs( x = "Average log2FC",
        y = "-Log(-adj. p-value)") + 
  theme(axis.text=element_text(size=8, face = 'bold'),
        axis.title.x = element_text(size=8, face="bold"),
        axis.title.y = element_text(size=8, face="bold"),
        legend.title = element_text(size=8, face = "bold"),
        legend.text = element_text(size=8),
        legend.position = 'none')
