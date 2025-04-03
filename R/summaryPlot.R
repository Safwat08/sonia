markers <- tabulamuris_ecatlas_pancreas_markers
meta <- "cluster"
cluster <- "Islet-like Capillary"
pval_cutoff <- 0.05
logfc_cutoff <- 0.25
label_n <- 5

# Preprocess data
# Get metadata from object and convert to data frame and order by frequency
data <- object@meta.data %>%
  rename("clusters" = !!sym(meta)) %>% # !!sym to rlang
  group_by(clusters) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(percent = 100*n/sum(n),
         x = "x") %>%
  arrange(desc(percent)) %>%
  arrange(clusters != cluster) %>% #Clever hack to reorder
  mutate(cols = do_ColorPalette("steelblue", nrow(.)))


# Change factor levels of active ident
my_levels <- data$clusters
object <- SetIdent(object, value = meta)
object@active.ident<- factor(x = object@active.ident, levels = my_levels)

# Make relevant 
markers <- markers %>%
  mutate(feature = rownames(.),
         neglogpval = -log(p_val_adj),
         color = ifelse(avg_log2FC > logfc_cutoff & p_val_adj < pval_cutoff, "yes", "no"),
         label = ifelse(rank(-avg_log2FC) <= label_n & p_val_adj < pval_cutoff, feature, NA))


max_notinf <- max(markers$neglogpval[!is.infinite(markers$neglogpval)])

markers <- markers %>%
  mutate(neglogpval = ifelse(is.infinite(neglogpval), max_notinf, neglogpval))

# DimPlot of object
colors <- data$cols
names(colors) <- data$clusters
cluster_n <- length(unique(data$clusters))
leg_nrow <- ceiling(cluster_n/2)

dPlot <- do_DimPlot(object,
                    pt.size = 0.05,
                    legend.icon.size = 2.5,
                    font.size = 10, 
                    colors.use = colors, 
                    legend.ncol = 2, 
                    legend.nrow = leg_nrow)

# Percentage stacked bar chart showing population percentages

# Get data from object and convert to dataframe
data <- as.data.frame(table(object@meta.data$cluster, object@meta.data$dataset))

# Change colnames
colnames(data) <- c("Cluster", "Dataset", "Freq")

# Change order of levels 
data$Cluster <- factor(data$Cluster, levels= my_levels)

# Plot DimPlot
dPlot <- do_DimPlot(int_atlas, colors.use = color.scheme, legend.icon.size = 3, font.size = 8)

# Plot percentage stacked bar chart
bPlot <- ggplot(data[], 
                aes(fill=Cluster, 
                    y=Freq, 
                    x=Dataset)) + 
  ylab("Cell Ratios") +
  geom_bar(position="fill", 
           stat="identity", 
           color = "black", 
           show.legend = F) + 
  theme_classic() + 
  scale_fill_discrete(breaks= my_levels) + 
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = c("Descartes-Fetal", 
                              "Tosti-Adult", 
                              "Tosti-Neonatal")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   size=8, 
                                   face="bold", color = "black"),
        axis.text.y = element_text(size=8, 
                                   face="bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8, 
                                    face="bold", color = "black"))

bPlot <- ggplot(data, 
                aes(fill=factor(clusters, levels = clusters), # rearrange
                    y=percent,
                    x=x)) + 
  ylab("Celltype Ratios") +
  geom_bar(position="fill", stat="identity", color = 'black') + 
  scale_fill_manual(values = data$cols) + 
  labs(fill = "Cluster") +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=8, face="bold"),
        axis.text.y = element_text(size = 8, face = "bold", colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = 'none')

# Volcano plot of DE genes

vPlot <- ggplot(data = markers,  
                aes(x = avg_log2FC, 
                    y = neglogpval,
                    color = color,
                    label = label), 
                show.legend = F) +
  geom_point(size = 1*1.5,
             color = 'black') +
  geom_point(size = 1) + 
  scale_color_manual(values = c("yes" = "#B44656", "no" = "steelblue")) +
  geom_hline(yintercept=c(-log(0.05)), 
             col="black", 
             linetype="dotted", 
             size = 1) +
  geom_vline(xintercept=c(0.25), 
             col="black", 
             linetype="dotted", 
             size = 1) + 
  geom_label_repel(fill="white", 
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


# Use plot_grid to combine all plots
p <- plot_grid(dPlot, bPlot, vPlot, nrow=1, ncol=3, rel_widths = c(3,1,3))

ggsave(p, 
       filename = "../figures/tabulamuris_ecatlas_umap_bar_volcano.png", 
       width = 8, 
       height = 4)
```