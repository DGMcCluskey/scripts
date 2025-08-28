library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(clustree)
library(SingleR)
library(tidyr)
library(viridis)
library(dplyr)
library(data.table)
library(scCustomize)
library(broom)
library(purrr)
library(reticulate)
library(SeuratWrappers)
library(gridExtra)
library(cowplot)
library(patchwork)
library(readr)
library(GeneOverlap)
library(VennDiagram)
library(ggVennDiagram)
library(ggupset)
library(scDblFinder)
library(celda)#decontx package apparently
library(ggbump)
library(slingshot)

#slingshot pseudotime on ali flumap data from Claire Smith
setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Claire_Smith_data")

ali <- readRDS("ALI_Flu_250626_seurat.rds")

table(ali$Ziegler_majority_voting)
table(ali$manual_annotation)


cluster.palette <- c("indianred2", "brown3", "hotpink", "coral", "darkred", "red", "violet", "mediumvioletred", "darkorchid", "darkorchid4", "orange", "goldenrod2", "grey50", "black",
                     "mediumseagreen", "forestgreen", "darkgreen", "green", "steelblue", "blue", "darkblue", "royalblue", "cyan", "cyan4")





p <- DimPlot(ali, group.by = "manual_annotation")+scale_colour_manual(values = cluster.palette)
p

ggsave(plot = p, filename = "ali_umap_by_experiment.png", dpi = 500,
       height = 6, width = 14)
DimPlot(ali, group.by = "Experiment")

p <- FeaturePlot(ali, features = c("H5N1_Turkey", "H5N1_Texas_Cattle", "H1N1_England"), pt.size = 2, ncol = 3)&
  scale_colour_viridis(option = "F", direction = -1)
p

ggsave(plot = p, filename = "viral_gene_expression.png", dpi = 500,
       height = 4, width = 12)


#extract elements for slingshot input
umap <- Embeddings(ali, reduction = "X_umap")
clusters <- ali$manual_annotation  # or seu$your_cluster_column

sce <- SingleCellExperiment(
  assays = list(logcounts = as.matrix(ali@assays$originalexp@data)),
  reducedDims = SimpleList(UMAP = umap)
)

# Add cluster labels as a colData column
colData(sce)$cluster <- clusters

rm(ali)
#run slingshot
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

SlingshotDataSet(sce)
slingLineages(sce)

# Extract UMAP coordinates
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- sce$cluster
umap_df$pseudotime <- slingPseudotime(sce)[, 1]  # First lineage
umap_df$cell <- rownames(umap_df)

# Extract lineage curves
curves <- slingCurves(sce)

# For plotting curves, convert curve lines to data.frame
curve_list <- lapply(curves, function(curve) {
  data.frame(curve$s[, 1:2])  # assuming 2D embedding
})
names(curve_list) <- paste0("Lineage", seq_along(curves))

# Plot with ggplot
#with all curves
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", na.value = "grey") +
  theme_bw() +
  lapply(seq_along(curve_list), function(i) {
    geom_path(data = curve_list[[i]], aes(x = Xumap_1, y = Xumap_2),
              color = "black", size = 1, inherit.aes = FALSE)
  }) +
  ggtitle("Slingshot Trajectories with Pseudotime")

#with one curve
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", na.value = "grey") +
  theme_bw() +
    geom_path(data = curve_list[[1]], aes(x = Xumap_1, y = Xumap_2),
              color = "black", size = 1, inherit.aes = FALSE) +
  ggtitle("Slingshot Trajectories with Pseudotime")


cluster.palette <- c("indianred2", "brown3", "hotpink", "coral", "darkred", "red", "violet", "mediumvioletred", "darkorchid", "darkorchid4", "orange", "goldenrod2", "grey50", "black",
                    "mediumseagreen", "forestgreen", "darkgreen", "green", "steelblue", "blue", "darkblue", "royalblue", "cyan", "cyan4")

cluster.palette2 <- c("mediumseagreen", "indianred2", "darkorchid3", "darkred", "cyan4", "goldenrod2",
                      "blue", "orange", "green", "steelblue2", "hotpink")

#colour by cluster
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.8)+
  theme_bw() + scale_colour_manual(values = cluster.palette)+
  lapply(seq_along(curve_list), function(i) {
    geom_path(data = curve_list[[i]], aes(x = Xumap_1, y = Xumap_2),
              color = "black", size = 1, inherit.aes = FALSE)
  }) +
  ggtitle("Slingshot Trajectories with Pseudotime")

p

ggsave(plot = p, filename = "pseudotime_all.png", dpi = 500,
       height = 5, width = 10)

saveRDS(sce, file = "ALI_sce_slingshot.rds", compress = F)
sce <- readRDS(file = "ALI_sce_slingshot.rds")

#here i am getting cell ids and making the umap and curve dataframes again. this allows me to plot the curves but only colour the clusters of relevance in each trajectory
library(SingleCellExperiment)
library(slingshot)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Extract pseudotime and cell_id
pseudotime_df <- as.data.frame(slingPseudotime(sce))
pseudotime_df$cell_id <- rownames(pseudotime_df)

# 2. Reshape pseudotime to long format (one row per cell-curve combo)
pseudotime_long <- pseudotime_df %>%
  pivot_longer(cols = starts_with("Lineage"), names_to = "curve", values_to = "pseudotime")

# 3. Extract UMAP and metadata
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- rownames(umap_df)
umap_df$cluster <- sce$cluster  # add cluster directly

# 4. Join all data into one plotting dataframe
plot_df <- pseudotime_long %>%
  left_join(umap_df, by = "cell_id") %>%
  mutate(highlight = ifelse(is.na(pseudotime), "Other", as.character(cluster)))

# 5. Get curve paths from slingshot object
curve_list <- slingCurves(sce)
curve_df <- bind_rows(
  lapply(seq_along(curve_list), function(i) {
    crv <- curve_list[[i]]
    data.frame(
      Xumap_1 = crv$s[, 1],
      Xumap_2 = crv$s[, 2],
      curve = paste0("Lineage", i)
    )
  })
)

curve_df$curve <- factor(curve_df$curve, levels = c("Lineage1", "Lineage2", "Lineage3", "Lineage4", "Lineage5", "Lineage6",
                                                    "Lineage7", "Lineage8", "Lineage9", "Lineage10", "Lineage11", "Lineage12"))


plot_df <- plot_df %>%
  mutate(highlight = ifelse(is.na(pseudotime), "Other", as.character(cluster)))

cluster_order <- levels(plot_df$cluster)
cluster_order <- c(cluster_order, "Other")

plot_df$highlight <- factor(plot_df$highlight, levels = cluster_order)

cluster.palette2 <- c(cluster.palette, "lightgrey")

# 6. Plot: facet by curve and color only involved clusters
p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = highlight), size = 0.6, alpha = 0.8) +
  geom_path(data = curve_df, aes(x = Xumap_1, y = Xumap_2, group = curve),
            inherit.aes = FALSE, color = "black", size = 1) +
  scale_color_manual(values = cluster.palette2) +
  facet_wrap(~curve) +
  theme_bw() +theme(
    legend.text = element_text(size = 14),     # Increase legend text size
    legend.title = element_text(size = 16),    # Increase legend title size
    legend.key.size = unit(1.5, "lines")       # Increase size of legend keys (dots/squares)
  )+
  ggtitle("Slingshot Lineages Faceted by Curve (Cluster-Colored)")
p

ggsave(plot=p, filename = "pseudotime_umap_by_lineage.png", dpi = 500,
       height = 10, width = 16)


#gene expression over pseudotime=========================================
# Extract pseudotime
pseudotime_df <- as.data.frame(slingPseudotime(sce))
pseudotime_df$cell_id <- rownames(pseudotime_df)

# Choose genes to plot
genes_of_interest <- c("TNF", "IFNG")  # replace with your genes

# Extract gene expression
expr_df <- as.data.frame(t(logcounts(sce)[genes_of_interest, ]))
expr_df$cell_id <- rownames(expr_df)

# Combine data
plot_df <- pseudotime_df %>%
  pivot_longer(cols = starts_with("Lineage"), names_to = "curve", values_to = "pseudotime") %>%
  filter(!is.na(pseudotime)) %>%  # only keep cells in the trajectory
  left_join(expr_df, by = "cell_id")

ggplot(filter(plot_df, curve == "Lineage11"), aes(x = pseudotime, y = TNF)) +
  #geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "TNF expression along Lineage11 pseudotime",
       x = "Pseudotime", y = "Log-normalized expression") +
  theme_bw()
