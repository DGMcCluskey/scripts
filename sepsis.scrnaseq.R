#Analysis of resolving vs non-resolving sepsis in human PBMCs
#data from https://pubmed.ncbi.nlm.nih.gov/34558746/

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(clustree)
library(SingleR)
library(tidyr)
library(viridis)
library(scCustomize)
library(dplyr)

setwd("C:/Users/dan94/OneDrive - University College London/.UCL Senior Research Fellow/Public datasets/GSE167363_Sepsis_PBMC/Cell_Ranger")

files <- list.files(path = getwd())

files <- files[1:12]

data <- list()


for (i in files) {
  seurat_data.i <- Read10X(data.dir = paste0(i))
  data[[i]] <- CreateSeuratObject(counts = seurat_data.i, min.features = 100, project = i)
}


list2env(data, .GlobalEnv)

rm(data)

sepsis <- merge(S1_T0, y = c(S1_T6, S2_T0, S2_T6, S3_T0, S3_T6, NS1_T0, NS1_T6, NS2_T0, NS2_T6,
                              HC1, HC2), 
                 add.cell.ids = c("S1_T0", "S1_T6", "S2_T0", "S2_T6", "S3_T0", "S3_T6", "NS1_T0", "NS1_T6", 
                                  "NS2_T0", "NS2_T6",
                                  "HC1", "HC2"), project = "Sepsis")

table(sepsis$orig.ident)
sepsis$sample <- sepsis$orig.ident
table(sepsis$sample)

#merging layers of object (Seurat v5 leaves each sample in their own layer)
sepsis <- JoinLayers(sepsis)
sepsis

rm(list = c("S1_T0", "S1_T6", "S2_T0", "S2_T6", "S3_T0", "S3_T6", "NS1_T0", "NS1_T6", 
            "NS2_T0", "NS2_T6",
            "HC1", "HC2"))

#Quality control and processing----------------------------------------------------------------------------------------------------

sepsis[["percent.mt"]] <- PercentageFeatureSet(sepsis, pattern = "^MT-")
VlnPlot(sepsis, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, group.by = "sample",
        ncol = 1)
sepsis <- subset(sepsis, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

sepsis <- NormalizeData(sepsis, normalization.method = "LogNormalize", scale.factor = 10000)

sepsis <- FindVariableFeatures(sepsis, selection.method = "vst", nfeatures = 2000)

patient <- c("S1_T0" = "S1", "S1_T6" = "S1", "S2_T0" = "S2", "S2_T6" = "S2", 
             "S3_T0" = "S3", "S3_T6" = "S3", "NS1_T0" = "NS1", "NS1_T6" = "NS1", 
             "NS2_T0" = "NS2", "NS2_T6" = "NS2",
             "HC1" = "HC1", "HC2" = "HC2")

Idents(sepsis) <- "sample"

sepsis <- RenameIdents(sepsis, patient)

sepsis@active.ident

sepsis <- StashIdent(sepsis, save.name = "patient")

#all.genes <- rownames(sepsis)
sepsis <- ScaleData(sepsis, features = VariableFeatures(sepsis), 
                    vars.to.regress = c("nCount_RNA", "percent.mt", "patient"))

saveRDS(sepsis, file = "sepsis.data.processed.rds")

#PCA

sepsis <- readRDS("sepsis.data.processed.rds")
sepsis <- RunPCA(sepsis, features = VariableFeatures(object = sepsis))

palette <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
before.correction <- DimPlot(sepsis, reduction = "pca", split.by = "sample", ncol = 4)+
  scale_colour_manual(values = palette)
before.correction

after.correction <- DimPlot(sepsis, reduction = "harmony", group.by = "patient", 
                            split.by = "sample", ncol = 4)+
  scale_colour_manual(values = palette)
after.correction

before.correction+after.correction

#harmony correction
sepsis <- RunHarmony(sepsis, c("patient"))


ElbowPlot(sepsis, ndims = 50, reduction = "harmony")

#Clustering

sepsis <- FindNeighbors(sepsis, dims = 1:15, reduction = "harmony")
sepsis <- FindClusters(sepsis, resolution = 1)
sepsis <- FindClusters(sepsis, resolution = 0.9)
sepsis <- FindClusters(sepsis, resolution = 0.8)
sepsis <- FindClusters(sepsis, resolution = 0.7)
sepsis <- FindClusters(sepsis, resolution = 0.6)
sepsis <- FindClusters(sepsis, resolution = 0.5)
sepsis <- FindClusters(sepsis, resolution = 0.4)
sepsis <- FindClusters(sepsis, resolution = 0.3)
sepsis <- FindClusters(sepsis, resolution = 0.2)

sepsis <- RunUMAP(sepsis, dims = 1:15, reduction = "harmony")

clustree(sepsis)

Idents(sepsis) <- "RNA_snn_res.0.5"

p <- DimPlot(sepsis, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, label.box = T)
p

ggsave(plot = p, filename = "UMAP.res.0.5.png", dpi = 500,
       height = 4, width = 6)

#marker inspection
DimPlot(sepsis, reduction = "umap", split.by = "sample", ncol = 4)
DotPlot(sepsis, features = c("CD3D", "CD4", "LTB", "SELL", "CCR7", "IL7R", "CD8A", "NKG7", "GNLY", 
                                 "MS4A1", "CD79A", "IGHD", "CD27", "CD38", "IGHG1", 
                             "LYZ", "CD68", "CLEC9A", "CLEC10A", "IRF7", "CD14", "FCGR3A", "CSF3R", "JAML", "PPBP", "HBB", "CD34", "MKI67"), 
        scale.by = "size", scale = 10)+
  scale_colour_viridis(option = "F", direction = -1)+
  #geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme(legend.text = element_text(size = 20), axis.text = element_text(size =24, colour = "black"),
        legend.title = element_text(size = 24), 
        axis.text.x = element_text(angle = 0))+ 
  theme(axis.text = element_text(size = 14, colour = "black"))+ylab(NULL)+xlab(NULL)+
  coord_flip()

FeaturePlot(sepsis, features = "PPBP")

#annotation
annotation <- c("0" = "Naive T", "1" = "B", "2" = "CD14+mono", "3" = "NK", "4" = "Platelets",
                "5" = "CD16+mono", "6" = "Memory T", "7" = "CD14+mono", "8" = "Neutrophil",
                "9" = "Plasmablast", "10" = "RBC", "11" = "HSC", "12" = "DC", "13" = "Debris",
                "14" = "CD8+T")

Idents(sepsis) <- "RNA_snn_res.0.5"

sepsis <- RenameIdents(sepsis, annotation)

sepsis@active.ident

sepsis <- StashIdent(sepsis, save.name = "labels")

sepsis$labels <- factor(sepsis$labels, levels = c("Naive T", "Memory T", "CD8+T", "NK",
                                                  "B", "Plasmablast", "CD14+mono", "CD16+mono", "DC", 
                                                  "Neutrophil", "Platelets",
                                                  "HSC", "RBC", "Debris"))

palette2 <- c(palette, "grey50", "black", "bisque3")

p <- DimPlot(sepsis, reduction = "umap", group.by = "labels")+
  scale_colour_manual(values = palette2)

p

ggsave(plot = p, filename = "sepsis.UMAP.all.clusters.png", dpi = 500,
       height = 4, width = 6)

table(sepsis$labels, sepsis$sample)

DotPlot(sepsis, features = c("CD3D", "CD4", "LTB", "SELL", "CCR7", "IL7R", "CD8A", "NKG7", "GNLY", 
                             "MS4A1", "CD79A", "IGHD", "CD27", "CD38", "IGHG1", 
                             "LYZ", "CD68", "CLEC9A", "CLEC10A", "IRF7", "CD14", "FCGR3A", "CSF3R", "JAML", "PPBP", "HBB", "CD34", "MKI67"), 
        scale.by = "size", scale = 10)+
  scale_colour_viridis(option = "F", direction = -1)+
  #geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme(legend.text = element_text(size = 20), axis.text = element_text(size =24, colour = "black"),
        legend.title = element_text(size = 24), 
        axis.text.x = element_text(angle = 0))+ 
  theme(axis.text = element_text(size = 14, colour = "black"))+ylab(NULL)+xlab(NULL)+
  coord_flip()

saveRDS(sepsis, file = "sepsis.clustered.rds", compress = F)

#removing platelets, RBC, HSC, neutro and debris and re-running processing
sepsis.clean <- subset(sepsis, idents = c("Naive T", "Memory T", "CD8+T", "NK",
                                           "B", "Plasmablast", "CD14+mono", "CD16+mono", "DC"))


sepsis.clean <- FindVariableFeatures(sepsis.clean, selection.method = "vst", nfeatures = 2000)
sepsis.clean <- ScaleData(sepsis.clean, features = VariableFeatures(sepsis.clean), 
                    vars.to.regress = c("nCount_RNA", "percent.mt", "patient"))
sepsis.clean <- RunPCA(sepsis.clean, features = VariableFeatures(object = sepsis.clean))
sepsis.clean <- RunHarmony(sepsis.clean, c("patient"))
ElbowPlot(sepsis.clean, ndims = 50, reduction = "harmony")

sepsis.clean <- FindNeighbors(sepsis.clean, dims = 1:22, reduction = "harmony")
sepsis.clean <- FindClusters(sepsis.clean, resolution = 1)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.9)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.8)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.7)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.6)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.5)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.4)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.3)
sepsis.clean <- FindClusters(sepsis.clean, resolution = 0.2)

sepsis.clean <- RunUMAP(sepsis.clean, dims = 1:22, reduction = "harmony")

clustree(sepsis.clean)

Idents(sepsis.clean) <- "RNA_snn_res.0.8"

saveRDS(sepsis.clean, file = "sepsis.cleaned.rds", compress = F)

p <- DimPlot(sepsis.clean, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, label.box = T)
p
ggsave(plot = p, filename = "UMAP.cleaned.res.0.8.png", dpi = 500,
       height = 4, width = 6)

p <- DimPlot(sepsis.clean, reduction = "umap", group.by = "patient", split.by = "sample", ncol = 4)+
  scale_colour_manual(values = palette)
p

FeaturePlot(sepsis.clean, features = c("CD3D", "CD8A", "CCR7", "SELL", "IL7R", "LTB"), ncol = 3)&
  scale_colour_viridis(option = "F", direction = -1)

table(sepsis.clean$RNA_snn_res.0.6)

DotPlot(sepsis.clean, features = c("CD3D", "CD4", "CD40LG", "LTB", "SELL", "CCR7", "IL7R", "CCL5",
                                   "CD8A", "NKG7", "GNLY", 
                             "MS4A1", "CD79A", "IGHD", "CD27", "CD38", "CD24", "CR2", "IGHG1", "IGHM", 
                             "LYZ", "CD68", "CD1C", "CLEC10A", "IRF7", "CD14", "FCGR3A", "S100A12", 
                             "CSF3R", "PPBP", "HBB"), 
        scale.by = "size", scale = 10)+
  scale_colour_viridis(option = "F", direction = -1)+
  #geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme(legend.text = element_text(size = 20), axis.text = element_text(size =24, colour = "black"),
        legend.title = element_text(size = 24), 
        axis.text.x = element_text(angle = 0))+ 
  theme(axis.text = element_text(size = 14, colour = "black"))+ylab(NULL)+xlab(NULL)+
  coord_flip()

annotation <- c("0" = "Naive T", "1" = "Naive B", "2" = "CD14+ mono", "3" = "NK", "4" = "CD14+CD16+ mono", 
                "5" = "Cytotoxic T",
                "6" = "Memory B", "7" = "Memory T", "8" = "Memory T", "9" = "CD14+ mono", 
                "10" = "Transitional B", "11" = "Doublets",
                "12" = "DC", "13" = "Cytotoxic T", "14" = "NK", "15" = "pDC", "16" = "Plasmablasts")

Idents(sepsis.clean) <- "RNA_snn_res.0.6"

sepsis.clean <- RenameIdents(sepsis.clean, annotation)

sepsis.clean@active.ident

sepsis.clean <- StashIdent(sepsis.clean, save.name = "labels")

sepsis.clean$labels <- factor(sepsis.clean$labels, levels = c("Naive T", "Memory T", "Cytotoxic T",
                                                              "NK", "Transitional B", "Naive B",
                                                              "Memory B", "Plasmablasts",
                                                              "CD14+ mono", "CD14+CD16+ mono",
                                                              "DC", "pDC", "Doublets"))

#removing debris cluster from dataset, but not re-clustering again
Idents(sepsis.clean) <- "labels"
sepsis.clean <- subset(sepsis.clean, idents = c("Naive T", "Memory T", "Cytotoxic T",
                                                "NK", "Transitional B", "Naive B",
                                                "Memory B", "Plasmablasts",
                                                "CD14+ mono", "CD14+CD16+ mono",
                                                "DC", "pDC"))


palette <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
palette2 <- c(palette, "grey50", "black", "bisque3")
p <- DimPlot(sepsis.clean, reduction = "umap", group.by = "labels")+
  scale_colour_manual(values = palette2)

p
ggsave(plot = p, filename = "sepsis.clean.annotated.png", dpi = 500,
       height = 4, width = 6)

saveRDS(sepsis.clean, file = "sepsis.cleaned.annotated.rds", compress = F)
#using singler to support annotation

ref.data <- MonacoImmuneData()

ref.data$label.fine

input <- as.SingleCellExperiment(sepsis.clean)

predictions <- SingleR(test = input, ref = ref.data, labels = ref.data$label.fine)

table(predictions$labels)

sepsis.clean[["SingleR.labels"]] <- predictions$labels

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

p <- DimPlot(sepsis.clean, reduction = "umap", group.by = "SingleR.labels")+
  scale_colour_manual(values = col_vector)
p

ggsave(plot = p, filename = "sepsis.clean.singler.annotation.png", dpi = 500,
       height = 4, width = 8)


#assessing confidence of the singler annotation
predictions$scores
plotScoreHeatmap(predictions)

predictions$scores

plotDeltaDistribution(predictions, show = "delta.next")

sepsis.clean[["Singler.confidence"]] <- predictions$delta.next


to.remove <- is.na(predictions$pruned.labels)
table(label = predictions$labels, removed = to.remove)


#inspection of ETS2 (and other RIP TFs)

p <- FeaturePlot(sepsis.clean, features = c("ETS2"), split.by = "sample",
                     keep.scale = "feature", combine = F)

cowplot::plot_grid(plotlist = p, ncol = 4)

myeloid <- subset(sepsis.clean, idents = c("DC", "CD14+CD16+ mono", "CD14+ mono"))
data <- DotPlot(myeloid, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEC",
                                      "JDP2", "MYC"), group.by = "labels", split.by = "sample",
        cols = col_vector)

data <- data$data

data <- separate(data, col = id, into = c("celltype", "sample"), sep = "_", extra = "merge")

ggplot(data, aes(x = features.plot, y = sample, fill = avg.exp.scaled,
                 size = pct.exp))+geom_point(shape = 21, colour = "black")+
  facet_wrap(~celltype)+
  scale_fill_viridis(option = "F", direction = -1)+
  theme(legend.text = element_text(size = 20), axis.text = element_text(size =24, colour = "black"),
        legend.title = element_text(size = 24), 
        axis.text.x = element_text(angle = 0))+ylab("Sample")+xlab(NULL)+theme_bw()

data2 <- separate(data, col = sample, into = c("patient", "time"), sep = "_")
data2[is.na(data2)] <- 0

data2$time <- factor(data2$time, levels = c("0", "T0", "T6"),
                     labels = c("0", "0", "6"))

data2$group <- factor(data2$patient, levels = c("NS1", "NS2", "S1", "S2", "S3", "HC1", "HC2"),
                      labels = c("NS", "NS", "S", "S", "S", "HC", "HC"))

#data2$time <- as.character(data2$time) 
#data2$time <- as.numeric(data2$time)

data2$gene.cell <- paste0(data2$celltype, "-", data2$features.plot)

ggplot(data2, aes(x = time, y = avg.exp, colour = group, group = patient))+geom_point(size = 3)+
  facet_wrap(~gene.cell, scales = "free_y", ncol = 7)+geom_line()+theme_bw()+
  scale_colour_manual(values = c("royalblue3", "brown2", "grey30"))

#dotplot version
data2$cell.time <- paste0(data2$celltype, "-", data2$time)
ggplot(data2, aes(x = cell.time, y = patient, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(shape = 21, colour = "black")+
  facet_wrap(~features.plot, ncol = 7)+
  scale_fill_viridis(option = "F", direction = -1)+ylab("Sample")+xlab(NULL)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, size = 10, colour = "black", hjust = 1))

#getting change in expression between time 0 and 6

data3 <- filter(data2, patient %in% c("NS1", "NS2", "S1", "S2", "S3"))

data3 <- data3[, -c(8,10,11)]

data3 <- pivot_wider(data3, names_from = time, values_from = c(avg.exp, pct.exp, avg.exp.scaled))

data3$exp.change <- data3$avg.exp_6-data3$avg.exp_0
data3$pct.change <- data3$pct.exp_6-data3$pct.exp_0
data3$scaled.change <- data3$avg.exp.scaled_6-data3$avg.exp.scaled_0

limit <- max(abs(data3$scaled.change)) * c(-1, 1)

limit2 <- max(abs(data3$pct.change)) * c(-1, 1)

ggplot(data3, aes(x = celltype, y = patient, fill = scaled.change, size = pct.change))+
  geom_point(shape = 21, colour = "black")+
  facet_wrap(~features.plot, ncol = 7)+
  scale_fill_distiller(palette = "RdBu", limit = limit)+ylab("Sample")+xlab(NULL)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, size = 10, colour = "black", hjust = 1))+
  scale_size_manual(limit = limit2)

#adding ETS2+ cells as metadata

sepsis.clean[["ETS2+"]] <- 
