#investigating resolution inhibitory pathways (RIPs) in sepsis looking at multiple publically available sepsis datsets
#in collaboration with Vincent
#script by Dan McCluskey
setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/RIPs_Vincent_project/samples")

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

#LOADING IN DATA===============================================================================================================================================
#first will load in two public datasets at once (Darden et al and Qiu et al) as these have the barcode/feature/matrix file structure
files <- list.files(path = getwd())
files <- files[1:21]
files
data <- list()

for (i in files) {
  # Read 10X data
  seurat_data <- Read10X(data.dir = i)
  # Handle both list (multi-assay) and matrix (single-assay) formats
  if (is.list(seurat_data)) {
    # Prefer "Gene Expression" if present
    if ("Gene Expression" %in% names(seurat_data)) {
      counts <- seurat_data[["Gene Expression"]]
    } else {
      warning(paste("No 'Gene Expression' assay found for", i, "- skipping this sample"))
      next  # skip this sample
    }
  } else {
    # Single assay: assume it's gene expression
    counts <- seurat_data
  }
  
  # Rename cells to include sample ID to avoid barcode collision
  colnames(counts) <- paste(i, colnames(counts), sep = "_")
  # Create Seurat object
  data[[i]] <- CreateSeuratObject(counts = counts, 
                                  min.features = 100, 
                                  project = i)
}

sample_names <- names(data)

#give each sample its own name
for (i in seq_along(data)) {
  data[[i]]$orig.ident <- sample_names[i]
}

#merge all samples into one seurat object called sepsis
sepsis <- merge(
  x = data[[1]],
  y = data[-1],
  add.cell.ids = sample_names,
  project = "Sepsis"
)

rm(data)

table(sepsis$orig.ident)
sepsis$sample <- sepsis$orig.ident
table(sepsis$sample)

sepsis

#merging layers of object (Seurat v5 leaves each sample in their own layer)
sepsis <- JoinLayers(sepsis)
sepsis

saveRDS(sepsis, file = "sepsis.darden.qiu.rds", compress = F)
sepsis <- readRDS("sepsis.darden.qiu.rds")

#read in the Reyes et al dataset (expression data is in a single csv, so this needs to be treated differently)
reyes <- fread("scp_gex_matrix_raw.csv")
reyes_meta <- read.csv("scp_meta_updated.csv")

# Assuming genes are rows, first column is gene names, rest are numeric
exp <- as.data.frame(reyes)
rownames(exp) <- exp[[1]]
exp[[1]] <- NULL

# Convert to matrix then to sparse matrix
exp <- as(as.matrix(exp), "dgCMatrix")
rm(reyes)

rownames(reyes_meta) <- reyes_meta$cell_id
reyes_meta <- reyes_meta[, -c(1)]

#i only want to keep samples that are total pbmc (not DC enriched) - i already subsetted the metadat excel before loading in, so now subsetting the expression matrix
keep <- rownames(reyes_meta)
exp_final <- exp[, colnames(exp) %in% keep]

#make seurat object, check numbers and make a sample column that is the condition + donor
reyes_seurat <- CreateSeuratObject(counts = exp_final, meta.data = reyes_meta)
table(reyes_seurat$condition, reyes_seurat$donor_id)
reyes_seurat$sample <- paste0(reyes_seurat$condition, "_", reyes_seurat$donor_id)
table(reyes_seurat$sample)
reyes_seurat

saveRDS(reyes_seurat, "sepsis.reyes.rds", compress = F)

#MERGE OBJECTS================================================================================================================================================================
sepsis$source <- ifelse(grepl("Qiu", sepsis$sample),
                        "Qiu",
                        "Darden")

table(sepsis$source)

sepsis$condition <- ifelse(grepl("HC", sepsis$sample),
                           "HC", "Sepsis")

table(sepsis$condition)

table(sepsis$sample)

sepsis$donor_id <- sepsis$sample

reyes_seurat$condition[reyes_seurat$condition == "Control"] <- "HC"

table(reyes_seurat$condition)

reyes_seurat$source <- "Reyes"

sepsis_all <- merge(sepsis, reyes_seurat)
sepsis_all

table(sepsis_all$condition, sepsis_all$donor_id)

table(sepsis_all$source)

sepsis_all <- JoinLayers(sepsis_all)
sepsis_all

rm(sepsis)
rm(reyes_seurat)
rm(exp)
rm(exp_final)
rm(counts)
rm(seurat_data)

#INTEGRATION================================================================================================================================================================

#splitting object into layers depending on experimental source
sepsis_all[["RNA"]] <- split(sepsis_all[["RNA"]], f = sepsis_all$source)
sepsis_all

#run standard pipeline to view un-integrated outcome
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sepsis_all[["percent.mt"]] <- PercentageFeatureSet(sepsis_all, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(sepsis_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
sepsis_all <- subset(sepsis_all, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

sepsis_all <- NormalizeData(sepsis_all)
sepsis_all <- FindVariableFeatures(sepsis_all)
sepsis_all <- ScaleData(sepsis_all)
sepsis_all <- RunPCA(sepsis_all)

sepsis_all <- FindNeighbors(sepsis_all, dims = 1:30, reduction = "pca")
sepsis_all <- FindClusters(sepsis_all, resolution = 1, cluster.name = "unintegrated_clusters")

sepsis_all <- RunUMAP(sepsis_all, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
p <- DimPlot(sepsis_all, reduction = "umap.unintegrated", group.by = c("source", "condition"))
p

ggsave(plot = p, filename = "uncorrected.sepsis.umaps.png", dpi = 500,
       height = 4, width = 10, path = "../plots")

sepsis_all

sepsis_all <- IntegrateLayers(sepsis_all, method = HarmonyIntegration,
                              layers = "counts", scale.layer = "scale.data")

sepsis_all

saveRDS(sepsis_all, file = "sepsis_combined.rds", compress = F)

#run clustering on harmony embedding
sepsis_all <- FindNeighbors(sepsis_all, dims = 1:30, reduction = "harmony")
sepsis_all <- FindClusters(sepsis_all, resolution = 1, cluster.name = "clusters")

sepsis_all <- RunUMAP(sepsis_all, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
p <- DimPlot(sepsis_all, reduction = "umap.harmony", group.by = c("source", "condition"))
p

ggsave(plot = p, filename = "harmony.sepsis.umaps.png", dpi = 500,
       height = 4, width = 10, path = "../plots")



#adding a new metadata column with tentative disease groups
#based on reyes manuscript, will consider bac-sep and icu-sep as sepsis and all else as other (and then maybe exclude other)

sepsis_all$condition2 <- ifelse(grepl("Sepsis", sepsis_all$condition),
                                "Sepsis", ifelse(grepl("-SEP", sepsis_all$condition), "Sepsis",
                                                 ifelse(grepl("HC", sepsis_all$condition), "HC",
                                                        "Other")))

table(sepsis_all$condition2, sepsis_all$condition)

#MARKER EXPRESSION===============================================================================================================================================
FeaturePlot(sepsis_all, reduction = "umap.harmony", 
            features = c("CD3D", "CD79A", "LYZ", "NKG7", "IRF7", "ETS2"))&scale_colour_viridis(option = "B")

DimPlot(sepsis_all, reduction = "umap.harmony", group.by = "source", split.by = "condition2")

#run a number of resolutions to compare clusters
sepsis_all <- FindClusters(sepsis_all, resolution = 0.3, cluster.name = "clusters_0.3")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.4, cluster.name = "clusters_0.4")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.5, cluster.name = "clusters_0.5")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.6, cluster.name = "clusters_0.6")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.7, cluster.name = "clusters_0.7")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.8, cluster.name = "clusters_0.8")
sepsis_all <- FindClusters(sepsis_all, resolution = 0.9, cluster.name = "clusters_0.9")
sepsis_all <- FindClusters(sepsis_all, resolution = 1, cluster.name = "clusters_1.0")



DimPlot(sepsis_all, reduction = "umap.harmony", group.by = c("clusters_0.3", "clusters_0.4", "clusters_0.5"), pt.size = 2)

VlnPlot(sepsis_all, features = c("CD3D", "CD8A", "CD4", "FOXP3", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                 "CD1C", "IRF7", "PPBP", "HBB"), pt.size = 0,
        group.by = "clusters_0.3", stack = T)
#violin plot above indicates that clusters 8 and 13 have high platelet marker expression, whilst 16 and 11 are RBCs. will exclude and re-cluster
Idents(sepsis_all) <- "clusters_0.3"
sepsis <- subset(sepsis_all, idents = c("8", "13", "16", "11"), invert = T)

DimPlot(sepsis, reduction = "umap.harmony", group.by = c("clusters_0.3", "clusters_0.4", "clusters_0.5"), pt.size = 2)

#RE-CLUSTERING AND RE-INTEGRATING====================================================================================================================================
sepsis[["RNA"]] <- split(sepsis[["RNA"]], f = sepsis$source)
sepsis

sepsis <- NormalizeData(sepsis)
sepsis <- FindVariableFeatures(sepsis)
sepsis <- ScaleData(sepsis)
sepsis <- RunPCA(sepsis)

sepsis <- FindNeighbors(sepsis, dims = 1:30, reduction = "pca")
sepsis <- FindClusters(sepsis, resolution = 0.5, cluster.name = "unintegrated_clusters")

sepsis <- RunUMAP(sepsis, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

rm(sepsis_all)
sepsis_final <- IntegrateLayers(sepsis, method = HarmonyIntegration,
                              layers = "counts", scale.layer = "scale.data")

saveRDS(sepsis_final, file = "sepsis_final.rds", compress = F)

sepsis_final <- FindNeighbors(sepsis_final, dims = 1:30, reduction = "harmony")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.5, cluster.name = "clusters")

sepsis_final <- RunUMAP(sepsis_final, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
p <- DimPlot(sepsis_final, reduction = "umap.harmony", group.by = c("source", "condition2"))
p

sepsis_final <- JoinLayers(sepsis_final)

#also save as an anndata object to be used in python (for celltypist)
library(zellkonverter)
sce <- as.SingleCellExperiment(myeloid, assay = "RNA")
table(sce$sepsis_state, sce$labels)

assayNames(sce)
assays(sce, withDimnames = FALSE)$counts_layer <- assays(sce)$counts

writeH5AD(sce, "myeloid.h5ad")
rm(sce)

#also need to manually extract the nieghbor graph so that exact clustering can be used with celltypist (and to avoid running pca again which is causing memory trouble in python)
# Extract the shared nearest neighbor graph (usually used for clustering)
graph <- sepsis_final[["RNA_snn"]]

# Convert to a matrix or sparse matrix
graph_matrix <- as(graph, "dgCMatrix")

# Save to file (Matrix Market format is compatible with Python)
Matrix::writeMM(graph_matrix, file = "snn_graph.mtx")

#also save cell names
write.csv(colnames(sepsis_final), file = "cell_names.csv", row.names = FALSE)


#MARKER EXPRESSION AND CLUSTER ANNOTATION============================================================================================================================
FeaturePlot(sepsis_final, reduction = "umap.harmony", 
            features = c("CD3D", "CD79A", "LYZ", "NKG7", "IRF7", "ETS2"))&scale_colour_viridis(option = "B")

#run a number of resolutions to compare clusters
sepsis_final <- FindClusters(sepsis_final, resolution = 0.3, cluster.name = "clusters_0.3")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.4, cluster.name = "clusters_0.4")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.5, cluster.name = "clusters_0.5")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.6, cluster.name = "clusters_0.6")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.7, cluster.name = "clusters_0.7")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.8, cluster.name = "clusters_0.8")
sepsis_final <- FindClusters(sepsis_final, resolution = 0.9, cluster.name = "clusters_0.9")
sepsis_final <- FindClusters(sepsis_final, resolution = 1, cluster.name = "clusters_1.0")

clustree(sepsis_final, prefix = "clusters_")


p1 <- DimPlot(sepsis_final, reduction = "umap.harmony", group.by = c("clusters_0.5"), pt.size = 2, label = T)
p1

ggsave(plot = p1, filename = "umap.clusters.res.0.5.png", dpi = 500,
       height = 5, width = 7, path = "../plots/")



FeaturePlot(sepsis_final, features = c("LYZ", "CD14", "FCGR3A", "TYROBP", "C1QC", "CSF3R", "CD68", "CSF1R", "MKI67", "TOP2A", "CX3CR1",
                                       "VCAN"), reduction = "umap.harmony")&
  scale_colour_viridis(option = "F", direction = -1)
DotPlot(sepsis_final, features = c("LYZ", "CD14", "FCGR3A", "TYROBP", "C1QC", "CSF3R", "CD68", "CSF1R", "MKI67", "TOP2A", "CX3CR1",
                                   "VCAN"), group.by = "clusters_0.5", 
        scale.by = "size",
        dot.min = 0.01)+ scale_colour_viridis(option = "F", direction = -1)



p1 <- DimPlot(sepsis_final, reduction = "umap.harmony", group.by = c("labels", "condition2", "source"), pt.size = 2)&
  scale_colour_brewer(palette = "Set1")
p1

p2 <- DimPlot(sepsis_final, reduction = "umap.harmony", split.by = "source", pt.size = 2)&
  scale_colour_brewer(palette = "Set1")
p2

table(sepsis_final$condition2, sepsis_final$donor_id)


ggsave(plot = p, filename = "sepsis.umaps.png", dpi = 500,
       height = 5, width = 14, path = "../plots/")


VlnPlot(sepsis_final, features = c("CD3D", "CD8A", "CD4", "FOXP3", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                 "CD1C", "IRF7", "ETS2", "CxCL8", "IL6", "IL1B", "IFNG"), pt.size = 0,
        group.by = "clusters_0.3", stack = T)


Proportion_Plot(sepsis_final, group_by_var = "clusters_0.3", split.by = "source")

DotPlot(sepsis_final, features = c("CD3D", "CD8A", "CD4", "FOXP3", "TRDV2", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                   "CD1C", "IRF7", "C1QA", "CLEC10A", "CD68", "CSF3R",
                                   "MKI67", "CxCL8", "IL6", "IL1B", "IFNG"), group.by = "clusters_0.5", scale.by = "size",
        dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)



Idents(sepsis_final) <- "clusters_0.5"
degs_res_0.5 <- FindAllMarkers(sepsis_final, logfc.threshold = log(1.5), min.pct = 0.2, only.pos = T)
degs_res_0.5 <- filter(degs_res_0.5, p_val_adj < 0.05)

annotation <- c("0" = "CD14+mono", "1" = "CD14+mono", "2" = "CD4+T", "3" = "CD8+T", "4" = "CD14+mono", 
                "5" = "Mac/mono", "6" = "NK",
                "7" = "B", "8" = "CD4+T", "9" = "CD14+mono", "10" = "CD8+T", "11" = "Neutrophil", "12" = "CD14+mono", 
                "13" = "DC",
                "14" = "pDC", "15" = "CD8+T", "16" = "CD14+mono", "17" = "B", "18" = "B")

Idents(sepsis_final) <- "clusters_0.5"

sepsis_final <- RenameIdents(sepsis_final, annotation)

sepsis_final@active.ident

sepsis_final <- StashIdent(sepsis_final, save.name = "labels")

p <- DimPlot(sepsis_final, reduction = "umap.harmony", group.by = c("labels"), pt.size = 2)+
  scale_colour_brewer(palette = "Set1")
p

ggsave(plot = p, filename = "sepsis.umap.png", dpi = 500,
       height = 5, width = 7, path = "../plots/")

p <- DotPlot(sepsis_final, features = c("CD3D", "CD8A", "CD4", "FOXP3", "TRDV2", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                   "CD1C", "IRF7", "C1QA", "CLEC10A", "CD68", "CSF3R",
                                   "MKI67", "ETS2"), group.by = "labels", scale.by = "size",
        dot.min = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)
p

ggsave(plot = p, filename = "annotation.dotplot.png", dpi = 500,
       height = 5, width = 12, path = "../plots/", bg="white")

prop.table(table(sepsis_final$source, sepsis_final$labels), margin = 1)*100

p <- Proportion_Plot(sepsis_final, group_by_var = "labels", split.by = "source")+scale_fill_brewer(palette = "Set1")
p

ggsave(plot = p, filename = "proportions.by.source.png", dpi = 500,
       height = 5, width = 6, path = "../plots/")

Proportion_Plot(sepsis_final, group_by_var = "labels", split.by = "condition2")


saveRDS(sepsis_final, file = "sepsis_final.rds", compress = F)
table(sepsis_final$sepsis_state)

sepsis_final <- readRDS("sepsis_final.rds")

#RIP gene expression============================================================================================
rips <- c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE")
p <- DotPlot(sepsis_final, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE"),
        group.by = "labels", scale.by = "size", dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(plot = p, filename = "RIP.expression.by.celltype.png", dpi = 500,
       height = 5, width = 7, path = "../plots/", bg="white")

DotPlot(sepsis_final, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE"),
        group.by = "condition2", scale.by = "size", dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)

rip_genes <- c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE")

sepsis_final <- AddModuleScore(sepsis_final, features = list(rip_genes),
               name = "RIP_genes", ctrl = 100)


sepsis_final@meta.data <- sepsis_final@meta.data[, !(names(sepsis_final@meta.data) %in% c("RIP_genes1", "RIP_genes2",
                                                                                          "RIP_genes3", "RIP_genes4",
                                                                                          "RIP_genes5", "RIP_genes6",
                                                                                          "RIP_genes7", "RIP_genes8",
                                                                                          "RIP_genes9"))]

DotPlot(sepsis_final, features = "RIP_genes1",
        group.by = "labels", scale.by = "size", dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)

VlnPlot(sepsis_final, features = "RIP_genes1", pt.size = 0, stack = F, group.by = "condition2")

FeaturePlot(sepsis_final, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", 
                                       "JDP2", "MYC", "CEBPE", "RIP_genes1"), 
            reduction = "umap.harmony", ncol = 5, 
            pt.size = 0.5, order = T, raster = F)&
  scale_color_viridis(option = "F", direction = -1)

table(sepsis_final$donor_id, sepsis_final$condition2)

#attempting to bin sepsis into mild vs severe

severe <- c("Darden_S_Late_2", "Darden_S_Late_3", "Darden_S_Late_4",
            "Qiu_NS1_T0", "Qiu_NS1_T6", "Qiu_NS2_T0", "Qiu_NS2_T6",
            "P633", "P636", "P640", "P662", "P669", "P670", "P671", "P672")
mild <- c("Darden_S_Late_1", "Qiu_S1_T0", "Qiu_S1_T6", "Qiu_S2_T0", "Qiu_S2_T6", "Qiu_S3_T0", "Qiu_S3_T6",
          "E1", "E12", "E16", "E25")

sepsis_final$status <- NA

# Assign 'severe' and 'mild' based on donor_id
sepsis_final$status[sepsis_final$donor_id %in% severe] <- "severe"
sepsis_final$status[sepsis_final$donor_id %in% mild]   <- "mild"
# For the rest, use 'condition' (i.e., "HC" or "other")
sepsis_final$status[is.na(sepsis_final$status)] <- sepsis_final$condition2[is.na(sepsis_final$status)]

table(sepsis_final$status)

Idents(sepsis_final) <- "status"
sepsis_only <- subset(sepsis_final, idents = c("HC", "mild", "severe"))

sepsis_only$label_status <- paste0(sepsis_only$labels, "-", sepsis_only$status, "-")
sepsis_only$source_status <- paste0(sepsis_only$source, "-", sepsis_only$status, "-")
sepsis_only$label_condition <- paste0(sepsis_only$labels, "-", sepsis_only$condition2)

dp <- DotPlot(sepsis_only, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE", "RIP_genes1"),
        group.by = "label_condition", scale.by = "size", dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)
dp

plotdata <- dp$data

plotdata <- separate(plotdata, col = id, into = c("label", "condition"), sep = "-")

p <- ggplot(plotdata, aes(x = condition, y = features.plot, colour = avg.exp.scaled, size = pct.exp))+geom_point()+
  theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~label, ncol = 5)+
  scale_y_discrete(limits = rev)+
  scale_colour_viridis(option = "F", direction = -1)+
  theme(axis.text = element_text(size = 12, colour = "black"), strip.text = element_text(size = 12, colour = "black",
                                                                                         face = "bold"))+
  ylab(NULL)+xlab("Disease status")+theme(axis.title = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust =1))
p

ggsave(plot = p, filename = "RIPs.dotplot.png", dpi = 500,
       height = 5, width = 8, path = "../plots/")


Idents(sepsis_only) <- "labels"
myeloid <- subset(sepsis_only, idents = c("CD14+mono", "Mac/mono"))
myeloid$cell_id <- colnames(myeloid)

myeloid$label_status

RidgePlot(myeloid, features = rips, group.by = "label_status", stack = F)&NoLegend()

dp <- DotPlot(myeloid, features = c("ETS2", "MITF", "TP53", "TFE3", "TFEB", "TFEC", "JDP2", "MYC", "CEBPE", "RIP_genes1"),
              group.by = "label_condition", scale.by = "size", dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)
dp

plotdata <- dp$data

plotdata <- separate(plotdata, col = id, into = c("label", "condition"), sep = "-")

p <- ggplot(plotdata, aes(x = condition, y = features.plot, colour = avg.exp.scaled, size = pct.exp))+geom_point()+
  theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~label, ncol = 5)+
  scale_y_discrete(limits = rev)+
  scale_colour_viridis(option = "F", direction = -1)+
  theme(axis.text = element_text(size = 12, colour = "black"), strip.text = element_text(size = 12, colour = "black",
                                                                                         face = "bold"))+
  ylab(NULL)+xlab("Disease status")+theme(axis.title = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust =1))
p

ggsave(plot = p, filename = "RIPs.dotplot.myeloid.png", dpi = 500,
       height = 5, width = 6, path = "../plots/")

#here use either condition2 (for comparing HCv sepsis) or status for severity
avg.exp <- AverageExpression(myeloid, features = rips, group.by = "labels", 
                             add.ident = c("donor_id", "status", "source"))

avg.exp <- as.data.frame(avg.exp$RNA)

avg.exp$gene <- rownames(avg.exp)

avg.exp <- pivot_longer(avg.exp, names_to = "group", values_to = "exp", cols = 1:102)

avg.exp <- separate(avg.exp, col = group, into = c("label", "patient", "status", "source"), sep = "_")

p <- ggplot(filter(avg.exp, gene %in% "JDP2"), aes(x = status, y = exp, fill = status))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+facet_wrap(~label+source, scales = "free", ncol = 3)+geom_jitter(width= 0.1)+
  scale_fill_manual(values = c("grey50", "mediumseagreen", "darkorchid2"))+xlab(NULL)+ylab("ETS2 expression")+
  theme(axis.text = element_text(size = 12, colour = "black"))

p

ggsave(plot = p, filename = "JDP2.by.severity.png", dpi = 500,
       height = 5, width = 7, path = "../plots/")

#IDENTIFYING GENES THAT CAN BE USED AS AN INDICATOR OF DISEASE SEVERITY==================================================================================
FeaturePlot(myeloid, features = c("HLA-DRA"), split.by = ("condition"), reduction = "umap.harmony")
VlnPlot(myeloid, features = "HLA-DRA", pt.size = 0, group.by = "donor_id", split.by = "labels")

#add more specific disease info to metadata
#attempting to bin sepsis into mild vs severe
non_survival <- c("Qiu_NS1_T0", "Qiu_NS1_T6", "Qiu_NS2_T0", "Qiu_NS2_T6")
no_icu <- c("E1", "E12", "E16", "E25", "Darden_S_Late_1")
icu <- c("P633", "P636", "P640", "P662", "P669", "P670", "P671", "P672", "Darden_S_Late_2", "Darden_S_Late_3", "Darden_S_Late_4",
                "Qiu_S1_T0", "Qiu_S1_T6", "Qiu_S2_T0", "Qiu_S2_T6", "Qiu_S3_T0", "Qiu_S3_T6")

sepsis_final$sepsis_state <- NA

# Assign 'severe' and 'mild' based on donor_id
sepsis_final$sepsis_state[sepsis_final$donor_id %in% non_survival] <- "non_survival"
sepsis_final$sepsis_state[sepsis_final$donor_id %in% no_icu] <- "no_icu"
sepsis_final$sepsis_state[sepsis_final$donor_id %in% icu] <- "icu"
# For the rest, use 'condition' (i.e., "HC" or "other")
sepsis_final$sepsis_state[is.na(sepsis_final$sepsis_state)] <- sepsis_final$condition2[is.na(sepsis_final$sepsis_state)]

table(sepsis_final$sepsis_state)

Idents(sepsis_final) <- "sepsis_state"

sepsis_only <- subset(sepsis_final, idents =c("non_survival", "no_icu", "icu", "HC"))

Idents(sepsis_only) <- "labels"
myeloid <- subset(sepsis_only, idents = c("CD14+mono", "Mac/mono"))

DotPlot(myeloid, features = c("HLA-DRA", "HLA-DRB1"), group.by = "sepsis_state")

#genes obtained from prior knowledge or from specific papers
#https://pmc.ncbi.nlm.nih.gov/articles/PMC9271908/
#https://pmc.ncbi.nlm.nih.gov/articles/PMC9968838/
#https://www.mdpi.com/2077-0383/9/1/127
#https://www.sciencedirect.com/science/article/pii/S107476131500076X <- immunity

input <- AverageExpression(myeloid, features = c("LTB",
                                                 "IL1R2", "S100A12", "HIF1A", "CD86", "HLA-DRA", "HLA-DPB1", "ALOX5AP", "ARHGEF10L",
                                                 "TGFB1", "CXCL8", "RETN", "CD36",
                                                 "IRAK3", "ETS2", "JDP2"), group.by = "donor_id", 
                           add.ident = c("source", "donor_id", "sepsis_state"))
input <- as.data.frame(input$RNA)
input$gene <- rownames(input)
input <- pivot_longer(input, names_to = "group", values_to = "exp", 1:102)
input <- separate(input, col = group, into = c("label", "source", "donor", "sepsis_state"), sep = "_")

input$sepsis_state <- factor(input$sepsis_state, levels = c("HC", "no-icu", "icu", "non-survival"),
                             labels = c("HC", "No ICU", "ICU", "Non-survival"))

#get average expression of each gene for the HCs in each source to use a reference to calculate relative expression

hc_means <- input %>%
  filter(sepsis_state == "HC") %>%
  group_by(gene, source) %>%
  summarise(hc_mean = mean(exp), .groups = "drop")

input_norm <- input %>%
  left_join(hc_means, by = c("gene", "source")) %>%
  mutate(relative_exp = (exp-hc_mean))


p <- ggplot(input_norm, aes(x = sepsis_state, y = relative_exp, fill = sepsis_state))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  scale_fill_manual(values = c("grey40", "indianred2", "steelblue2", "darkorchid2"))+facet_wrap(~gene, scales = "free")
p

ggsave(plot = p, filename = "sepsis_severity_markers.png", dpi = 500,
       height = 6, width = 8, path = "../plots/")

#stats
run_tukey <- function(df) {
  aov_model <- aov(exp ~ sepsis_state, data = df)
  tukey_result <- TukeyHSD(aov_model, "sepsis_state")
  # Convert to data frame and add gene name
  tidy_result <- as.data.frame(tukey_result$sepsis_state)
  tidy_result$comparison <- rownames(tidy_result)
  tidy_result$gene <- unique(df$gene)
  return(tidy_result)
}

tukey_results <- input %>%
  group_by(gene) %>%
  group_split() %>%
  map_df(run_tukey)

tukey_results <- tukey_results %>%
  rename(
    diff = diff,
    lwr = lwr,
    upr = upr,
    p_adj = `p adj`
  ) %>%
  select(gene, comparison, diff, lwr, upr, p_adj)

tukey_results <- filter(tukey_results, p_adj < 0.05)

#doing module score (sepsis severity) and then correct for batch by deducting average healthy control score
# Run AddModuleScore
myeloid <- AddModuleScore(myeloid, features = list(c("HLA-DPB1", "HLA-DRA", "CD86", "TGFB1")), name = "negative_correlation")

# Extract metadata
scores_df <- myeloid@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, source, sepsis_state, negative_correlation1)

# Compute mean module score for healthy controls per batch
hc_means <- scores_df %>%
  filter(sepsis_state == "HC") %>%
  group_by(source) %>%
  summarise(hc_mean = mean(negative_correlation1, na.rm = TRUE))

# Join and subtract
scores_df <- left_join(scores_df, hc_means, by = "source") %>%
  mutate(relative_score = negative_correlation1 - hc_mean)

# Add back to Seurat object
myeloid$relative_score <- scores_df$relative_score

myeloid$sepsis_state <- factor(myeloid$sepsis_state, levels = c("HC", "no_icu", "icu", "non_survival"))

VlnPlot(myeloid, features = "relative_score", group.by = "sepsis_state", pt.size = 0)+
  scale_fill_manual(values = c("grey50", "indianred2", "steelblue2", "darkorchid2"))

length(unique(myeloid$donor_id))

DimPlot(myeloid, reduction = "umap.harmony", group.by = "source")

myeloid$sepsis_state_sample <- paste0(myeloid$sepsis_state, "_", myeloid$donor_id)

colors_use <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

p <- FeaturePlot_scCustom(
  myeloid, colors_use = colors_use,
  features = "relative_score",
  reduction = "umap.harmony",
  split.by = "sepsis_state_sample",
  na_cutoff = NA,
  combine = TRUE, num_columns = 7)&theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )& scale_color_gradientn(
    colours = colors_use,
    limits = c(-max(abs(myeloid$relative_score), na.rm = TRUE), max(abs(myeloid$relative_score), na.rm = TRUE)),
    oob = scales::squish,
    na.value = "grey70"
  )
p
ggsave(plot = p, filename = "umaps.donor.module.score.png", dpi = 300,
       height = 20,width = 26, bg = "white", path = "../plots/")

#differential gene expression between each sepsis category
#pseudobulking first
Idents(myeloid) <- "labels" 
classical <- subset(myeloid, idents = "CD14+mono")
non_classical <- subset(myeloid, idents = "Mac/mono")

pseudo_classical <- AggregateExpression(classical, assays = "RNA", return.seurat = T, group.by = c("sepsis_state", "donor_id", "source"))
pseudo_non_classical <- AggregateExpression(non_classical, assays = "RNA", return.seurat = T, group.by = c("sepsis_state", "donor_id", "source"))

Idents(pseudo_classical) <- "sepsis_state"
classical_non_surv_vs_hc <- FindMarkers(pseudo_classical, only.pos = F, ident.1 = "non-survival", ident.2 = "HC",
                                     min.pct = 0, logfc.threshold = 0, test.use = "MAST", latent.vars = "source")
classical_non_surv_vs_hc$sig <- ifelse(classical_non_surv_vs_hc$avg_log2FC > 0.5 & classical_non_surv_vs_hc$p_val_adj < 0.05, "Up",
                                   ifelse(classical_non_surv_vs_hc$avg_log2FC < -0.5 & classical_non_surv_vs_hc$p_val_adj < 0.05, "Down",
                                   "NS"))
classical_non_surv_vs_hc$sig <- factor(classical_non_surv_vs_hc$sig, levels = c("Up", "Down", "NS"))

classical_non_surv_vs_hc_sig <- filter(classical_non_surv_vs_hc, p_val_adj < 0.05)
write.csv(classical_sepsis_degs, file = "classical_monocytes_unfiltered_degs_non_surv_vs_hc.csv")
classical_sepsis_degs_sig <- filter(classical_sepsis_degs, p_val_adj < 0.05)





Idents(pseudo_non_classical) <- "sepsis_state"
non_classical_sepsis_degs <- FindAllMarkers(pseudo_non_classical, only.pos = F, min.pct = 0.05, logfc.threshold = 0, test.use = "MAST", latent.vars = "source")
non_classical_sepsis_degs_sig <- filter(non_classical_sepsis_degs, p_val_adj < 0.05)

#volcanos

ggplot(classical_non_surv_vs_hc, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")+
  scale_colour_manual(values = c("indianred2", "royalblue3", "grey50"))+theme_bw()


#SUBCLUSTERING MYELOID CELLS (also excluding "other" conditions)=============================================================================================
sepsis_final <- readRDS(file = "sepsis_final.rds")
Idents(sepsis_final) <- "condition2"
table(sepsis_final$condition2)
sepsis_only <- subset(sepsis_final, idents = c("HC", "Sepsis"))

Idents(sepsis_only) <- "labels"
myeloid <- subset(sepsis_only, idents = c("CD14+mono", "Mac/mono", "DC"))

rm(sepsis_final, sepsis_only)

myeloid[["RNA"]] <- split(myeloid[["RNA"]], f = myeloid$source)
myeloid

myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid)
myeloid <- ScaleData(myeloid, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA", "source"))
myeloid <- RunPCA(myeloid)

ElbowPlot(myeloid, ndims = 50)

myeloid <- FindNeighbors(myeloid, dims = 1:16, reduction = "pca")
myeloid <- FindClusters(myeloid, resolution = 0.5, cluster.name = "unintegrated_clusters")

myeloid <- RunUMAP(myeloid, dims = 1:16, reduction = "pca", reduction.name = "umap.unintegrated")

saveRDS(myeloid, file = "meyloid.pre-subclustering.rds", compress = F)
myeloid <- readRDS(file = "meyloid.pre-subclustering.rds")

#tring both harmony and scvi (harmony didn't work very well)
myeloid_integrated <- IntegrateLayers(myeloid, method = scVIIntegration,
                                layers = "counts", scale.layer = "scale.data", new.reduction = "integrated.scvi",
                                conda_env = "C://Users/dan94/miniconda3/envs/scrnaseq/")

saveRDS(myeloid_integrated, file = "myeloid_integrated.rds", compress = F)


#run clustering on harmony embedding
myeloid_integrated <- FindNeighbors(myeloid_integrated, dims = 1:16, reduction = "integrated.scvi")

#run a number of resolutions to compare clusters
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.3, cluster.name = "clusters_0.3")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.4, cluster.name = "clusters_0.4")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.5, cluster.name = "clusters_0.5")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.6, cluster.name = "clusters_0.6")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.7, cluster.name = "clusters_0.7")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.8, cluster.name = "clusters_0.8")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 0.9, cluster.name = "clusters_0.9")
myeloid_integrated <- FindClusters(myeloid_integrated, resolution = 1, cluster.name = "clusters_1.0")

myeloid_integrated <- RunUMAP(myeloid_integrated, dims = 1:16, reduction = "integrated.scvi", reduction.name = "umap.scvi")
p <- DimPlot(myeloid_integrated, reduction = "umap.scvi", group.by = c("source", "sepsis_state", "clusters_0.3"))
p

saveRDS(myeloid_integrated, file = "myeloid_integrated.rds", compress = F)
myeloid_integrated <- readRDS(file = "myeloid_integrated.rds")

myeloid_integrated <- JoinLayers(myeloid_integrated)
myeloid_integrated
rm(myeloid)

FeaturePlot(myeloid_integrated, features = c("LYZ", "CD14", "FCGR3A", "CD1C", "CD68", 
                                             "IL10", "CD163", "IL1B", "BIRC3", "CLEC10A"), ncol = 3, reduction = "umap.scvi")&
  scale_colour_viridis(option = "F", direction = -1)

p1 <- VlnPlot(myeloid_integrated, features = c("LYZ", "CD14", "FCGR3A", "CD1C", "CD68", 
                                         "FCER1A", "S100A12", 
                                         "CD163", "IL1B", "BIRC3", "CLEC10A", "LGALS2", "ITGAX",
                                         "IL7R", "CD69", "CLU", "RETN", "MT1X", "CD63", "CD3D"),
        pt.size = 0, group.by = "clusters_0.3", stack = T, flip = T)+
  theme(legend.position = "none")
p1


Idents(myeloid_integrated) <- "clusters_0.3"
myeloid.degs.0.3 <- FindAllMarkers(myeloid_integrated, only.pos = T, min.pct = 0.1)
myeloid.degs.0.3 <- filter(myeloid.degs.0.3, p_val_adj < 0.05)

#cluster 4 is contamination (high CD3 and platelet markers)
Idents(myeloid_integrated) <- "clusters_0.3"
myeloid_integrated <- subset(myeloid_integrated, idents = c("0", "1", "2", "3", "5", "6"))

p2 <- DimPlot(myeloid_integrated, group.by = "clusters_0.3", split.by = "source", reduction = "umap.scvi")+
  ggtitle(NULL)
p3 <- DimPlot(myeloid_integrated, group.by = "clusters_0.3", split.by = "sepsis_state", reduction = "umap.scvi")+
  ggtitle(NULL)+theme(legend.position = "none")
(p2/p3) | p1

labels <- c("0" = "LGALS2+classical", "1" = "LGALS2+classical", "2" = "LGALS2-classical", "3" = "Nonclassical",
            "5" = "DC", "6" = "Activated_classical")

myeloid_integrated <- RenameIdents(myeloid_integrated, labels)

table(myeloid_integrated@active.ident)

myeloid_integrated <- StashIdent(myeloid_integrated, save.name = "labels")


p1 <- VlnPlot(myeloid_integrated, features = c("CD14", "FCGR3A", "CD1C", "CD68", 
                                               "FCER1A", "S100A12", 
                                               "CD163", "IL1B", "CLEC10A", "LGALS2",
                                               "IL7R", "CD69", "CLU", "RETN", "CD63"),
              pt.size = 0, group.by = "labels", stack = T, flip = T)+
  theme(legend.position = "none")
p1
p2 <- DimPlot(myeloid_integrated, group.by = "labels", split.by = "source", reduction = "umap.scvi")+
  ggtitle(NULL)
p3 <- DimPlot(myeloid_integrated, group.by = "labels", split.by = "sepsis_state", reduction = "umap.scvi")+
  ggtitle(NULL)+theme(legend.position = "none")
combo <- (p2/p3) | p1
combo
ggsave(plot = combo, filename = "myeloid_subclustered_annotation.png", dpi = 500,
       height = 6, width = 12, path = "../plots/")

myeloid_integrated$sepsis_state <- factor(myeloid_integrated$sepsis_state, 
                                          levels = c("HC", "no_icu", "icu", "non_survival"))
p <- Proportion_Plot(myeloid_integrated, group_by_var = "labels", split.by = "donor_id")+
  scale_fill_manual(values = c("indianred2", "darkred", "steelblue", "goldenrod2", "mediumseagreen"))+
  coord_flip()
p

ggsave(plot = p, filename = "myeloid.proportions.by.sepsis.state.png", dpi = 500,
       height = 5, width = 7, path = "../plots/")
#export to h5ad for TF enrichment in python
table(myeloid_integrated$sepsis_state)
sce <- as.SingleCellExperiment(myeloid_integrated, assay = "RNA")
table(myeloid_integrated$sepsis_state, myeloid_integrated$labels)

assayNames(sce)
assays(sce, withDimnames = FALSE)$counts_layer <- assays(sce)$counts

writeH5AD(sce, "myeloid_integrated.h5ad")
rm(sce)


#from enrichment analysis in python, plotting number of TFs with increased/decreased activity
#focusing on nonclassical monocytes
tf_data <- data.frame(condition = c("No-ICU", "ICU", "Non-survival"),
                      increased = c(14, 7, 2),
                      decreased = c(1, 5, 14))
tf_data <- pivot_longer(tf_data, names_to = "Activity", values_to = "number", 2:3)

tf_data$condition <- factor(tf_data$condition, levels = c("No-ICU", "ICU", "Non-survival"))

p <- ggplot(tf_data, aes(x = condition, y = number, fill = Activity))+geom_col(position = position_dodge())+
  theme_bw()+scale_fill_manual(values = c("steelblue", "darkred"))+
  ggtitle("Nonclassical monocytes transcription factor enrichment comapred to HC")+
  ylab("Number of TFs")+
  xlab(NULL)+theme(axis.text = element_text(size = 14, colour = "black"))
p

ggsave(plot = p, filename = "nonclassical.tf.enrichment.number.png", dpi = 500,
       height = 5, width = 6, path = "../plots/")

#sepsis severity markers again but split by population
Idents(myeloid_integrated) <- "labels"
nonclassical <- subset(myeloid_integrated, idents = "Nonclassical") 
input <- AverageExpression(nonclassical, features = c("LTB",
                                                 "IL1R2", "S100A12", "HIF1A", "CD86", "HLA-DRA", "HLA-DPB1", "ALOX5AP", "ARHGEF10L",
                                                 "TGFB1", "CXCL8", "RETN", "CD36",
                                                 "IRAK3", "ETS2", "JDP2"), group.by = "donor_id", 
                           add.ident = c("source", "donor_id", "sepsis_state"))
input <- as.data.frame(input$RNA)
input$gene <- rownames(input)
input <- pivot_longer(input, names_to = "group", values_to = "exp", 1:50)
input <- separate(input, col = group, into = c("source", "donor", "sepsis_state"), sep = "_")

input$sepsis_state <- factor(input$sepsis_state, levels = c("HC", "no-icu", "icu", "non-survival"),
                             labels = c("HC", "No ICU", "ICU", "Non-survival"))

#get average expression of each gene for the HCs in each source to use a reference to calculate relative expression

hc_means <- input %>%
  filter(sepsis_state == "HC") %>%
  group_by(gene, source) %>%
  summarise(hc_mean = mean(exp), .groups = "drop")

input_norm <- input %>%
  left_join(hc_means, by = c("gene", "source")) %>%
  mutate(relative_exp = (exp-hc_mean))


p <- ggplot(input_norm, aes(x = sepsis_state, y = relative_exp, fill = sepsis_state))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  scale_fill_manual(values = c("grey40", "indianred2", "steelblue2", "darkorchid2"))+facet_wrap(~gene, scales = "free")
p

ggsave(plot = p, filename = "sepsis_severity_markers.png", dpi = 500,
       height = 6, width = 8, path = "../plots/")

#plotting number of DEGs (generated by Deseq2 in python) in each celltype by comparison
test.df <- data.frame(a = c("hi", "bye", "hello"), b = c(3,5,8), c = c("yellow", "red", "blue"))


deg.df <- data.frame(label = c("LGALS2+_classical", "LGALS2+_classical", "LGALS2+_classical", "LGALS2+_classical", 
                     "LGALS2-_classical", "LGALS2-_classical", "LGALS2-_classical", "LGALS2-_classical", 
                     "non_classical", "non_classical", "non_classical", "non_classical"),
                     comparison = c("non_survival_vs_HC", "icu_vs_HC", "no_icu_vs_HC", "icu_vs_no_icu",
                                    "non_survival_vs_HC", "icu_vs_HC", "no_icu_vs_HC", "icu_vs_no_icu",
                                    "non_survival_vs_HC", "icu_vs_HC", "no_icu_vs_HC", "icu_vs_no_icu"),
                     degs = c(0, 417, 6, 0,
                              172, 0, 0, 0,
                              676, 91, 15, 0))

deg.df$comparison <- factor(deg.df$comparison, levels = c("no_icu_vs_HC", "icu_vs_HC", "non_survival_vs_HC", "icu_vs_no_icu"),
                            labels = c("Non-ICU vs HC", "ICU vs HC", "Non-survivor vs HC",
                                       "Non-ICU vs ICU"))

deg.df[deg.df==0] <- NA


p <- ggplot(deg.df, aes(x = label, y = degs, fill = comparison))+geom_col(position = position_dodge())+theme_bw()+
  xlab(NULL)+ylab("Number of DEGs")+scale_fill_manual(values = c("steelblue2", "darkorchid2", "darkred",
                                                                   "mediumseagreen"))
p

ggsave(plot = p, filename = "DEGs.number.png", dpi = 500,
       height = 4, width = 6, path = "../plots/")
#random forest to understand what genes are associated with disease state
#first want pseudobulk level data
Idents(myeloid_integrated) <- "labels"
nonclassical <- subset(myeloid_integrated, idents = "Nonclassical") 

hvgs <- VariableFeatures(nonclassical)
nonclassical <- subset(nonclassical, features = hvgs)

nonclassical$meta <- paste(nonclassical$donor_id, 
                            nonclassical$sepsis_state,
                            nonclassical$source, sep = "|")

table(nonclassical$meta)

input <- AverageExpression(nonclassical, group.by = "meta")
input <- as.data.frame(input$RNA)
input$gene <- rownames(input)
head(colnames(input))
input <- pivot_longer(input, names_to = "meta", values_to = "expression", cols = matches("\\|"))
input <- separate(input, col = "meta", into = c("donor", "condition", "source"),
                  sep = "\\|")

input <- pivot_wider(input, names_from = "gene", values_from = "expression")

#clean gene names up
names(input) <- make.names(names(input), unique = TRUE)

#make sure condition is a factor
input$condition <- as.factor(input$condition)

library(ranger)
set.seed(61)

rf_model <- ranger(
  formula = condition ~ .,
  data = input[, -c(1,3)],  # remove donor_id or other non-predictors
  num.trees = 10000,
  importance = "impurity",
  classification = TRUE,
  probability = TRUE,
)

rf_model

# View variable importance
importance_scores <- rf_model$variable.importance
importance_scores <- sort(importance_scores, decreasing = TRUE)

head(importance_scores, 10)

imp_df <- data.frame(Gene = names(importance_scores),
                     Importance = importance_scores)

ggplot(imp_df[1:40, ], aes(y = reorder(Gene, Importance), x = Importance)) +
  geom_bar(stat = "identity", fill = "tomato") +
  labs(title = "Top 20 Genes Associated with Condition", y = "Gene", x = "Importance Score")+
  theme_bw()+scale_x_continuous(expand = c(0,0))+theme(axis.text = element_text(size = 14, colour = "black"))

#get predictions of each sample
pred <- rf_model$predictions
rownames(pred) <- input$donor

#ANALYSING DATASETS SEPERATELY=================================================================================================================
#based on feedback from nikolic lab meeting, dataset is far too messy (especially with lack of solid metadata)
#therefore if a look at datsets separately and identify common trends this may be better

#also, to try and make this work more focused, i need a hypothesis

#hypothesis: there is a unique transcription factor signature in patients with unresolving sepsis compared to resolving sepsis
#to answer this i have the following aims:
#1. identify and distribution of myeloid cell subsets in resolving and unresolving sepsis
#2. identify the differentially expressed genes in unresolving sepsis
#3. generate a signature of enriched transcription factors in unresolving sepsis
sepsis <- read_rds("sepsis_final.rds")

sepsis$labels <- factor(sepsis$labels, levels = c("CD4+T", "CD8+T", "NK", "B", "CD14+mono", "Mac/mono",
                                                      "DC", "pDC", "Neutrophil"),
                          labels = c("CD4+T", "CD8+T", "NK", "B", "Classical_mono", "Non-classical_mono",
                                     "DC", "pDC", "Neutrophil"))

p1 <- DimPlot(sepsis, group.by = "labels", reduction = "umap.harmony")+scale_colour_manual(values = umap_cols)
p1

ggsave(plot = p, filename = "integrated_sepsis_annotated_umap.png", dpi = 500,
       height = 5, width = 7)

p2 <- DotPlot(sepsis, group.by = "labels", features = c("CD3D", "CD40LG", "CD8A", "NKG7", "GZMB", "CD79A", "IGHD",
                                  "LYZ", "CD14", "FCGR3A", "CD1C", "IRF7", "CSF3R"), 
             scale.by = "size", dot.min = 0.05)+
  scale_colour_viridis(option = "F", direction = -1)+xlab(NULL)+ylab(NULL)+
  scale_y_discrete(limits = rev)
p2

p3 <- DimPlot(sepsis, split.by = "source", reduction = "umap.harmony")+scale_colour_manual(values = umap_cols)
p3


combo <- (p1|p2) / p3
combo

ggsave(plot = combo, filename = "sepsis_integrated_key_umaps_and_dotplot.png",
       dpi = 500, height = 10, width = 16)


p <- Proportion_Plot(sepsis, group_by_var = "labels", split.by = "condition2")+
  scale_fill_manual(values = umap_cols)
p

ggsave(plot = p, filename = "Qiu_data_proportions.png", dpi = 500,
       height = 4, width = 6)


Idents(sepsis) <- "source"
qiu_data <- subset(sepsis, idents = "Qiu")
table(qiu_data$source)

#combining time 0 and time 6 for each patient to increase numbers
qiu_data$patient <- ifelse(grepl("HC1", qiu_data$donor_id), "HC1",
                           ifelse(grepl("HC2", qiu_data$donor_id), "HC2",
                                  ifelse(grepl("NS1", qiu_data$donor_id), "NS1",
                                         ifelse(grepl("NS2", qiu_data$donor_id), "NS2",
                                                ifelse(grepl("_S1", qiu_data$donor_id), "S1",
                                                       ifelse(grepl("_S2", qiu_data$donor_id), "S2",
                                                              ifelse(grepl("_S3", qiu_data$donor_id), "S3",
                                                              "A")))))))
table(qiu_data$patient)

qiu_data$condition <- ifelse(grepl("HC", qiu_data$patient), "HC",
                             ifelse(grepl("NS", qiu_data$patient), "NS", "S"))

table(qiu_data$condition)

qiu_data$labels <- factor(qiu_data$labels, levels = c("CD4+T", "CD8+T", "NK", "B", "CD14+mono", "Mac/mono",
                                                      "DC", "pDC", "Neutrophil"),
                          labels = c("CD4+T", "CD8+T", "NK", "B", "Classical_mono", "Non-classical_mono",
                                     "DC", "pDC", "Neutrophil"))

qiu_data$condition <- factor(qiu_data$condition, levels = c("HC", "S", "NS"))


DimPlot(qiu_data, group.by = "labels", reduction = "umap.harmony")+scale_colour_manual(values = umap_cols)


p <- Proportion_Plot(qiu_data, group_by_var = "labels", split.by = "condition")+
  scale_fill_manual(values = umap_cols)
p

ggsave(plot = p, filename = "Qiu_data_proportions.png", dpi = 500,
       height = 4, width = 6)

#differential expression analysis
qiu_pseudo <- AggregateExpression(qiu_data, assays = "RNA", return.seurat = T, group.by = c("patient", "labels", "condition"))

qiu_pseudo$condition_label <- paste0(qiu_pseudo$condition, "_", qiu_pseudo$labels)
table(qiu_pseudo$condition_label)

Idents(qiu_pseudo) <- "condition_label"

mono_degs <- FindMarkers(object = qiu_pseudo, ident.1 = "NS_Mac/mono", ident.2 = "S_Mac/mono", test.use = "DESeq2")

#deseq2 needs n = 3 in each group so will run on non-pseudobulk just to see what signature i get
qiu_data$condition_label <- paste0(qiu_data$condition, "_", qiu_data$labels)
Idents(qiu_data) <- "condition_label"
table(qiu_data$condition_label)

mono_ns_vs_s_degs <- FindMarkers(object = qiu_data, ident.1 = "NS_CD14+mono", ident.2 = "S_CD14+mono", test.use = "MAST",
                         min.pct = 0.1, logfc.threshold = 0,
                         latent.vars = c("percent.mt"))

mono_ns_vs_hc_degs_sig <- filter(mono_ns_vs_hc_degs, p_val_adj < 0.05)


#attempting to run DEA on all cell subsets
celltypes <- unique(qiu_data$labels)

de_results <- map(celltypes, function(ct) {
  message("Processing: ", ct)
  
  # Subset using base R condition directly
  subset_obj <- subset(qiu_data, subset = labels == ct)
  
  # Check condition balance
  subtype_counts <- table(subset_obj$condition)
  if (length(subtype_counts) < 2 || any(subtype_counts < 10)) {
    message(paste("Skipping", ct, "due to insufficient cells per condition"))
    return(NULL)
  }
  
  # Run DE with MAST
  markers <- FindMarkers(subset_obj,
                         ident.1 = "NS",
                         ident.2 = "S",
                         group.by = "condition",
                         test.use = "MAST",
                         logfc.threshold = 0.1,
                         min.pct = 0.1)
  
  markers$celltype <- ct
  return(markers)
})

# Combine all results
names(de_results) <- celltypes[!map_lgl(de_results, is.null)]
#if some celltypes don't have degs due to lack of numbers then this will mess order above up
#thefore need to manually alter them
names(de_results)[9] <- "DC"

de_results_df <- bind_rows(de_results, .id = "celltype")

de_results_filtered <- map(de_results, function(df) {
  if (is.null(df)) return(NULL)
  dplyr::filter(df, p_val_adj < 0.05)
})
de_results_filtered_df <- bind_rows(de_results_filtered, .id = "celltype")

write.csv(de_results_df, file = "qiu_ns_vs_s_unfiltered_negbinom.csv")
write.csv(de_results_filtered_df, file = "qiu_ns_vs_s_significant_negbinom.csv")


de_results_df <- read.csv("qiu_ns_vs_s_unfiltered.csv")
de_results_filtered_df <- read.csv("qiu_ns_vs_s_significant.csv")


#add a heuristic stat column
#replace 0 values with a tiny number first
de_results_df$p_val_adj[de_results_df$p_val_adj == 0] <- 1e-300
de_results_df$stat <- -log10(de_results_df$p_val_adj) * sign(de_results_df$avg_log2FC)

de_results_filtered_df$p_val_adj[de_results_filtered_df$p_val_adj == 0] <- 1e-300
de_results_filtered_df$stat <- -log10(de_results_filtered_df$p_val_adj) * sign(de_results_filtered_df$avg_log2FC)

#plot number of degs
de_results_filtered_df %>%
  group_by(celltype) %>%
  summarise(n_DEGs = n())

deg_numbers <- data.frame( "Label" = c("CD4+T", "CD8+T", "NK", "B",
                                          "CD14+mono", "Non_classical_mono", "DC", "pDC", "Neutrophils"),
                          "DEGs" = c(1622, 170, 374, 772, 1114, 575, 0,0,0))


deg_numbers$Label <- factor(deg_numbers$Label, levels = unique(deg_numbers$Label))

umap_cols <- brewer.pal(name = "Paired", n = 12)
umap_cols <- c(umap_cols,"grey50", "black", "cyan4")


p <- ggplot(deg_numbers, aes(y = Label, x = DEGs, fill = Label))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = umap_cols)+
  ggtitle("Number of DEGs in non-survival vs survival sepsis")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)
p

ggsave(plot = p, filename = "qiu_data_degs_number.png", dpi = 500,
       height = 5, width = 6)


#only degs above or below a log2 fold change of 2 (so absolute FC of 4) to be more stringent

#plot number of degs
de_results_filtered_df %>%
  group_by(celltype) %>%
  filter(avg_log2FC > 2 | avg_log2FC < -2) %>%
  summarise(n_DEGs = n())

deg_numbers <- data.frame( "Label" = c("CD4+T", "CD8+T", "NK", "B",
                                       "Classical_mono", "Non_classical_mono", "DC", "pDC", "Neutrophils"),
                           "DEGs" = c(86, 28, 100, 241, 273, 222, 0,0,0))



#separated by up and down
deg_numbers <- de_results_filtered_df %>%
  mutate(Direction = case_when(
    avg_log2FC > 2  ~ "Up",
    avg_log2FC < -2 ~ "Down",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Direction)) %>%
  group_by(celltype, Direction) %>%
  summarise(n_DEGs = n(), .groups = "drop")
deg_numbers

deg_numbers$Label <- factor(deg_numbers$celltype, levels = c("CD4+T", "CD8+T", "NK", "B",
                            "CD14+mono", "Non_classical_mono", "DC", "pDC", "Neutrophils"),
                            labels = c("CD4+T", "CD8+T", "NK", "B",
                                       "Classical_mono", "Non_classical_mono", "DC", "pDC", "Neutrophils"))

umap_cols <- brewer.pal(name = "Paired", n = 12)
umap_cols <- c(umap_cols,"grey50", "black", "cyan4")


p <- ggplot(deg_numbers, aes(y = Label, x = n_DEGs, fill = Direction))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen"))+
  ggtitle("Number of DEGs in non-survival vs survival sepsis \n
          (Fold change > |4|)")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)
p

ggsave(plot = p, filename = "qiu_data_degs_number_log2FC_2.png", dpi = 500,
       height = 5, width = 6)

#get top 25 DEGs from each celltype to inspect
top25_degs <- de_results_filtered_df %>%
  group_by(celltype) %>%
  slice_max(order_by = abs(avg_log2FC), n = 25, with_ties = FALSE) %>%
  ungroup()

top25_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", top25_degs$X)


#load back in if needed
mono_ns_vs_hc_degs <- read.csv("../plots/nonclassical_ns_vs_hc_degs.csv")
mono_s_vs_hc_degs <- read.csv("../plots/nonclassical_s_vs_hc_degs.csv")


mono_ns_vs_hc_degs$significance <- "NS"
mono_ns_vs_hc_degs$significance[mono_ns_vs_hc_degs$p_val_adj < 0.05 & mono_ns_vs_hc_degs$avg_log2FC > 0.5] <- "Upregulated"
mono_ns_vs_hc_degs$significance[mono_ns_vs_hc_degs$p_val_adj < 0.05 & mono_ns_vs_hc_degs$avg_log2FC < -0.5] <- "Downregulated"

mono_s_vs_hc_degs$significance <- "NS"
mono_s_vs_hc_degs$significance[mono_s_vs_hc_degs$p_val_adj < 0.05 & mono_s_vs_hc_degs$avg_log2FC > 0.5] <- "Upregulated"
mono_s_vs_hc_degs$significance[mono_s_vs_hc_degs$p_val_adj < 0.05 & mono_s_vs_hc_degs$avg_log2FC < -0.5] <- "Downregulated"


ggplot(mono_s_vs_hc_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significance))+theme_bw()+
  geom_point() +
  scale_color_manual(values = c("NS" = "grey50",
                                "Upregulated" = "darkred", 
                                "Downregulated" = "royalblue2"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")

#testing the overlap of DEGs
#Two lists and genome size
genomesize <- 5504

ns_sig <- filter(mono_ns_vs_hc_degs, significance %in% c("Upregulated", "Downregulated"))
s_sig <- filter(mono_s_vs_hc_degs, significance %in% c("Upregulated", "Downregulated"))

#be careful here - if its a new list made in this script then rownames(data) will work
#if loading in from a previously made DEG list then it'll likely be data$X
go.obj <- newGeneOverlap(ns_sig$X, s_sig$X, genomesize)
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

overlap <- as.data.frame(go.obj@intersection)
ns_only <- as.data.frame(setdiff(go.obj@listA, go.obj@intersection))
names(ns_only)[1] <- "gene"
s_only <- as.data.frame(setdiff(go.obj@listB, go.obj@intersection))

input <- list(Non_survival = rownames(ns_sig), Survival = rownames(s_sig))

venn.diagram(input, filename = "../plots/classical_mono_degs_venn.png", fill = c("red", "blue"),
             resolution = 500, col = "white", cat.pos = c(180,180), cat.cex = 2,
             cex = 2, alpha = 0.5, label.col = "white", fontface = "bold")

#manually viewing and ranking exclusive genes (with all info not just names)
ns_only_top <- filter(ns_sig, X %in% ns_only$gene)

#exporting nonclassical monocytes (just from qiu et al) for decoupler analysis in python
Idents(qiu_data) <- "labels"
to_export <- subset(qiu_data, idents = "Mac/mono")
library(zellkonverter)
sce <- as.SingleCellExperiment(to_export, assay = "RNA")
assayNames(sce)
assays(sce, withDimnames = FALSE)$counts_layer <- assays(sce)$counts
writeH5AD(sce, "qiu_nonclassical.h5ad")
rm(sce)

#or whole dataset
sce <- as.SingleCellExperiment(qiu_data, assay = "RNA")
assayNames(sce)
assays(sce, withDimnames = FALSE)$counts_layer <- assays(sce)$counts
writeH5AD(sce, "qiu_dataset.h5ad")
rm(sce)


#looking at TF expression based on results from decoupler
myeloid <- subset(qiu_data, idents = c("Mac/mono", "CD14+mono"))
VlnPlot(myeloid, features = c("ATF1", "ATF3", "HIF1A", "JDP2"), 
        group.by = "labels", split.by = "condition", pt.size = 0, ncol = 2)+
  theme(legend.position = "right")



#running gsea on degs list (the one with all celltypes) - i'll just subset by celltype as i go
ns_vs_s_degs <- read.csv("qiu_ns_vs_s_unfiltered.csv")
#run GSEA on DEGs

ns_vs_s_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", ns_vs_s_degs$X)

#soi = subset of interest
#this will change as i analyse different populations - be careful to save plots with correct names
soi_degs <- filter(ns_vs_s_degs, celltype %in% "Non_classical_mono")

#volcano first
soi_degs$sig <- ifelse(soi_degs$avg_log2FC > 1 & soi_degs$p_val_adj < 0.05, "Up",
                                 ifelse(soi_degs$avg_log2FC < -1 & soi_degs$p_val_adj < 0.05, "Down",
                                        "NS"))
soi_degs$sig <- factor(soi_degs$sig, levels = c("Up", "Down", "NS"))

# Top 10 by log2FC
top_fc <- soi_degs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

# Top 10 by adjusted p-value
top_pval <- soi_degs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

# Combine and keep only unique genes
top_union <- bind_rows(top_fc, top_pval) %>%
  distinct(gene, .keep_all = TRUE)

bottom_fc <- soi_degs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%  # ascending order
  slice_head(n = 10)

# Get bottom 10 by adjusted p-value (i.e., worst p-values still under 0.05)
bottom_pval <- soi_degs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(p_val_adj)) %>%
  slice_head(n = 10)

# Combine and keep unique genes
bottom_union <- bind_rows(bottom_fc, bottom_pval) %>%
  distinct(gene, .keep_all = TRUE)

# Combine for labeling
label_genes <- bind_rows(top_union, bottom_union)

volcano <- ggplot(soi_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = -1, linetype = "dashed")+
  scale_colour_manual(values = c("darkorchid3", "mediumseagreen", "grey50"))+theme_bw()+
  ggtitle("Non-classical monocytes (Non-survivor vs survivor)")+
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    segment.color = "grey50"
  )

volcano

ggsave(plot = volcano, filename = "non_classical_mono_ns_vs_s.png", dpi = 500,
       height = 5, width = 6)

#back to gsea


soi_genes <- soi_degs$avg_log2FC
names(soi_genes) <- soi_degs$gene


soi_entrez <- bitr(names(soi_genes), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

soi_genes <- soi_genes[soi_entrez$SYMBOL]
names(soi_genes) <- soi_entrez$ENTREZID

#save genes here if needed to then use in comparecluster
non_classical_genes <- soi_genes


gene_list <- list(non_classical = non_classical_genes, classical= classical_genes)
gene_list <- lapply(gene_list, function(x) sort(x, decreasing = TRUE))

cc_results <- compareCluster(
  geneClusters = gene_list,
  fun = "gseGO",
  OrgDb = org.Hs.eg.db
  #orgdb needed for gseGO but not KEGG
  #reactome needs organism =
  #organism = "human"
)


dotplot(cc_results)
cc_results <- clusterProfiler::simplify(cc_results, cutoff=0.7, by="p.adjust",
                                        select_fun=min)

comparison_results <- cc_results@compareClusterResult

# Top 10 enriched (high NES) per cluster
top_enriched <- comparison_results %>%
  group_by(Cluster) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

# Top 10 unenriched (low NES) per cluster
top_unenriched <- comparison_results %>%
  filter(NES < 0) %>%
  group_by(Cluster) %>%
  arrange(NES) %>%
  slice_head(n = 10) %>%
  ungroup()

#or by p value
# Top 20 enriched by FDR per cluster
top_enriched_p <- comparison_results %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)


# Combine enriched and unenriched
top_pathways <- bind_rows(top_enriched, top_unenriched) %>%
  ungroup()

# Get shared pathways (present in both clusters)
shared_pathways <- comparison_results %>%
  group_by(Description) %>%
  filter(n_distinct(Cluster) == 2) %>%
  ungroup()

# Rank shared pathways by average absolute NES (or use max NES depending on preference)
top_shared <- shared_pathways %>%
  group_by(Description) %>%
  summarize(avg_abs_NES = mean(abs(NES)), .groups = "drop") %>%
  arrange(desc(avg_abs_NES)) %>%
  slice_head(n = 10) %>%
  inner_join(shared_pathways, by = "Description")

# Combine all
final_pathways <- bind_rows(top_pathways, top_shared) %>%
  distinct(Description, Cluster, .keep_all = TRUE)

# Identify presence of each pathway across clusters
cluster_counts <- final_pathways %>%
  distinct(Description, Cluster) %>%
  group_by(Description) %>%
  summarize(cluster_combo = paste(sort(unique(Cluster)), collapse = "_"), .groups = "drop")

# Map combo to group label
cluster_counts <- cluster_counts %>%
  mutate(group = case_when(
    cluster_combo == "non_classical_classical" ~ "shared",
    cluster_combo == "non_classical" ~ "non_classical_only",
    cluster_combo == "classical" ~ "classical_only",
    TRUE ~ "other"
  ))

# Merge back into final_pathways
final_pathways <- final_pathways %>%
  left_join(cluster_counts, by = "Description")


p <- ggplot(final_pathways, aes(x = Cluster, y = reorder(Description, NES), 
                                fill = NES, size = -log10(p.adjust)))+
  geom_point(shape = 21)+scale_size_continuous(range = c(3,8))+
  theme_bw()+scale_fill_distiller(palette = "RdBu")+
  ylab(NULL)+xlab(NULL)+
  theme(axis.text = element_text(size = 16, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_grid(group ~ ., space = "free", scales = "free")+
  theme(axis.text.y = element_text(size = 16, colour = "black"))
p

ggsave(plot = p, filename = "GSEA_monocytes_non_survivor_vs_survivial.png", dpi = 500,
       height = 16, width = 10)

top_pathways$direction <- ifelse(top_pathways$NES > 0, "Enriched in non-survivors",
                                   ifelse(top_pathways$NES < 0, "Enriched in survivors", "Other"))

#visualising by bar instead
p <- ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES), fill = direction))+
  geom_col()+theme_bw()+
  facet_wrap(~Cluster, scales = "free_y", ncol = 1)+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen", "grey50"))+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  theme(strip.text = element_text(size = 14, colour = "black"))+
  ylab(NULL)+xlab("Normalised enrichment score")

p

ggsave(plot = p, filename = "GSEA_monocytes_non_surviro_vs_survivor_bars.png", dpi = 500,
       height = 12, width = 12)
#plot by p value
max_abs <- max(abs(top_enriched_p$NES), na.rm = TRUE)


ggplot(top_enriched_p, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), fill = NES))+
         geom_col()+theme_bw()+
  facet_wrap(~Cluster, scales = "free_y", ncol = 1)+
  scale_fill_distiller(palette = "RdBu", limits = c(-max_abs, max_abs))


#get genes that make up the pathways
gene_long <- comparison_results %>%
  dplyr::select(ID, Description, Cluster, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") %>%
  rename(ENTREZID = core_enrichment)

# Step 2: Map Entrez IDs to gene symbols
gene_symbols <- bitr(gene_long$ENTREZID, fromType = "ENTREZID",
                     toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Step 3: Merge back
gene_long_named <- gene_long %>%
  left_join(gene_symbols, by = "ENTREZID")


