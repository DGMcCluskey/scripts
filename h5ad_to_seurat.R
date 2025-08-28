library(zellkonverter)
library(scCustomize)
library(Seurat)


setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Claire_Smith_data")
ali <- readH5AD("ALI_Flu_250626.h5ad")
ali <- readRDS("ALI_Flu_250626_seurat.rds")

DotPlot(ali, features = c("CD3D", "LYZ", "COL1A1"), group.by = "Infection_status")

ali@assays@data$counts <- ali@assays@data$X

VlnPlot(ali, features = c("KRT5"), pt.size = 0, group.by = "Ziegler_majority_voting")


ali_seurat <- as.Seurat(ali)
