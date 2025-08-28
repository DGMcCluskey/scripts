#converting h5 (normally anndata from scanpy workflow) files to seurat objects
setwd("/Users/dan94/rips_project/")
library(zellkonverter)
library(Seurat)
data <- readH5AD("GSE174188_CLUES1_adjusted.h5ad")
data.seurat <- as.Seurat(data)
rm(data)