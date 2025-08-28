#converting h5 (normally anndata from scanpy workflow) files to seurat objects
library(zellkonverter)
library(Seurat)
poms <- readH5AD("poms.annotated.dan.h5ad")
poms.seurat <- as.Seurat(poms)
rm(poms)