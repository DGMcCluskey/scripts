#save as an anndata object to be used in python
library(zellkonverter)
sce <- as.SingleCellExperiment(myeloid, assay = "RNA")
assayNames(sce)
assays(sce, withDimnames = FALSE)$counts_layer <- assays(sce)$counts

writeH5AD(sce, "myeloid.h5ad")
rm(sce)

#also need to manually extract the nieghbor graph so that exact clustering can be used with celltypist 
#(and to avoid running pca again which is causing memory trouble in python)
# Extract the shared nearest neighbor graph (usually used for clustering)
graph <- sepsis_final[["RNA_snn"]]

# Convert to a matrix or sparse matrix
graph_matrix <- as(graph, "dgCMatrix")

# Save to file (Matrix Market format is compatible with Python)
Matrix::writeMM(graph_matrix, file = "snn_graph.mtx")

#also save cell names
write.csv(colnames(sepsis_final), file = "cell_names.csv", row.names = FALSE)

