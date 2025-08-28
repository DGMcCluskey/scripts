setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Cutaneous lupus public data/GSE186476_RAW/")

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
library(DoubletFinder)

files <- list.files(path = getwd())
#files <- files[1:21]
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
cle <- merge(
  x = data[[1]],
  y = data[-1],
  add.cell.ids = sample_names,
  project = "CLE"
)

rm(data)

table(cle$orig.ident)
cle$sample <- cle$orig.ident
table(cle$sample)

cle

cle <- JoinLayers(cle)
cle

#filtering and processing

cle[["percent.mt"]] <- PercentageFeatureSet(cle, pattern = "^MT-")
VlnPlot(cle, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
cle <- subset(cle, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)

cle <- NormalizeData(cle)
cle <- FindVariableFeatures(cle)
cle <- ScaleData(cle, vars.to.regress = c("nCount_RNA", "percent.mt"))
cle <- RunPCA(cle)

ElbowPlot(cle, ndims = 30)

cle <- FindNeighbors(cle, dims = 1:14, reduction = "pca")

cle <- RunUMAP(cle, dims = 1:14, reduction = "pca", seed.use =61)

cle <- FindClusters(cle, resolution = 0.3, cluster.name = "resolution_0.3")
cle <- FindClusters(cle, resolution = 0.4, cluster.name = "resolution_0.4")
cle <- FindClusters(cle, resolution = 0.5, cluster.name = "resolution_0.5")
cle <- FindClusters(cle, resolution = 0.6, cluster.name = "resolution_0.6")
cle <- FindClusters(cle, resolution = 0.7, cluster.name = "resolution_0.7")
cle <- FindClusters(cle, resolution = 0.8, cluster.name = "resolution_0.8")
cle <- FindClusters(cle, resolution = 0.9, cluster.name = "resolution_0.9")
cle <- FindClusters(cle, resolution = 1.0, cluster.name = "resolution_1.0")

p <- DimPlot(cle, group.by = c("resolution_0.5"))
p

#add metadata
table(cle$sample)

cle$condition <- ifelse(grepl("Healthy", cle$sample), "HC", "CLE")
cle$location <- ifelse(grepl("Healthy", cle$sample), "HC",
                      ifelse(grepl("NonLes", cle$sample), "Non_lesional", "Lesional"))

# Extract everything from "D" onwards (D + number)
cle$patient_id <- sub(".*_(D\\d+)_.*", "\\1", cle$orig.ident)
cle$patient_id <- paste0(cle$condition, "_", cle$patient_id)
table(cle$patient_id)
table(cle$sample)
scle <- c("CLE_D1", "CLE_D2", "CLE_D3", "CLE_D5")
dle <- c("CLE_D4", "CLE_D6", "CLE_D7")
cle$subtype <-NA


# Assign SCLE and DLE based on donor_id
cle$subtype[cle$patient_id %in% scle] <- "SCLE"
cle$subtype[cle$patient_id %in% dle]   <- "DLE"
# For the rest, use 'condition' (i.e., "HC" or "other")
cle$subtype[is.na(cle$subtype)] <- cle$condition[is.na(cle$subtype)]

table(cle$subtype)


p <- DimPlot(cle, group.by = c("condition", "location", "subtype"))
p


#pca at patient level to quickly see if dle and scle differ
Idents(cle) <-"sample"

patient_avg <- AverageExpression(cle, group.by = "sample", return.seurat = T)
patient_avg <- FindVariableFeatures(patient_avg)
patient_avg <- ScaleData(patient_avg)
patient_avg <- RunPCA(patient_avg, features = VariableFeatures(patient_avg), npcs = 27)

patient_metadata <- cle@meta.data %>%
  group_by(sample) %>%
  slice(1) %>%  # Take first row per patient
  dplyr::select(sample, condition, patient_id, subtype, location) %>%  # Add other patient-level variables as needed
  as.data.frame()

rownames(patient_metadata) <- patient_metadata$sample
patient_avg <- AddMetaData(patient_avg, patient_metadata)

DimPlot(patient_avg, reduction = "pca", group.by = c("condition", "subtype",
                                                     "patient_id", "location", pt.size = 2), ncol = 2)

#back to annotation
clustree(cle, prefix = "resolution_")

DimPlot(cle, group.by = "label", label = T)+NoLegend()
table(cle$resolution_0.8)
FeaturePlot(cle, features = c("CD3D", "CD8A", "NKG7", "CD79A", "IGHD", "LYZ", "CD14", "FCGR3A", "IRF7", "PPBP", "HBB"))&
  scale_colour_viridis(option = "F", direction = -1)

DotPlot(cle, features = c("PTPRC", "CD3D", "CD8A", "CD40LG",
                          "MS4A1", "CD79A", "LYZ", "CLEC9A", "CLEC10A",
                          "IRF7", "JCHAIN",
                          "CD68", "CD163", "CD14", "FCGR3A",
                          "TPSB2", "TPSAB1",
                          "COL3A1", "DCN", "ACTA2", "RGS5",
                          "CDH5", "PROX1", "SELE",
                          "DCD", "PIP", "KRT14", "KRT5", "KRT17", "KRT15",
                          "CDH19", "SOX10"), 
        group.by = "resolution_0.8", scale.by = "size", scale = 10)+
  scale_colour_viridis(option = "F", direction = -1)+theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))

saveRDS(cle, file = "bili_etal_cle_seurat.rds", compress = F)

#pull out immune cells for subclustering and analysis
Idents(cle) <- "resolution_0.8"
immune_cle <- subset(cle, idents = c("0", "10", "14", "22"))

immune_cle <- NormalizeData(immune_cle)
immune_cle <- FindVariableFeatures(immune_cle)
immune_cle <- ScaleData(immune_cle, vars.to.regress = c("nCount_RNA", "percent.mt"))
immune_cle <- RunPCA(immune_cle)

ElbowPlot(immune_cle, ndims = 30)

immune_cle <- FindNeighbors(immune_cle, dims = 1:15, reduction = "pca")

immune_cle <- RunUMAP(immune_cle, dims = 1:15, reduction = "pca", seed.use =61)

immune_cle <- FindClusters(immune_cle, resolution = 0.3, cluster.name = "resolution_0.3")
immune_cle <- FindClusters(immune_cle, resolution = 0.4, cluster.name = "resolution_0.4")
immune_cle <- FindClusters(immune_cle, resolution = 0.5, cluster.name = "resolution_0.5")
immune_cle <- FindClusters(immune_cle, resolution = 0.6, cluster.name = "resolution_0.6")
immune_cle <- FindClusters(immune_cle, resolution = 0.7, cluster.name = "resolution_0.7")
immune_cle <- FindClusters(immune_cle, resolution = 0.8, cluster.name = "resolution_0.8")
immune_cle <- FindClusters(immune_cle, resolution = 0.9, cluster.name = "resolution_0.9")
immune_cle <- FindClusters(immune_cle, resolution = 1.0, cluster.name = "resolution_1.0")

p <- DimPlot(immune_cle, group.by = c("resolution_0.5"))
p

#doublet removal using DoubletFinder
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_cle <- paramSweep(immune_cle, PCs = 1:15, sct = FALSE)
sweep.stats_cle <- summarizeSweep(sweep.res.list_cle, GT = FALSE)
bcmvn_cle <- find.pK(sweep.stats_cle)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- immune_cle@meta.data$labels
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(immune_cle@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
immune_cle <- doubletFinder(immune_cle, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

DimPlot(immune_cle, group.by = "DF.classifications_0.25_0.09_624")

FeaturePlot(immune_cle, features = c("CD3D", "CD8A", "NKG7", "CD79A", "IGHD", "LYZ", "CD14", "FCGR3A", "IRF7"))

DotPlot(cle_immune, features = c("CD3D", "CD8A", "CD40LG", "CCR7", "SELL", "LEF1", 
                                 "KLF2", "FOXP1", "LMNA", "ANXA1", "ATP1B1", "GEM", "SELPLG", "CD83", "CXCL13",
                                 "NKG7", "GZMB", "GZMK", "GNLY", "FCER1G",
                                 "CCL5",
                                 "FOXP3", "KLRB1", "IFIT1", "STAT1",
                                 "IFIT3",
                                 "CD79A", "IGHD", "LYZ", "CD14", "FCGR3A", "IL3RA", "C1QA", "APOE", "CD207",
                                 "CD1C", "IRF7", "TPSB2"), 
        group.by = "big_labels", scale.by = "size", scale = 10, dot.min = 0.1)+
  scale_colour_viridis(option = "F", direction = -1)+theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))

annotation <- c("0" = "CD4+T_em", "1" = "CD4+T_cm", "2" = "CD8+T_em",
                "3" = "Treg", "4" = "CD4+T_em", "5" = "Tcm",
                "6" = "CD8+T_em", "7" = "CD4+T_em", "8" = "CD8+T_cm",
                "9" = "Mast", "10" = "cDC", "11" = "CD8+T_IFN",
                "12" = "Mono", "13" = "pDC", "14" = "Mono_IFN", "15" = "CD8+T_em", 
                "16" = "LC", "17" = "CD8+T_em", "18" = "Doublets", "19" = "CD4+T_em", "20" = "B")

Idents(immune_cle) <- "resolution_1.0"

immune_cle <- RenameIdents(immune_cle, annotation)
immune_cle@active.ident

immune_cle <- StashIdent(immune_cle, save.name = "labels")

table(immune_cle$labels)

table(immune_cle$labels, immune_cle$condition)

saveRDS(immune_cle, file = "bili_etal_cle_immune_only.rds", compress = F)
cle_immune <- readRDS("bili_etal_cle_immune_only.rds")

umap_cols <- brewer.pal(name = "Paired", n = 12)
umap_cols <- c(umap_cols,"grey50", "black", "cyan4")
p <- DimPlot(cle_immune, group.by = "labels", pt.size = 0.5)+scale_colour_manual(values = umap_cols)
p

ggsave(plot = p, filename = "CLE_UMAP_annotated.png", dpi = 500,
       height = 5, width = 8)

p <- DimPlot(immune_cle, group.by = c("big_labels", "subtype", "location") ,
             pt.size = 0.5)&scale_colour_manual(values = umap_cols)
p 

ggsave(plot = p, filename = "CLE_UMAPs_combined.png", dpi = 500,
       height = 5, width = 16)


immune_cle$subtype <- factor(immune_cle$subtype, levels = c("HC", "SCLE", "DLE"))


Idents(cle_immune) <- "big_labels"
cle_immune <- subset(cle_immune, idents = "Doublets", invert = T)
Idents(cle_immune) <- "location"
nonlesional <- subset(cle_immune, idents = "Non_lesional")
lesional <- subset(cle_immune, idents = "Lesional")

p1 <- Proportion_Plot(lesional, group_by_var = "big_labels", 
                split.by = "subtype")+scale_fill_manual(values = umap_cols)+
  ggtitle("Lesional skin")
p2 <- Proportion_Plot(nonlesional, group_by_var = "big_labels", 
                      split.by = "subtype")+scale_fill_manual(values = umap_cols)+
  ggtitle("Non-lesional skin")

p1+theme(legend.position = "none")+
  p2+ylab(NULL)+plot_layout(ncol = 2)

ggsave(plot = prop.data, filename = "proportion.bars.png", dpi = 500,
       height = 6, width = 6)



prop.data <- prop.data$data

patient_metadata <- cle_immune@meta.data %>%
  group_by(sample) %>%
  summarise(
    condition = first(condition),
    # Add any other metadata you want
    subtype = first(subtype),           # if available
    location = first(location),
    patient_id = first(patient_id),# if available
    .groups = 'drop'
  )

names(prop.data)[2] <- "sample"

prop.data <- prop.data %>%
  left_join(patient_metadata, by = "sample")

prop.data$group <- paste0(prop.data$subtype, "_", prop.data$location)

prop.data$group <- factor(prop.data$group, 
                          levels = c("HC_HC", "SCLE_Non_lesional",
                                     "DLE_Non_lesional",
                                     "SCLE_Lesional",
                                     "DLE_Lesional"))

prop_check <- prop.data %>%
  group_by(sample) %>%
  summarise(total_prop = sum(value))

print(prop_check)

ggplot(prop.data, aes(x = group, y = value, fill = group))+
  scale_fill_manual(values = c("grey50", "indianred2", "steelblue2", 
                               "darkred", "cyan4"))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Cluster, scales = "free_y", ncol = 5)+
  theme_bw()+xlab(NULL)

anova_res <- aov(value ~ Cluster * group, data = prop.data)
tukey_res <- TukeyHSD(anova_res)
tukey_df <- data.frame(tukey_res$`Cluster:group`)
tukey_df$comparison <- rownames(tukey_df)
#tukey_df <- filter(tukey_df, p.adj < 0.05)

# Extract cluster from each side of the comparison
tukey_df <- tukey_df %>%
  mutate(
    # Split the comparison and extract clusters
    comparison_split = strsplit(comparison, "-"),
    cluster1 = sapply(comparison_split, function(x) sub(":.*", "", x[1])),
    cluster2 = sapply(comparison_split, function(x) sub(":.*", "", x[2]))
  ) %>%
  # Keep only rows where clusters match
  filter(cluster1 == cluster2)

#adding broader label names
big_labels <- c("CD4+T_em" = "CD4+T", "CD4+T_cm" = "CD4+T", "CD8+T_em" = "CD8+T", "CD8+T_cm" = "CD8+T",
                "Treg" = "Treg", "B" = "B", "pDC" = "pDC", "Mono" = "Mono", "Mono_IFN" = "Mono",
                "CD8+T_IFN" = "CD8+T", "Mast" = "Mast", "Doublets" = "Doublets", "LC" = "DC", "cDC" = "DC")

Idents(immune_cle) <- "labels"

immune_cle <- RenameIdents(immune_cle, big_labels)
immune_cle@active.ident

immune_cle <- StashIdent(immune_cle, save.name = "big_labels")
table(immune_cle$big_labels)

immune_cle$biglabel_subtype <- paste0(immune_cle$big_labels, "_", immune_cle$subtype)
table(immune_cle$biglabel_subtype)

#account for ribosomal genes
ribosomal_genes <- grep("^RPS|^RPL", rownames(immune_cle), value = TRUE)
immune_cle[["percent.ribo"]] <- PercentageFeatureSet(immune_cle, pattern = "^RP[SL]")

table(immune_cle$location)

Idents(immune_cle) <- "location"
immune_cle_lesional <- subset(immune_cle, idents = c("Lesional"))
immune_cle_nonlesional <- subset(immune_cle, idents = c("Non_lesional"))
table(immune_cle_nonlesional$subtype)
#not pseudobulkig due to low numbers (at most it will be 4 vs 3 but reality is most populations not represented in 3 patients)
celltypes <- unique(immune_cle_lesional$big_labels)

de_results <- map(celltypes, function(ct) {
  message("Processing: ", ct)
  
  # Subset using base R condition directly
  subset_obj <- subset(immune_cle_lesional, subset = big_labels == ct)
  
  # Check condition balance
  subtype_counts <- table(subset_obj$subtype)
  if (length(subtype_counts) < 2 || any(subtype_counts < 10)) {
    message(paste("Skipping", ct, "due to insufficient cells per condition"))
    return(NULL)
  }
  
  # Run DE with MAST
  markers <- FindMarkers(subset_obj,
                         ident.1 = "SCLE",
                         ident.2 = "DLE",
                         group.by = "subtype",
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
names(de_results)[5] <- "Tcm"
de_results_df <- bind_rows(de_results, .id = "celltype")

de_results_filtered <- map(de_results, function(df) {
  if (is.null(df)) return(NULL)
  dplyr::filter(df, p_val_adj < 0.05)
})
de_results_filtered_df <- bind_rows(de_results_filtered, .id = "celltype")

l_de_results <- de_results_df
l_de_results_filtered <- de_results_filtered_df

write.csv(de_results_df, file = "lesional_SCLE_vs_DLE_unfiltered.csv")
write.csv(de_results_filtered_df, file = "lesional_SCLE_vs_DLE_significant.csv")

#plot number of degs
l_de_results_filtered %>%
  group_by(celltype) %>%
  summarise(n_DEGs = n())

deg_numbers <- data.frame("Location" = c(rep("Non_lesional", 9),
                                         rep("Lesional", 9)),
                          "Label" = rep(c("CD4+T", "CD8+T", "Treg", "Tcm",
                                      "B", "Mono", "DC", "pDC", "Mast"), 2),
                          "DEGs" = c(99, 1231, 159, 82, 0, 302, 22, 1754,
                                     0,
                                     7, 107, 0,0,0,741,109,0,0))


deg_numbers$Label <- factor(deg_numbers$Label, levels = rev(unique(deg_numbers$Label)))

p <- ggplot(deg_numbers, aes(y = Label, x = DEGs, fill = Location))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen"))+
  ggtitle("Number of DEGs in SCLE vs DLE")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))
p

ggsave(plot = p, filename = "DEG_numbers_SCLE_vs_DLE.png", dpi = 500,
       height = 4, width = 6)

#get top 25 DEGs from each celltype to inspect
top25_degs <- non_lesional_degs %>%
  group_by(celltype) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = abs(avg_log2FC), n = 25, with_ties = FALSE) %>%
  ungroup()

top25_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", top25_degs$X)

#next steps and questions
#what pathways are degs involved in
#how are these cell types functional or phenotypically different between
#scle and dle

#run GSEA on DEGs
non_lesional_degs <- read.csv("non_lesional_SCLE_vs_DLE_unfiltered.csv")
lesional_degs <- read.csv("lesional_SCLE_vs_DLE_unfiltered.csv")

non_lesional_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", non_lesional_degs$X)
lesional_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", lesional_degs$X)

mono_non_lesional_degs <- filter(non_lesional_degs, celltype %in% "Mono")
mono_lesional_degs <- filter(lesional_degs, celltype %in% "Mono")

#volcano first
mono_lesional_degs$sig <- ifelse(mono_lesional_degs$avg_log2FC > 0.5 & mono_lesional_degs$p_val_adj < 0.05, "Up",
                                     ifelse(mono_lesional_degs$avg_log2FC < -0.5 & mono_lesional_degs$p_val_adj < 0.05, "Down",
                                            "NS"))
mono_lesional_degs$sig <- factor(mono_lesional_degs$sig, levels = c("Up", "Down", "NS"))
l_volcano <- ggplot(mono_lesional_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")+
  scale_colour_manual(values = c("darkorchid", "mediumseagreen", "grey50"))+theme_bw()+
  ggtitle("DEGs in lesional monocytes (SCLE vs DLE)")

l_volcano

mono_non_lesional_degs$sig <- ifelse(mono_non_lesional_degs$avg_log2FC > 0.5 & mono_non_lesional_degs$p_val_adj < 0.05, "Up",
                                       ifelse(mono_non_lesional_degs$avg_log2FC < -0.5 & mono_non_lesional_degs$p_val_adj < 0.05, "Down",
                                              "NS"))
mono_non_lesional_degs$sig <- factor(mono_non_lesional_degs$sig, levels = c("Up", "Down", "NS"))
nl_volcano <- ggplot(mono_non_lesional_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")+
  scale_colour_manual(values = c("darkorchid", "mediumseagreen", "grey50"))+theme_bw()+
  ggtitle("DEGs in non-lesional monocytes (SCLE vs DLE)")

nl_volcano

combo <- l_volcano+nl_volcano+plot_layout(ncol = 1)
combo

ggsave(plot = combo, filename = "monocyte_volcanos.png", dpi = 500,
       height = 8, width = 6)

#back to gsea


nl_genes <- mono_non_lesional_degs$avg_log2FC
names(nl_genes) <- mono_non_lesional_degs$gene


nl_entrez <- bitr(names(nl_genes), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

nl_genes <- nl_genes[nl_entrez$SYMBOL]
names(nl_genes) <- nl_entrez$ENTREZID


l_genes <- mono_lesional_degs$avg_log2FC
names(l_genes) <- mono_lesional_degs$gene


l_entrez <- bitr(names(l_genes), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

l_genes <- l_genes[l_entrez$SYMBOL]
names(l_genes) <- l_entrez$ENTREZID



gene_list <- list(non_lesional = nl_genes, lesional = l_genes)
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
    cluster_combo == "non_lesional_lesional" ~ "shared",
    cluster_combo == "lesional" ~ "lesional_only",
    cluster_combo == "non_lesional" ~ "non_lesional_only",
    TRUE ~ "other"
  ))

# Merge back into final_pathways
final_pathways <- final_pathways %>%
  left_join(cluster_counts, by = "Description")

library(forcats)
final_pathways <- final_pathways %>%
  group_by(Cluster) %>%
  arrange(desc(NES)) %>%
  mutate(Description_facet = factor(Description, levels = unique(Description))) %>%
  ungroup()


#final_pathways$group <- factor(final_pathways$group, levels = c("lesional_only", "shared", "non_lesional_only"))
p <- ggplot(final_pathways, aes(x = Cluster, y = reorder(Description, NES), 
                                fill = NES, size = -log10(p.adjust)))+
  geom_point(shape = 21)+scale_size_continuous(range = c(3,8))+
  theme_bw()+scale_fill_distiller(palette = "RdBu")+
  ylab(NULL)+xlab(NULL)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(group ~ ., space = "free", scales = "free")+
  theme(axis.text.y = element_text(size = 14, colour = "black"))
p


library(stringr)
final_pathways$Description <- str_wrap(final_pathways$Description, width = 60)

#pull out lesional an non lesional pathways separately to avoid the facet ordering issue
top_pathways$direction <- ifelse(top_pathways$NES > 0, "Up in SCLE", "Up in DLE")


lesional_pathways <- filter(top_pathways, Cluster %in% "lesional")
non_lesional_pathways <- filter(top_pathways, Cluster %in% "non_lesional")

p1 <- ggplot(lesional_pathways, aes(x = NES, y = reorder(Description, NES), fill = direction))+
  geom_col()+theme_bw()+scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+ylab(NULL)+
  ggtitle("Lesional pathway enrichment SCLE vs DLE")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  xlab(NULL)
p1

p2 <- ggplot(non_lesional_pathways, aes(x = NES, y = reorder(Description, NES), fill = direction))+
  geom_col()+theme_bw()+scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+ylab(NULL)+
  ggtitle("Non_lesional pathway enrichment SCLE vs DLE")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  xlab("Normalised enrichment score")
p2

combo <- p1+p2+plot_layout(ncol = 1)
combo

ggsave(plot = combo, filename = "monocyte_pathways_gsea.png", dpi = 500,
       height = 12, width = 12)

#CELLCHAT=========================================================================================
library(CellChat)
#cell-cell interactions
#first using cellchat
cle_immune <- read_rds("bili_etal_cle_immune_only.rds")
Idents(cle_immune) <- "big_labels"
cle_immune <- subset(cle_immune, idents = c("Doublets"), invert = T)
Idents(cle_immune) <- droplevels(Idents(cle_immune))
cle_immune$big_labels <- droplevels(cle_immune$big_labels)
table(cle_immune$big_labels)
Idents(cle_immune) <- "location"
cle_immune$samples <- cle_immune$sample
cle_immune <- subset(cle_immune, idents = "Non_lesional")
Idents(cle_immune) <- "subtype"
DLE <- subset(cle_immune, idents = "DLE")
SCLE <- subset(cle_immune, idents = "SCLE")

DLE.chat <- createCellChat(DLE, group.by = "big_labels", assay = "RNA")
SLE.chat <- createCellChat(SCLE, group.by = "big_labels", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#use of cellchatdb except for non-protein signalling
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
DLE.chat@DB <- CellChatDB.use
SLE.chat@DB <- CellChatDB.use
#HC.chat@DB <- CellChatDB.use

DLE.chat <- subsetData(DLE.chat) # This step is necessary even if using the whole database
SLE.chat <- subsetData(SLE.chat) # This step is necessary even if using the whole database
#HC.chat <- subsetData(HC.chat)


future::plan("multisession", workers = 4) # do parallel
DLE.chat <- identifyOverExpressedGenes(DLE.chat)
DLE.chat <- identifyOverExpressedInteractions(DLE.chat)
SLE.chat <- identifyOverExpressedGenes(SLE.chat)
SLE.chat <- identifyOverExpressedInteractions(SLE.chat)
#HC.chat <- identifyOverExpressedGenes(HC.chat)
#HC.chat <- identifyOverExpressedInteractions(HC.chat)

DLE.chat <- subsetData(DLE.chat)
SLE.chat <- subsetData(SLE.chat)

#DLE.chat <- projectData(DLE.chat, PPI.human)
DLE.chat <- computeCommunProb(DLE.chat, type = "triMean", population.size = F)
SLE.chat <- computeCommunProb(SLE.chat, type = "triMean", population.size = F)
#HC.chat <- computeCommunProb(HC.chat, type = "triMean")


DLE.chat <- filterCommunication(DLE.chat, min.cells = 10)
SLE.chat <- filterCommunication(SLE.chat, min.cells = 10)
#HCE.chat <- filterCommunication(HC.chat, min.cells = 10)

DLE.chat <- computeCommunProbPathway(DLE.chat)
SLE.chat <- computeCommunProbPathway(SLE.chat)
#HC.chat <- computeCommunProbPathway(HC.chat)

DLE.chat <- aggregateNet(DLE.chat)
SLE.chat <- aggregateNet(SLE.chat)
#HC.chat <- aggregateNet(HC.chat)

options(future.seed = TRUE)

DLE.chat <- netAnalysis_computeCentrality(DLE.chat, slot.name = "netP")
SLE.chat <- netAnalysis_computeCentrality(SLE.chat, slot.name = "netP")
#HC.chat <- netAnalysis_computeCentrality(HC.chat, slot.name = "netP")

#saving each cellchat object so i dont have to load again
saveRDS(DLE.chat, file = "DLE.cellchat.rds")
saveRDS(SLE.chat, file = "SCLE.cellchat.rds")
#saveRDS(HC.chat, file = "HC.cellchat.rds")

object.list <- list(SCLE = SLE.chat, DLE = DLE.chat)
cle.chat <- mergeCellChat(object.list, add.names = names(object.list))


netVisual_diffInteraction(cle.chat, comparison = c(1,2), vertex.size.max = 10, edge.width.max = 1,
                          label.edge = F)

p1 <- compareInteractions(cle.chat, show.legend = F, group = c("SCLE", "DLE"))+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+#scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 16, colour = "black"))
p2 <- compareInteractions(cle.chat, show.legend = F, group = c("SCLE", "DLE"), measure = "weight")+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+#scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 16, colour = "black"))

p <- p1+p2

p

ggsave(plot = p, filename = "SCLE_DLE_interactions.number.comparison.png",
       dpi = 500,
       height = 5, width = 6)

netVisual_heatmap(cle.chat, comparison = c(1,2))
netVisual_heatmap(cle.chat, measure = "weight")

rankNet(cle.chat, mode = "comparison", comparison = c(1,2), measure = "weight", stacked = F, do.stat = T)
pathways <- rankNet(cle.chat, mode = "comparison", comparison = c(1,2), measure = "weight", stacked = T, do.stat = T)
pathways
pathways <- pathways$data
pathways <- filter(pathways, pvalues < 0.05)

diff.pathways <- unique(pathways$name)

p <- ggplot(pathways, aes(x = contribution.scaled, y = name, fill = group))+
  geom_col(position = "fill", width = 0.5)+
  theme_bw()+scale_x_continuous(expand = c(0,0))+geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30")+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  ylab("Pathway")+xlab("Relative signalling strength")

p

ggsave(plot = p, filename = "relative.pathways.png", dpi = 500, height = 8, width = 8)

a <- netAnalysis_signalingRole_scatter(DLE.chat)
a <- a$data
b <- netAnalysis_signalingRole_scatter(SLE.chat)
b <- b$data

colnames(a) <- paste(colnames(a), "_DLE")
colnames(b) <- paste(colnames(b), "_SCLE")

c <- cbind(a,b)

c$outgoing.diff <- c$`x _DLE`-c$`x _SCLE`
c$incoming.diff <- c$`y _DLE`-c$`y _SCLE`
c$count.diff <- c$`Count _DLE`-c$`Count _SCLE`

p <- ggplot(c , aes(x = outgoing.diff, y = incoming.diff, fill = `labels _DLE`, label = `labels _DLE`))+
  geom_label_repel()+
  geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(shape = 21, size = 5, colour = "black")+theme_bw()+scale_fill_brewer(palette = "Paired")

p

ggsave(plot = p, filename = "ingoing.outgoing.change.png", dpi = 500,
       height = 5, width = 7)


netVisual_diffInteraction(cle.chat, comparison = c(1,2), vertex.size.max = 10, edge.width.max = 1,
                          label.edge = F)

netAnalysis_signalingRole_scatter(SLE.chat)
netAnalysis_signalingRole_scatter(DLE.chat)
netAnalysis_signalingChanges_scatter(cle.chat, idents.use = "pDC")

 
netVisual_chord_cell(DLE.chat, signaling = "BAFF", title.name = "DLE BAFF signalling")
netVisual_chord_cell(SLE.chat, signaling = "BAFF", title.name = "SCLE BAFF signalling")

netAnalysis_contribution(SLE.chat, signaling = "BAFF")

netVisual_chord_cell(SLE.chat, signaling = "SEMA7", title.name = "BAFF signalling \n in DLE")

cle.chat@idents

SLE.chat@idents
netVisual_bubble(cle.chat,remove.isolate = F, comparison = c(1,2), signaling = "BAFF", 
                 title.name = "Differential BAFF signalling", color.text = c("mediumseagreen", "darkorchid"),
                 dot.size.min = 10, color.heatmap = "viridis")+
  theme(axis.text = element_text(size = 16))

netAnalysis_signalingRole_heatmap(DLE.chat, pattern = "outgoing", slot.name = "netP")

#pull out strongest interactions in pDCs and CD8+ T cells



