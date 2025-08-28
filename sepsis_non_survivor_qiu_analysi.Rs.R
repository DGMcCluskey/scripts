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
library(scDblFinder)
library(celda)#decontx package apparently
library(ggbump)


#LOADING IN DATA===============================================================================================================================================
#first will load in two public datasets at once (Darden et al and Qiu et al) as these have the barcode/feature/matrix file structure
files <- list.files(path = getwd())
files
files <- files[c(49,50,53:62)]
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

sepsis[["percent.mt"]] <- PercentageFeatureSet(sepsis, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(sepsis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
sepsis <- subset(sepsis, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

#add important metadata
patient <- c("Qiu_NS1_T0" = "NS1", "Qiu_NS1_T6" = "NS1",
             "Qiu_NS2_T0" = "NS2", "Qiu_NS2_T6" = "NS2",
             "Qiu_S1_T0" = "S1", "Qiu_S1_T6" = "S1",
             "Qiu_S2_T0" = "S2", "Qiu_S2_T6" = "S2",
             "Qiu_S3_T0" = "S3", "Qiu_S3_T6" = "S3",
             "Qiu_HC1" = "HC1", "Qiu_HC2" = "HC2")

Idents(sepsis) <- "sample"
sepsis <- RenameIdents(sepsis, patient)
sepsis@active.ident
sepsis <- StashIdent(sepsis, save.name = "patient")

sex <- c("NS1" = "M", "NS2" = "F", "S1" = "M", "S2" = "F", "S3" = "F",
         "HC1" = "M", "HC2" = "F")

Idents(sepsis) <- "patient"
sepsis <- RenameIdents(sepsis, sex)
sepsis@active.ident
sepsis <- StashIdent(sepsis, save.name = "sex")
set.seed(61)

ncol(sepsis)

sce <- as.SingleCellExperiment(sepsis)
sce <- scDblFinder(sce)

sepsis$scDblFinder.class <- colData(sce)$scDblFinder.class
sepsis$scDblFinder.score <- colData(sce)$scDblFinder.score

VlnPlot(sepsis, features = "scDblFinder.score")
table(sepsis$scDblFinder.class)

sepsis <- subset(sepsis, scDblFinder.class == "singlet")

sepsis <- NormalizeData(sepsis)
sepsis <- FindVariableFeatures(sepsis)
sepsis <- ScaleData(sepsis, vars.to.regress = c("percent.mt", "nCount_RNA", "sex"))
sepsis <- RunPCA(sepsis)

ElbowPlot(sepsis, ndims = 30)

sepsis <- FindNeighbors(sepsis, dims = 1:16, reduction = "pca")

sepsis <- FindClusters(sepsis, resolution = 0.3, cluster.name = "resolution_0.3")
sepsis <- FindClusters(sepsis, resolution = 0.4, cluster.name = "resolution_0.4")
sepsis <- FindClusters(sepsis, resolution = 0.5, cluster.name = "resolution_0.5")
sepsis <- FindClusters(sepsis, resolution = 0.6, cluster.name = "resolution_0.6")
sepsis <- FindClusters(sepsis, resolution = 0.7, cluster.name = "resolution_0.7")
sepsis <- FindClusters(sepsis, resolution = 0.8, cluster.name = "resolution_0.8")

sepsis <- RunUMAP(sepsis, dims = 1:16, reduction = "pca")

status <- c("NS1" = "NS", "NS2" = "NS", "S1" = "S", "S2" = "S", "S3" = "S",
            "HC1" = "HC", "HC2" = "HC")

Idents(sepsis) <- "patient"
sepsis <- RenameIdents(sepsis, status)
sepsis@active.ident
sepsis <- StashIdent(sepsis, save.name = "status")

#age has to be binned due to poor metadata
age <- c("NS1" = 5, "NS2" = 3, "S1" = 2, "S2" = 3, "S3" = 4,
         "HC1" = 1, "HC2" = 2)

Idents(sepsis) <- "patient"
sepsis <- RenameIdents(sepsis, age)
sepsis@active.ident
sepsis <- StashIdent(sepsis, save.name = "age")

sepsis$age_factor <- factor(sepsis$age)

p <- DimPlot(sepsis, group.by = c("resolution_0.5"), raster = F, label = T)
p

DimPlot(sepsis, split.by = "patient", group.by = "resolution_0.5", raster = T)

p <- DotPlot(sepsis, features = c("CD3D", "CD8A", "CD4", "FOXP3", "TRDV2", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                        "CD1C", "IRF7", "C1QA", "CLEC10A", "CD68", "CSF3R",
                                        "PPBP", "HBB", "FCGR3B", "CEACAM8"), 
             group.by = "resolution_0.5", scale.by = "size",
             dot.min = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)
p

#based on dotplot at resolution 0.5, there are a number of platelet,rbc and neutrophil populations
#i will remove these and re-cluser
Idents(sepsis) <- "resolution_0.5"
sepsis_filtered <- subset(sepsis, idents = c("1","15","11"), invert = T)
sepsis <- sepsis_filtered
rm(sepsis_filtered)

sepsis[["RNA"]] <- split(sepsis[["RNA"]], f = sepsis$patient)
sepsis

sepsis <- NormalizeData(sepsis)
sepsis <- FindVariableFeatures(sepsis)
sepsis <- ScaleData(sepsis, vars.to.regress = c("percent.mt", "nCount_RNA", "sex"))
sepsis <- RunPCA(sepsis)
ElbowPlot(sepsis, ndims = 30, reduction = "harmony")

sepsis <- IntegrateLayers(sepsis, method = HarmonyIntegration,
                                layers = "counts", scale.layer = "scale.data")
sepsis <- JoinLayers(sepsis)

sepsis <- FindNeighbors(sepsis, dims = 1:12, reduction = "harmony")
sepsis <- FindClusters(sepsis, resolution = 0.3, cluster.name = "resolution_0.3")
sepsis <- FindClusters(sepsis, resolution = 0.4, cluster.name = "resolution_0.4")
sepsis <- FindClusters(sepsis, resolution = 0.5, cluster.name = "resolution_0.5")
sepsis <- FindClusters(sepsis, resolution = 0.6, cluster.name = "resolution_0.6")
sepsis <- FindClusters(sepsis, resolution = 0.7, cluster.name = "resolution_0.7")
sepsis <- FindClusters(sepsis, resolution = 0.8, cluster.name = "resolution_0.8")
sepsis
sepsis <- RunUMAP(sepsis, dims = 1:12, reduction = "harmony")
sepsis

p <- DimPlot(sepsis, group.by = c("resolution_0.5"), raster = F, label = T)
p

DimPlot(sepsis, group.by = "resolution_0.5", raster = F)

p <- DotPlot(sepsis, features = c("CD3D", "CD8A", "CD4", "FOXP3", "TRDV2", "NKG7", "GZMB", "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A",
                                  "CD1C", "IRF7", "C1QA", "CLEC10A", "CD68", "CSF3R",
                                  "PPBP", "HBB", "FCGR3B", "CEACAM8"), 
             group.by = "resolution_0.5", scale.by = "size",
             dot.min = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)
p
table(sepsis$resolution_0.5)
#still some dirty clusters, with platelet genes and also markers of multiple lineages. very small so will just remove
Idents(sepsis) <- "resolution_0."
sepsis_filtered <- subset(sepsis, idents = c("7","8","9"), invert = T)
sepsis <- sepsis_filtered
rm(sepsis_all)


#26/06/25 DecontX==============================================================================
#after running rest of script, noticed high heamoglobin, mostly in one patient, suggesting ambient rna
#therefore running DecontX to remove, then will re-cluster and analyse downstream
sepsis <- readRDS(file = "../qiu_reanalysis/qiu_reanalysis_seurat.rds")


# Extract raw counts and metadata
counts <- GetAssayData(sepsis, slot = "counts", assay = "RNA")
cell_metadata <- sepsis@meta.data

# Create SingleCellExperiment
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = cell_metadata
)

# Assign cluster labels
colData(sce)$cluster <- sepsis$labels

sce <- decontX(sce)

clean_counts <- assay(sce, "decontXcounts")

#diagnostics of how many cells are contaminated
table(round(colData(sce)$decontX_contamination, 2))  # distribution
hist(colData(sce)$decontX_contamination)

# Or create a new assay (safer)
sepsis[["decontX"]] <- CreateAssayObject(counts = clean_counts)
DefaultAssay(sepsis) <- "decontX" #CRITICAL
sepsis$decontX_contamination <- colData(sce)$decontX_contamination

class(sepsis[["decontX"]])

FeaturePlot(sepsis_clean, features = c("decontX_contamination"))
sepsis <- sepsis_clean
rm(sepsis_clean)
rm(sce)
rm(sepsis_hbb)
#now re-running standard pipeline
sepsis
sepsis[["decontX"]] <- split(sepsis[["decontX"]], f = sepsis$patient)
sepsis

# Access raw counts of the HBB gene (also allows me to see if decontX has worked)

hbb_counts <- FetchData(sepsis, vars = "HBB", layer = "counts", assay = "RNA")
hbb_counts_decontx <- FetchData(sepsis, vars = "HBB", layer = "counts", assay = "decontX")

# Subset cells with < 6 HBB UMIs
cells_to_keep <- rownames(hbb_counts_decontx)[hbb_counts_decontx$HBB < 6]

# Subset the Seurat object
sepsis <- subset(sepsis, cells = cells_to_keep)
sepsis <- NormalizeData(sepsis, assay = "decontX")
sepsis <- FindVariableFeatures(sepsis, assay = "decontX")
sepsis <- ScaleData(sepsis, vars.to.regress = c("percent.mt", "nCount_RNA", "sex"), assay = "decontX")
sepsis <- RunPCA(sepsis, assay = "decontX")
ElbowPlot(sepsis, ndims = 30, reduction = "pca")

sepsis <- IntegrateLayers(sepsis, method = HarmonyIntegration,
                          layers = "counts", scale.layer = "scale.data")
sepsis <- JoinLayers(sepsis)

sepsis <- FindNeighbors(sepsis, dims = 1:12, reduction = "harmony")
sepsis <- FindClusters(sepsis, resolution = 0.3, cluster.name = "resolution_0.3")
sepsis <- FindClusters(sepsis, resolution = 0.4, cluster.name = "resolution_0.4")
sepsis <- FindClusters(sepsis, resolution = 0.5, cluster.name = "resolution_0.5")
sepsis <- FindClusters(sepsis, resolution = 0.6, cluster.name = "resolution_0.6")
sepsis <- FindClusters(sepsis, resolution = 0.7, cluster.name = "resolution_0.7")
sepsis <- FindClusters(sepsis, resolution = 0.8, cluster.name = "resolution_0.8")
sepsis
sepsis <- RunUMAP(sepsis, dims = 1:12, reduction = "harmony")
sepsis

p <- DimPlot(sepsis, group.by = c("resolution_0.5"), raster = F, label = T)
p

DimPlot(sepsis, split.by = "patient", group.by = "resolution_0.5")


#ANNOTATION===========================================================
DimPlot(sepsis, group.by = "resolution_0.7", raster = F, label = T)

FeaturePlot(sepsis, features = c("CD14", "FCGR3A"))

p <- DotPlot(sepsis, features = c("CD3D", "CD8A", "CD4", "FOXP3", "SELL", "LEF1", "CCR7", "IL7R", 
                                  "NKG7", "GZMB", 
                                  "CD79A","IGHD", "CD27", "BLIMP1", "PRDM1", "LYZ", "CD14", "FCGR3A",
                                  "CD1C", "IRF7", "C1QA", "CLEC10A", "CD68", "CSF3R",
                                  "PPBP", "HBB"), 
             group.by = "resolution_0.7", scale.by = "size", assay = "RNA",
             dot.min = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)
p

annotation <- c("0" = "T_naive", "1" = "B", "2" = "Classical_mono",
                "3" = "NK", "4" = "T_cm", "5" = "T_cm", "6" = "B", "7" = "T_cytotoxic",
                "8" = "Intermediate_mono", "9" = "B")

Idents(sepsis) <- "resolution_0.7"
sepsis <- RenameIdents(sepsis, annotation)
sepsis@active.ident
sepsis <- StashIdent(sepsis, save.name = "labels")

ncol(sepsis)

sepsis$labels <- factor(sepsis$labels, levels = c("T_naive", "T_cm", "T_cytotoxic", "NK",
                                                  "B", "Classical_mono", "Intermediate_mono"))

p <- DimPlot(sepsis, group.by = "labels", raster = F, label = T, label.size = 5)+
  scale_colour_brewer(palette = "Paired")+
  theme_void()+NoLegend()+ggtitle("total cells = 29,226")
p

p2 <- DimPlot(sepsis, split.by = "patient", group.by = "labels", raster = F, label = F)+
  scale_colour_brewer(palette = "Paired")+theme_void()+ggtitle(NULL)

p2

p3 <- DotPlot(sepsis, features = c("CD3D", "CD8A", "CD4", "SELL", "LEF1", "CCR7", 
                                  "NCAM1", "GZMB",  
                                  "CD79A","IGHD", "CD27", "LYZ", "CD14", "FCGR3A"), 
             group.by = "labels", scale.by = "size",
             dot.min = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)+xlab(NULL)+ylab(NULL)
p3

combo <- (p | p3) / p2
combo

ggsave(plot = combo, filename = "../qiu_reanalysis/umaps_dotlot.png", dpi = 500,
       height = 8, width = 14)

saveRDS(sepsis, file = "../qiu_reanalysis/qiu_reanalysis_seurat.rds", compress = F)
sepsis <- readRDS(sepsis, file = "../qiu_reanalysis/qiu_reanalysis_seurat.rds", compress = F)

sample_colors <- c("indianred2", "steelblue")
p <- Stacked_VlnPlot(monocytes, features = c("LYZ", "CD14", "FCGR3A", "PPBP", "PF4", "SELP",
                                             "TBXA2R", "F2RL3"), pt.size = 0.5, group.by = "status", split.by = "labels",
                colors_use = sample_colors)

p

ggsave(plot = p, filename = "../qiu_reanalysis/violins.monocyte.platelet.markers.png", dpi = 500,
       height = 8, width = 6)

FeatureScatter(monocytes, feature1 = "LYZ", feature2 = "THBS1", group.by = "labels")
#PROPORTIONS======================================================================
sepsis$status <- factor(sepsis$status, levels = c("HC", "S", "NS"))

p <- Proportion_Plot(sepsis, group_by_var = "labels", split.by = "patient")+
  scale_fill_brewer(palette = "Paired")
p
#ratio of classical to intermediate monocytes
Idents(sepsis) <- "labels"
monocytes$status <- factor(monocytes$status, levels = c("HC", "S", "NS"))
monocytes <- subset(sepsis, idents = c("Classical_mono", "Intermediate_mono"))

#get same colours in palette for the monocytes
cols <- brewer.pal(name = "Paired", n = 7)
cols <- cols[6:7]

p2 <- Proportion_Plot(monocytes, group_by_var = "labels", split.by = "patient")+
  scale_fill_manual(values = cols)+ggtitle("Proportion of monocyte subsets")
p2

combo <- p+p2+plot_layout(ncol = 2)
combo

ggsave(plot = combo, filename = "../qiu_reanalysis/proportion.plots.png", dpi = 500,
       height = 6, width = 14)


#correlate monocyte proportions with serum cytokines provided in their paper
prop_data <- p2$data
prop_data$IL6 <- c(30.2, 30.2, 142.3, 142.3, 133, 133, 2.48, 2.48, 0.31, 0.31, 0, 0, 0.002, 0.002)
prop_data$IL8 <- c(6.65, 6.65, 27.2, 27.2, 41.7, 41.7, 0.61, 0.61, 0.4, 0.4, 0.03, 0.03, 0.026, 0.026)
prop_data$IL10 <- c(0.13, 0.13, 9.71, 9.71, 0.52, 0.52, 0.15, 0.15, 0.39, 0.39, 0, 0, 0, 0)
prop_data$SOFA <- c(11,11,16,16,11,11,15,15,7,7,NA,NA,NA,NA)
prop_data$APACHEII <- c(18,18,38,38,31,31,41,41,19,19,NA,NA,NA,NA)

ggplot(prop_data, aes(x = APACHEII, y = value, colour = split.by))+geom_point()+
  theme_bw()+facet_wrap(~Cluster)

#DIFFERENTIAL EXPRESSION ANALYSIS======================================================================
#attempting to run DEA on all cell subsets
celltypes <- unique(sepsis$labels)

de_results <- map(celltypes, function(ct) {
  message("Processing: ", ct)
  
  # Subset using base R condition directly
  subset_obj <- subset(sepsis, subset = labels == ct)
  
  # Check condition balance
  subtype_counts <- table(subset_obj$status)
  if (length(subtype_counts) < 2 || any(subtype_counts < 10)) {
    message(paste("Skipping", ct, "due to insufficient cells per condition"))
    return(NULL)
  }
  
  # Run DE with MAST
  markers <- FindMarkers(subset_obj,
                         ident.1 = "S",
                         ident.2 = "HC",
                         group.by = "status",
                         test.use = "MAST",
                         latent.vars = c("sex", "age_factor"),
                         logfc.threshold = 0.1,
                         min.pct = 0.1)
  
  markers$celltype <- ct
  return(markers)
})

# Combine all results
names(de_results) <- celltypes[!map_lgl(de_results, is.null)]

de_results_df <- bind_rows(de_results, .id = "celltype")
de_results_df <- de_results_df[!grepl("^HB[AB]", rownames(de_results_df)), ]

de_results_filtered <- map(de_results, function(df) {
  if (is.null(df)) return(NULL)
  dplyr::filter(df, p_val_adj < 0.05)
})
de_results_filtered_df <- bind_rows(de_results_filtered, .id = "celltype")
de_results_filtered_df <- de_results_filtered_df[!grepl("^HB[AB]", de_results_filtered_df$gene), ]


write.csv(de_results_df, file = "../qiu_reanalysis/qiu_reanalysis_s_vs_hc_unfiltered.csv")
write.csv(de_results_filtered_df, file = "../qiu_reanalysis/qiu_reanalysis_s_vs_hc_significant.csv")


#also runnig on total monocytes, to do this i'll subset them out first to avoid re-doing degs for other celltypes
Idents(sepsis) <- "labels"
monocytes <- subset(sepsis, idents = c("Classical_mono", "Intermediate_mono"))
Idents(monocytes) <- "status"

DimPlot(monocytes, group.by = "status")

monocyte_degs <- FindMarkers(monocytes, ident.1 = "S", ident.2 = "HC",
                             test.use = "MAST",
                             latent.vars = c("sex", "age_factor"),
                             logfc.threshold = 0.1,
                             min.pct = 0.1
                             )

monocyte_degs$gene <- rownames(monocyte_degs)

monocyte_degs <- monocyte_degs[!grepl("^HB[AB]", monocyte_degs$gene), ]

monocyte_degs_significant <- filter(monocyte_degs, p_val_adj < 0.05)


write.csv(monocyte_degs, file = "../qiu_reanalysis/total_mono_s_vs_hc_unfiltered.csv")
write.csv(monocyte_degs_significant, file = "../qiu_reanalysis/total_mono_s_vs_hc_filtered.csv")


#load in gene lists as needed
de_results_filtered_df <- read.csv("../qiu_reanalysis/qiu_reanalysis_ns_vs_s_significant.csv")
de_results_filtered_df$gene <- gsub("\\.\\.\\.[0-9]+$", "", de_results_filtered_df$X)
de_results_unfiltered_df <- read.csv("../qiu_reanalysis/qiu_reanalysis_ns_vs_s_unfiltered.csv")
de_results_unfiltered_df$gene <- gsub("\\.\\.\\.[0-9]+$", "", de_results_unfiltered_df$X)

deg_numbers <- de_results_filtered_df %>%
  group_by(celltype) %>%
  filter(avg_log2FC > 0.5 | avg_log2FC < -0.5) %>% 
  summarise(n_DEGs = n())

deg_numbers$celltype <- factor(deg_numbers$celltype, levels = c("T_naive", "T_cm", "T_cytotoxic", "NK",
                                                  "B", "Classical_mono", "Intermediate_mono"))



p <- ggplot(deg_numbers, aes(y = celltype, x = n_DEGs, fill = celltype))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  ggtitle("Number of DEGs in non-survival vs survival sepsis")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)+
  geom_text(aes(label = n_DEGs), 
            position = position_dodge(width = 0.9), 
            hjust = -0.2, size = 4, color = "black")+xlim(0,1200)+
  theme(legend.position = "none")
p

ggsave(plot = p, filename = "../qiu_reanalysis/number_degs.png", dpi = 500,
       height = 5, width = 6)


#volcanos
volcano_input <- filter(de_results_unfiltered_df, celltype %in% c("Classical_mono", "Intermediate_mono"))
volcano_input$sig <- ifelse(volcano_input$avg_log2FC > 0.5 & volcano_input$p_val_adj < 0.05, "Up",
                       ifelse(volcano_input$avg_log2FC < -0.5 & volcano_input$p_val_adj < 0.05, "Down",
                              "NS"))
volcano_input$sig <- factor(volcano_input$sig, levels = c("Up", "Down", "NS"))

top_genes <- volcano_input %>%
  group_by(celltype) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

bottom_genes <- volcano_input %>%
  group_by(celltype) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 10)

label_genes <- bind_rows(top_genes, bottom_genes)



volcano <- ggplot(volcano_input, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+
  geom_point(size = 2)+facet_wrap(~celltype, ncol = 1, scales = "free_y")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")+
  scale_colour_manual(values = c("darkorchid3", "mediumseagreen", "grey50"))+theme_bw()+
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = Inf)+
  ylim(0, 50)
volcano

ggsave(plot = volcano, filename = "../qiu_reanalysis/monocyte_volcanos.png", dpi = 500,
       height = 8, width = 8)

#degs up and down and with a logfc cutoff

deg_numbers <- de_results_filtered_df %>%
  mutate(Direction = case_when(
    avg_log2FC > 0.5  ~ "Increased in non-survivors",
    avg_log2FC < -0.5 ~ "Increased in survivors",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Direction)) %>%
  group_by(celltype, Direction) %>%
  summarise(n_DEGs = n(), .groups = "drop")
deg_numbers

p <- ggplot(deg_numbers, aes(y = celltype, x = n_DEGs, fill = Direction))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen"))+
  ggtitle("Number of DEGs in non-survival vs survival sepsis \n
          (Fold change > |2|)")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)
p

de_results_filtered_df$gene <- rownames(de_results_filtered_df)

#get top 25 DEGs from each celltype to inspect
top25_degs <- de_results_filtered_df %>%
  group_by(celltype) %>%
  slice_max(order_by = abs(avg_log2FC), n = 25, with_ties = FALSE) %>%
  ungroup()

top25_degs$gene <- gsub("\\.\\.\\.[0-9]+$", "", top25_degs$gene)


#GSEA============================================================================
#soi = subset of interest
#this will change as i analyse different populations - be careful to save plots with correct names
#to be extra safe i load in genes from csv again
de_results_df <-read.csv("../qiu_reanalysis/qiu_reanalysis_ns_vs_s_unfiltered.csv")
de_results_df$gene <- de_results_df$X
de_results_df$gene <- gsub("\\.\\.\\.[0-9]+$", "", de_results_df$gene)
soi_degs <- filter(de_results_df, celltype %in% "NK")

soi_genes <- soi_degs$avg_log2FC
names(soi_genes) <- soi_degs$gene


soi_entrez <- bitr(names(soi_genes), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

soi_genes <- soi_genes[soi_entrez$SYMBOL]
names(soi_genes) <- soi_entrez$ENTREZID

#save genes here if needed to then use in comparecluster
classical_genes <- soi_genes
intermediate_genes <- soi_genes
nk_genes <- soi_genes


gene_list <- list(classical = classical_genes, intermediate= intermediate_genes)
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

top_terms <- cc_results@geneClusters$classical[1:5]
for (term in top_terms) {
  print(gseaplot(cc_results, geneSetID = term, title = term))
}

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

#alternatively rank by p value
# Top 20 enriched by FDR per cluster
top_enriched_p <- comparison_results %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)


# Combine enriched and unenriched
top_pathways <- bind_rows(top_enriched, top_unenriched) %>%
  ungroup()

top_pathways$direction <- ifelse(top_pathways$NES > 0, "Enriched in non-survivors",
                                 ifelse(top_pathways$NES < 0, "Enriched in survivors", "Other"))

#visualising by bar instead
p <- ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES), fill = p.adjust))+
  geom_col()+theme_bw()+
  facet_wrap(~Cluster, scales = "free_y", ncol = 1)+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen", "grey50"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme(strip.text = element_text(size = 18, colour = "black"))+
  ylab(NULL)+xlab("Normalised enrichment score")

p

ggsave(plot = p, filename = "../qiu_reanalysis/GSEgo_monocytes_ns_vs_s.png", dpi = 500,
       height = 10, width = 16)

#p value version
p <- ggplot(top_enriched_p, aes(x = NES, y = reorder(Description, -log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+theme_bw()+
  facet_wrap(~Cluster, scales = "free_y", ncol = 1)+
  scale_fill_distiller(palette = "Reds")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  theme(strip.text = element_text(size = 14, colour = "black"))+
  ylab(NULL)+xlab("Normalised enrichment score")

p



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



#analysis when only running one cell subset
#if running on just one group without comparison
genes <- sort(soi_genes, decreasing = T)
up_NS_down_S_results <- gseGO(genes, OrgDb = org.Hs.eg.db)

up_NS_down_S_results <- clusterProfiler::simplify(up_NS_down_S_results, cutoff=0.7, by="p.adjust",
                                        select_fun=min)

dotplot(up_NS_down_S_results)

up_NS_down_S_results <- up_NS_down_S_resultss@result


dotplot(up_NS_down_S_results)

top_enriched <- up_NS_down_S_results %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_unenriched <- up_NS_down_S_results %>%
  arrange(NES) %>%
  slice_head(n = 10)

# Combine enriched and unenriched
top_pathways <- bind_rows(top_enriched, top_unenriched) %>%
  ungroup()

top_pathways$direction <- ifelse(top_pathways$NES > 0, "Enriched in non-survivors",
                                 ifelse(top_pathways$NES < 0, "Enriched in survivors", "Other"))

#visualising by bar instead
p <- ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES), fill = direction))+
  geom_col()+theme_bw()+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen", "grey50"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme(strip.text = element_text(size = 18, colour = "black"))+
  ylab(NULL)+xlab("Normalised enrichment score")

p



#load in top TFs after running in decoupler
classical_tf <- read.csv("../qiu_reanalysis/classical_mono_tfs_ns_vs_s.csv")
intermediate_tf <- read.csv("../qiu_reanalysis/intermediate_mono_tfs_ns_vs_s.csv")

up_ns_down_s_NSdirection <- read.csv("../qiu_reanalysis/monocyte_TFs_up_ns_down_s.csv")
up_ns_down_s_Sdirection <- read.csv("../qiu_reanalysis/monocyte_TFs_up_ns_down_s_Sdirection.csv")

# Top x by logFC
top <- up_ns_down_s_NSdirection %>%
  arrange(desc(logFC_sample)) %>%
  slice_head(n = 15)

# Bottom x by logFC
bottom <- classical_tf %>%
  arrange(logFC_sample) %>%
  slice_head(n = 15)

# Combine
top_bottom <- bind_rows(top, bottom)
top_bottom <- distinct(top_bottom)


top_bottom$direction <- ifelse(top_bottom$logFC_sample > 0, "Enriched in non-survivors",
                               "Enriched in survivors")

p <- ggplot(top_bottom, aes(y = reorder(TF, logFC_sample), x = logFC_sample, fill = direction))+geom_col()+
  theme_bw()+scale_fill_manual(values = c("darkorchid", "mediumseagreen"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  xlab("Enrichment score")+ylab(NULL)
p

ggsave(plot = p, filename = "../qiu_reanalysis/intermediate_mono_top_TFs_ns_vs_s.png", dpi = 500,
       height = 6, width = 6)

#pulling in downstream genes of TFs of interest - these were gathered from the decoupler results in python and exported to csv
USF1 <- read.csv("../qiu_reanalysis/USF1_target_genes.csv")
USF1 <- USF1[2]
USF2 <- read.csv("../qiu_reanalysis/USF2_target_genes.csv")
USF2 <- USF2[2]

de_results_df <-read.csv("../qiu_reanalysis/qiu_reanalysis_ns_vs_s_significant.csv")
de_results_df$gene <- de_results_df$X
de_results_df$gene <- gsub("\\.\\.\\.[0-9]+$", "", de_results_df$gene)
soi_degs <- filter(de_results_df, celltype %in% "Classical_mono")
soi_degs <- filter(soi_degs, gene %in% USF1$X0)

volcano_input <- soi_degs

volcano_input$sig <- ifelse(volcano_input$avg_log2FC > 0.5 & volcano_input$p_val_adj < 0.05, "Up",
                            ifelse(volcano_input$avg_log2FC < -0.5 & volcano_input$p_val_adj < 0.05, "Down",
                                   "NS"))
volcano_input$sig <- factor(volcano_input$sig, levels = c("Up", "Down", "NS"))

top_genes <- volcano_input %>%
  group_by(celltype) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

bottom_genes <- volcano_input %>%
  group_by(celltype) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 10)

label_genes <- bind_rows(top_genes, bottom_genes)
label_genes <- distinct(label_genes)

volcano <- ggplot(volcano_input, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig))+
  geom_point(size = 2)+facet_wrap(~celltype, ncol = 1, scales = "free_y")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+geom_vline(xintercept = 0.5, linetype = "dashed")+
  geom_vline(xintercept = -0.5, linetype = "dashed")+
  scale_colour_manual(values = c("darkorchid3", "mediumseagreen", "grey50"))+theme_bw()+
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = Inf) 
volcano


#pulling in all tfs with padj for inspection
TFs <- read.csv("../qiu_reanalysis/TFs_and_padj_mono_up_S_down_or_flat_NS.csv")
#TFs <- filter(TFs, padj < 0.05)

#TFs$direction <- ifelse(TFs$logFC_sample > 0, "Enriched non-survivor","Enriched survivor")

# Top
top <- TFs %>%
  arrange(desc(logfc_sample)) %>%
  slice_head(n = 10)

# Bottom 
bottom <- TFs %>%
  arrange(logfc_sample) %>%
  slice_head(n = 10)

# Combine
top_bottom <- bind_rows(top, bottom)
top_bottom <- distinct(top_bottom)


p <- ggplot(top_bottom, aes(x = logfc_sample, y = reorder(X, logfc_sample), fill = padj))+geom_col(colour = "black")+
  theme_bw()+ylab(NULL)+xlab("Enrichment")+theme(axis.text = element_text(size = 16, colour = "black"))+
  scale_fill_distiller(palette = "Blues")
p
ggsave(plot = p, filename = "../qiu_reanalysis/up_or_higher_NS_down_or_flat_S_TFs.png", dpi = 500,
       height = 4, width = 8)


#plotting downstream genes of interest
classical_mono <- subset(sepsis, idents = "Classical_mono")
Idents(classical_mono) <- "status"
classical_mono <- subset(classical_mono, idents = c("NS", "S"))
p <- DotPlot(classical_mono, features = c("IL10", "CCL7", "CXCL3", "CXCL8", "THBS1", "ATF4", "HIF1A", "HIF1AN", "HMOX1",
                                     "EIF4E", "ETS2", "IRAK3", "MAPK1", "PIK3CA"),
        group.by = "patient", scale.by = "size", dot.scale = 10, dot.min = 0.01)+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis(option = "F", direction = -1)+ylab("Sample")+xlab(NULL)+coord_flip()

p

ggsave(plot = p, filename = "../qiu_reanalysis/classical_mono_interesting_genes_dotplot.png", dpi = 500,
       height = 8, width = 8, bg = "white")


#VENN DEG OVERLAPS=================================================================================
ns_s <- read.csv("../qiu_reanalysis/total_mono_ns_vs_s_filtered.csv")
ns_s <- filter(ns_s, avg_log2FC > 0.5)
ns_hc <- read.csv("../qiu_reanalysis/total_mono_ns_vs_hc_filtered.csv")
ns_hc <- filter(ns_hc, avg_log2FC > 0.5)
s_hc <- read.csv("../qiu_reanalysis/total_mono_s_vs_hc_filtered.csv")
s_hc <- filter(s_hc, avg_log2FC > 0.5)

genomesize = nrow(monocytes) 

go.obj <- newGeneOverlap(ns_hc$X, s_hc$X, genomesize)
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

overlap <- as.data.frame(go.obj@intersection)
ns_only <- as.data.frame(setdiff(go.obj@listA, go.obj@intersection))
names(ns_only)[1] <- "gene"
s_only <- as.data.frame(setdiff(go.obj@listB, go.obj@intersection))
names(s_only)[1] <- "gene"

input <- list("Up non-survivors \n vs HC" = ns_hc$X, "Up survivors \n vs HC" = s_hc$X, 
              "Up non-survivors vs survivors" = ns_s$X)


venn.diagram(input, filename = "../qiu_reanalysis/mono_degs_venn.upregulated.png", fill = c("mediumseagreen", "darkorchid", 
                                                                                "goldenrod2"),
             resolution = 500, col = "white", cat.pos = c(340,20,180), cat.cex = 2,
             cex = 2, alpha = 0.5, label.col = "black", fontface = "bold")


intersect_ns_nss <- intersect(input[[1]], input[[3]])
intersect_ns_nss

intersect_s_nss <- intersect(input[[2]], input[[3]])
intersect_s_nss

ns_s_only <- setdiff(input[["Up non-survivors vs survivors"]], 
                     union(input[["Up non-survivors \n vs HC"]], 
                           input[["Up survivors \n vs HC"]]))
ns_s_only

genes_in_all_three <- Reduce(intersect, input)
genes_in_all_three

write.csv(genes_in_all_three, file = "../qiu_reanalysis/genes_up_all_comparisons.csv")
write.csv(ns_s_only, file = "../qiu_reanalysis/genes_up_NS_vs_S_only.csv")

ns_hc_only <- setdiff(input[["Up non-survivors \n vs HC"]], 
                     union(input[["Up non-survivors vs survivors"]], 
                           input[["Up survivors \n vs HC"]]))

write.csv(ns_hc_only, file = "../qiu_reanalysis/genes_up_NS_vs_HC_only.csv")

DotPlot(monocytes, features = genes_in_all_three, group.by = "status", 
        scale.by = "size", dot.scale = 10, dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)+ylab("Sample")+xlab(NULL)+coord_flip()

ns_s_only <- ns_s_only[!grepl("^RPL", ns_s_only)]
ns_s_only <- ns_s_only[!grepl("^RPS", ns_s_only)]


shared_ns_genes <- intersect(input[["Up non-survivors \n vs HC"]], 
                             input[["Up non-survivors vs survivors"]])

# Remove any that are also up in s_hc
ns_s_ns_hc_genes <- setdiff(shared_ns_genes, input[["Up survivors \n vs HC"]])



p <- DotPlot(monocytes, features = ns_s_ns_hc_genes, group.by = "status", 
        scale.by = "size", dot.scale = 10, dot.min = 0.01)+
  scale_colour_viridis(option = "F", direction = -1)+ylab("Sample")+xlab(NULL)+coord_flip()
p

ggsave(plot = p, filename = "../qiu_reanalysis/dotplot_genes_ns_s_ns_hc.png", dpi = 500,
       height = 5, width = 6, bg = "white")

VlnPlot(monocytes, features = c("CXCL8", "CXCL3", "HIF1A"), stack = T, group.by = "status", pt.size = 0)

#looking at directionality of all differentially expressed genes to try and find ones going in common direction/patterns
ns_hc <- read.csv("../qiu_reanalysis/total_mono_ns_vs_hc_unfiltered.csv")
s_hc <- read.csv("../qiu_reanalysis/total_mono_s_vs_hc_unfiltered.csv")


#ns_hc <- filter(ns_hc, p_val_adj < 0.05)
#s_hc <- filter(s_hc, p_val_adj < 0.05)

#merge them and add a suffix
mono_degs <- merge(ns_hc, s_hc, by = "gene", suffixes = c("-NS", "-S"), all = F)
#remove uneeded columns
mono_degs <- mono_degs[, c(1,4,7,10,13)]


#keep only genes where they are at least significant in one of the conditions
mono_degs$keep <- ifelse(mono_degs$`p_val_adj-NS` < 0.05 | mono_degs$`p_val_adj-S` < 0.05, "yes", "no")

mono_degs <- filter(mono_degs, keep == "yes")

#add a healthy fold change just for plotting
mono_degs$`avg_log2FC-HC` <- 0

mono_degs <- pivot_longer(mono_degs, names_to = c(".value", "condition"), names_sep = "-", cols = -gene)

#remove the NA condition introduced in the pivot_longer

mono_degs <- filter(mono_degs, condition %in% c("HC", "S", "NS"))


mono_degs$condition <- factor(mono_degs$condition, levels = c("HC", "S", "NS"))

ggplot(mono_degs, aes(x = condition, y = avg_log2FC, group = gene))+geom_line()+
  theme_bw()

#initial plot is too crowded/no meaningful grouping
#i will add groups based on whether a gene goes up in both groups, down in both groups, up in one or the other
df_labeled <- mono_degs %>%
  group_by(gene) %>%
  summarise(
    log2fc_NS = avg_log2FC[condition == "NS"],
    log2fc_S  = avg_log2FC[condition == "S"],
    padj_NS = p_val_adj[condition == "NS"],
    padj_S = p_val_adj[condition == "S"],
    log2fc_status = case_when(
      log2fc_NS > 0.5 & log2fc_S > 0.5 & padj_NS < 0.05 & padj_S < 0.05 ~ "up_in_both",
      log2fc_NS < -0.5 & log2fc_S < -0.5 & padj_NS < 0.05 & padj_S < 0.05 ~ "down_in_both",
      log2fc_NS > 0.5 & log2fc_S < -0.5 & padj_NS < 0.05 & padj_S < 0.05 ~ "up_in_NS_down_in_S",
      log2fc_NS < -0.5 & log2fc_S > 0.5 & padj_NS < 0.05 & padj_S < 0.05 ~ "down_in_NS_up_in_S",
      log2fc_NS > 0.5 & padj_NS < 0.05 & padj_S >= 0.05 ~ "up_in_NS_only",
      log2fc_S > 0.5 & padj_NS >= 0.05 & padj_S < 0.05 ~ "up_in_S_only",
      log2fc_NS < -0.5 & padj_NS < 0.05 & padj_S >= 0.05 ~ "down_in_NS_only",
      log2fc_S < -0.5 & padj_NS >= 0.05 & padj_S < 0.05 ~ "down_in_S_only",
      TRUE ~ "other"
    ),
    .groups = "drop"
  )

df_labeled <- filter(df_labeled, log2fc_status != "other")

df_labeled$log2fc_HC <- 0
df_labeled$padj_HC <- 1

#df_labeled <- df_labeled[, -c(4,5)]

mono_degs <- df_labeled %>%
  pivot_longer(
    cols = c(log2fc_NS, log2fc_S, log2fc_HC,  padj_NS, padj_S, padj_HC),
    names_to = c(".value", "condition"),
    names_pattern = "(log2fc|padj)_(NS|S|HC)"
  )


mono_degs$condition <- factor(mono_degs$condition, levels = c("S", "HC", "NS"),
                              labels = c("S", "HC", "NS"))

mono_degs$log2fc_status <- factor(mono_degs$log2fc_status,
                                  levels = c("up_in_both", "down_in_both",
                                             "up_in_NS_only", "up_in_S_only",
                                             "down_in_NS_only", "down_in_S_only",
                                             "up_in_NS_down_in_S", "down_in_NS_up_in_S"), 
                                  labels = c("Up in NS & S", "Down in NS & S",
                                             "Up in NS only","Up in S only",
                                             "Down in NS only","Down in S only",
                                             "Up in NS & down in S","Down in NS & up in S"))


mono_degs$sig <- ifelse(mono_degs$padj < 0.05, "p < 0.05", "p >= 0.05")

p <- ggplot(mono_degs, aes(x = condition, y = log2fc, fill = sig, 
                      group = gene))+
  geom_bump(colour = "grey30")+geom_point(shape = 21)+facet_wrap(~log2fc_status, scales = "free_y", ncol = 2)+theme_bw()+
  scale_fill_manual(values = c("goldenrod2", "grey30"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+ylab("Log2FC compared to healthy controls")+
  xlab(NULL)+theme(strip.text = element_text(size = 16, colour = "black"))
p

ggsave(plot = p, filename = "../qiu_reanalysis/divergent.DEGs.png", dpi = 500,
       height = 10, width = 10)

#get counts of genes
counts <- mono_degs %>% group_by(log2fc_status) %>% summarise(n = n_distinct(gene))


#pull out gene groups to look at GSEA  and TF enrichment
write.csv(mono_degs, file = "../qiu_reanalysis/monocyte_deg_directionality_log2fc_0.5.csv")
goi <- filter(mono_degs, log2fc_status %in% "Up in NS & down in S")
goi <- filter(goi, condition %in% "NS")

#filtering the up in ns & s list for genes significantly up in NS relative to S
#will then add these genes to the GSEA below and the TF analysis
ns_vs_s <- read.csv("../qiu_reanalysis/total_mono_ns_vs_s_unfiltered.csv")
ns_vs_s <- filter(ns_vs_s, p_val_adj < 0.05)
ns_vs_s <- filter(ns_vs_s, avg_log2FC > 0)

#before further analysis, will do a plot of these genes with the geom_bump visualisation
divergent <- read.csv("../qiu_reanalysis/monocyte_deg_directionality.csv")

input <- filter(divergent, gene %in% ns_vs_s$X)
input$sig <- ifelse(input$padj < 0.05, "p < 0.05", "p >= 0.05")
input <- filter(input, log2fc_status %in% "Up in NS & S")
input$condition <- factor(input$condition, levels = c("S", "HC", "NS"))
input$log2fc_status <- factor(input$log2fc_status, levels = "Up in NS & S", labels = "Significantly higher in NS than S")
p <- ggplot(input, aes(x = condition, y = log2fc, fill = sig, 
                           group = gene))+
  geom_bump(colour = "grey30")+geom_point(shape = 21)+facet_wrap(~log2fc_status, scales = "free_y", ncol = 2)+theme_bw()+
  scale_fill_manual(values = c("goldenrod2", "grey30"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+ylab("Log2FC compared to healthy controls")+
  xlab(NULL)+theme(strip.text = element_text(size = 16, colour = "black"))
p

ggsave(plot = p, filename = "../qiu_reanalysis/divergent.higher.ns.than.s.png", dpi = 500,
       height = 4, width = 5)




divergent <- filter(divergent, log2fc_status %in% c("Up in NS & S"))
divergent <- filter(divergent, condition %in% "NS")



#UpsetPlot on genes as a different/summary visualisation
library(UpSetR)


genes_upset <- mono_degs %>%
  mutate(
    up_in_NS_only = str_detect(log2fc_status, "Up in NS only"),
    down_in_NS_only = str_detect(log2fc_status, "Down in NS only"),
    up_in_S_only = str_detect(log2fc_status, "Up in S only"),
    down_in_S_only = str_detect(log2fc_status, "Down in S only"),
    up_in_ns_and_s = str_detect(log2fc_status, "Up in NS & S"),
    down_in_ns_and_s = str_detect(log2fc_status, "Down in NS & S"),
    up_in_ns_and_down_in_s = str_detect(log2fc_status, "Up in NS & down in S"),
    down_in_ns_and_up_in_s = str_detect(log2fc_status, "Down in NS & up in S")
  )

genes_upset <- filter(genes_upset, condition %in% "NS")


columns <- c(
  "up_in_NS_only", "down_in_NS_only", "up_in_S_only", "down_in_S_only",
  "up_in_ns_and_s", "down_in_ns_and_s", "up_in_ns_and_down_in_s", "down_in_ns_and_up_in_s"
)

genes_upset[columns] <- lapply(genes_upset[columns], as.integer)
genes_upset <- as.data.frame(genes_upset)
upset(genes_upset, sets = columns)




#GSEA on divergent genes
#here im bringing in genes from further down the script (where i split genes based on directionality)
divergent <- read.csv("../qiu_reanalysis/monocyte_deg_directionality.csv")

divergent <- filter(divergent, log2fc_status %in% c("Up in NS only", "Up in NS & down in S"))
divergent <- filter(divergent, condition %in% "NS")

#plus the genes that are higher in NS than S (from just above this section, called input)
input <- filter(input, condition %in% "NS")

divergent_input <- rbind(divergent, input)
#save this for TF input
write.csv(divergent_input, file = "../qiu_reanalysis/NS_up_or_exclusive_genes.csv")

soi_genes <- divergent_input$log2fc
names(soi_genes) <- divergent_input$gene
soi_genes
soi_entrez <- bitr(names(soi_genes), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

soi_genes <- soi_genes[soi_entrez$SYMBOL]
names(soi_genes) <- soi_entrez$ENTREZID

input <- names(soi_genes)

go_results <- enrichGO(input, OrgDb = org.Hs.eg.db)

pathways <- go_results@result
pathways <- filter(pathways, padj < 0.05)









#MACHINE LEARNING APPROACHES=======================================================================

#the goal of this is to go from correlations to CAUSAL relationships

#mediation analysis
library(mediation)

#for this i need to assess a mediator (e.g. gene expression, or activity) at the patient level

input <- as.data.frame(AverageExpression(monocytes, features = "JUN", group.by = "patient", assays = "RNA"))
input <- pivot_longer(input, names_to = "patient_id", values_to = "mediator", 1:7)
input$group <- c("NS", "NS", "S", "S", "S", "HC", "HC")
input$outcome <- c("Died", "Died", "Survivor", "Survivor", "Survivor", "Survivor", "Survivor")

input$outcome <- ifelse(input$outcome == "Died", 1, 0)

input <- filter(input, group %in% c("NS", "S"))

# Mediator model: TF activity depends on group
med.fit <- lm(mediator ~ group, data = input)
med.fit
# Outcome model: survival depends on TF activity and group
out.fit <- glm(outcome ~ group + mediator, data = input, family = binomial)
out.fit
str(input$group)
input$group <- factor(input$group)
unique(med.fit$model$group)
unique(out.fit$model$group)
str(input$mediator)
table(input$mediator)
str(input$outcome)
table(input$outcome)
sum(is.na(input$group))
sum(is.na(input$mediator))
sum(is.na(input$outcome))
contrasts(input$group)

# Mediation
med.out <- mediate(med.fit, out.fit, treat = "group", mediator = "mediator", boot = TRUE)
summary(med.out)