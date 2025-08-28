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
library(readr)

download.packages(c("Seurat", "SeuratData", "SeuratObject", 
                    "ggplot2" "tidyr", "dplyr"), destdir = "r_pkgs", type = "source")
#attempting to run DEA on all cell subsets
setwd("/myriadfs/home/rmhadgm/Scratch")

sepsis <- readRDS("qiu_reanalysis_seurat.rds")

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
                         ident.1 = "NS",
                         ident.2 = "S",
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

de_results_filtered <- map(de_results, function(df) {
  if (is.null(df)) return(NULL)
  dplyr::filter(df, p_val_adj < 0.05)
})
de_results_filtered_df <- bind_rows(de_results_filtered, .id = "celltype")

write.csv(de_results_df, file = "qiu_reanalysis_ns_vs_s_unfiltered.csv")
write.csv(de_results_filtered_df, file = "qiu_reanalysis_ns_vs_s_significant.csv")
