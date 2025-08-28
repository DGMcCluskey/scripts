#glmer on poms data to assess changes in abundance over time in foetal and maternal cells
#code adapted from rik lindeboom's scripts from nature covid challenge paper

library(zellkonverter)
library(SingleCellExperiment)
library(lme4)
library(dplyr)
library(numDeriv)
library(msigdbr)
library(Rsamtools)
library(metafor)
library(tidyverse)
library(DESeq2)
library(scCustomize)
library(Seurat)


setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Tim_data")

poms <- readH5AD("poms.annotated.dan.h5ad")

poms_seurat <- as.Seurat(poms)

poms.seurat$gestation_bin <- ifelse(poms.seurat$gestation_week > 28, "moderate", "extreme")

DimPlot(poms.seurat, group.by = "basic_labels")

FeaturePlot(poms.seurat, features = c("PRDM1", "IRF4", "XBP1",
                                      "CD27", "CD38", "CD79A"))


poms <- poms.seurat


table(poms$basic_labels, poms$patient_nice)
poms$cell_id <- paste0("Cell_", seq_len(ncol(poms)))

poms.df <- data.frame(cell_id = poms$cell_id, sample_id = poms$patient_nice,
                      fetmat = poms$fetmat, time = poms$gestation_week,
                      batch = poms$poolID, labels = poms$basic_labels,
                      gestation_bin = poms$gestaion_bin)

table(poms.df$labels)

poms.df$labels <- factor(poms.df$labels, levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                         "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                         "Mono/DC", "pDC", "HSC"))


poms.df$batch <- factor(poms.df$batch, levels = c("sc12", "sc13", "sc14", "sc15", "sc16", "sc17", "sc18", "sc19"),
                        labels = c("sc12_13", "sc12_13", "sc14_15", "sc14_15", "sc16_17", "sc16_17", "sc18_19", "sc18_19"))

table(poms.df$sample_id, poms.df$batch)

poms.df$bin.time <- factor(poms.df$time, levels = c(23.4, 24.1, 25,25.3,25.4,27.3,29.6,32.7,34,35.9),
                          labels = c("23-24(n=2)", "23-24(n=2)", "25-26(n=3)", "25-26(n=3)", "25-26(n=3)","27-30(n=2)", "27-30(n=2)",
                                     "32-36(n=3)", "32-36(n=3)", "32-36(n=3)"))


poms.df$group <- ifelse(grepl("Maternal", poms.df$sample_id), "Maternal",
                        ifelse(poms.df$sample_id == "Fetal_9", "Covid_vaccine",
                               "Fetal"))

poms.df$group <- factor(poms.df$group, levels = c("Fetal", "Covid_vaccine", "Maternal"))


umap_cols <- c("steelblue2", "blue", "cyan4", "darkorchid", "darkorchid4", "mediumvioletred", "hotpink",
               "mediumseagreen", "forestgreen", "green2", "red", "indianred2", "darkred", "orange", "grey50")


poms.df <- poms.df %>% group_by(sample_id, fetmat, group, bin.time, batch, labels, time) %>% tally()
p <- ggplot(poms.df, aes(x = sample_id, y = n, fill = labels))+geom_col(stat = "identity", position = "fill", colour = "white")+theme_bw()+
  facet_grid(~group, scales = "free_x", space = "free")+scale_y_continuous(expand =(c(0,0)))+
  scale_fill_manual(values = umap_cols)+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  ylab("Proportion of cells")+xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave(plot = p, filename = "poms.stacked.bar.proportions.png", dpi = 500,
       height = 6, width = 14)

p <- ggplot(poms.df, aes(x = sample_id, y = n, fill = labels))+geom_col(stat = "identity", position = "dodge", colour = "white")+theme_bw()+
  facet_grid(~group, scales = "free_x", space = "free")+scale_y_continuous(expand =(c(0,0)))+
  scale_fill_manual(values = umap_cols)+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  ylab("Proportion of cells")+xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

p


#now just b cell subsets
p <- ggplot(filter(poms.df, labels %in% c("Transitional_B", "Naive_B", "Memory_B")), aes(x = sample_id, y = n, fill = labels))+geom_col(stat = "identity", position = "fill", colour = "white")+theme_bw()+
  facet_grid(~group, scales = "free_x", space = "free")+scale_y_continuous(expand =(c(0,0)))+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid", "goldenrod2"))+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  ylab("Proportion of cells")+xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("B cell subsets")

p

ggsave(plot = p, filename = "poms.b.cells.stacked.bar.proportions.png", dpi = 500,
       height = 6, width = 14)

#get combined proportion of non-transitional b cells and then plot these
bcell_df <- poms.df %>% filter(fetmat %in% "Fetal")

# Calculate proportions per sample
bcell_props <- bcell_df %>%
  group_by(sample_id) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  mutate(proportion = (n / total_n)*100)

#keep only non-transitional, add proportions and plot
input <- filter(bcell_props, labels %in% c("Naive_B", "Memory_B"))
input <- input[, c(1,3,6,10)]
input <- pivot_wider(input, names_from = "labels", values_from = "proportion")
input$prop <- input$Naive_B+input$Memory_B

# Extract only one cell type (e.g., "Memory_B")

p <- ggplot(input, aes(x = reorder(sample_id, prop), y = prop))+geom_col(position = position_dodge(), fill = "indianred2")+
  theme_bw()+xlab(NULL)+ylab("Proportion differentiated (naive+memory) B cells \n of total PBMCs")+theme(axis.text = element_text(size = 14, colour = "black"))+
  theme(axis.title.x = element_text(size = 16, colour = "black"))+coord_flip()
p

ggsave(plot = p, filename = "proportion.non.transitional.b.png", dpi = 500,
       height = 3, width = 6)
#proportions for boxplot
poms.df <- pivot_wider(poms.df, names_from = "labels", values_from = "n")
poms.df[is.na(poms.df)] <- 0

poms.df$total <- rowSums(poms.df[6:20])

poms.df <- pivot_longer(poms.df, names_to = "label", values_to = "n", 6:20) 

poms.df$proportion <- (poms.df$n/poms.df$total)*100

p <- ggplot(poms.df, aes(x = fetmat, y = proportion, fill = fetmat))+geom_boxplot(outlier.shape = NA)+
  theme_bw()+facet_wrap(~label, scales = "free_y", ncol = 5)+geom_point(shape = 21)+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen"))+
  xlab(NULL)+ylab("Proportion of cells")+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  theme(axis.text = element_text(size = 12, colour = "black"))
p

ggsave(plot = p, filename = "proportion.boxplots.png", dpi = 500,
       height = 6, width = 12)

#longitudinal proportions

#do this with Fetal_9 removed as they are clear outlier (and there was a covid infection in the mother)

p <- ggplot(filter(poms.df, !sample_id %in% "Fetal_9"), aes(x = time, y = proportion, colour = fetmat, linetype = fetmat))+geom_smooth(se = F)+
  geom_point()+facet_wrap(~label, scales = "free_y", ncol = 5)+theme_bw()+scale_colour_manual(values = c("darkorchid2", "mediumseagreen"))+
  xlab("Gestational weeks")+ylab("Proportion of cells")+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  theme(axis.text = element_text(size = 12, colour = "black"))
p

ggsave(plot = p, filename = "proportions.longitudinal.png", dpi = 500,
       height = 6, width = 18)


fet.df <- filter(poms.df, fetmat == "Fetal")

fet.df <- fet.df %>% group_by(sample_id, fetmat, time, batch, labels) %>% tally()

fet.df$sample_id <- droplevels(fet.df$sample_id)

table(fet.df$sample_id, fet.df$labels)

insert <- data.frame(sample_id	= c("Fetal_1", "Fetal_9", "Fetal_12"),
                     fetmat = c("Fetal", "Fetal", "Fetal"),
                     time = c(29.6, 27.3, 23.4),
                     batch = c("sc12_13", "sc14_15", "sc18_19"),
                     labels = c("CD16+_mono", "CD16+_mono", "CD16+_mono"), n = c(0,0,0))
                     
                     
fet.df <- rbind(fet.df, insert)

fet.df <- pivot_wider(fet.df, names_from = "labels", values_from = "n")
fet.df$total <- rowSums(fet.df[, 5:19])
fet.df <- pivot_longer(fet.df, names_to = "labels", values_to = "number", 5:19)
fet.df$prop <- (fet.df$number/fet.df$total)*100

fet.df$time

fet.df$bin.time <- factor(fet.df$time, levels = c(23.4, 24.1, 25,25.3,25.4,27.3,29.6,32.7,34,35.9),
                          labels = c("23-24(n=2)", "23-24(n=2)", "25-26(n=3)", "25-26(n=3)", "25-26(n=3)","27-30(n=2)", "27-30(n=2)",
                                     "32-36(n=3)", "32-36(n=3)", "32-36(n=3)"))

table(fet.df$bin.time)


#EXTREME vs MODERATE PRETERM COMPARISON============================================================================================

poms.fet <- filter(poms.df, fetmat %in% "Fetal")
poms.fet <- filter(poms.fet, !sample_id %in% "Fetal_9")


poms.fet <- poms.fet %>% group_by(sample_id, fetmat, gestation_bin, batch, labels, time) %>% tally()
p <- ggplot(poms.fet, aes(x = sample_id, y = n, fill = labels))+geom_col(stat = "identity", position = "fill", colour = "white")+theme_bw()+
  facet_grid(~gestation_bin, scales = "free_x", space = "free")+scale_y_continuous(expand =(c(0,0)))+
  scale_fill_manual(values = umap_cols)+
  theme(strip.text = element_text(size = 12, colour = "black"))+
  ylab("Proportion of cells")+xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave(plot = p, filename = "poms.gestation.stacked.bar.proportions.png", dpi = 500,
       height = 6, width = 14)

#boxplots
poms.fet <- pivot_wider(poms.fet, names_from = "labels", values_from = "n")
poms.fet[is.na(poms.fet)] <- 0

poms.fet$total <- rowSums(poms.fet[6:20])

poms.fet <- pivot_longer(poms.fet, names_to = "label", values_to = "n", 6:20) 

poms.fet$proportion <- (poms.fet$n/poms.fet$total)*100


poms.fet$label <- factor(poms.fet$label, levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                                                    "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                                                    "Mono/DC", "pDC", "HSC"))

p <- ggplot(poms.fet, aes(x = gestation_bin, y = proportion, fill = gestation_bin))+
  geom_boxplot(outlier.shape = NA, alpha = 0.4, colour = "black")+
  geom_jitter(width = 0.1, shape = 21, colour = "black", size = 2)+
  theme_bw()+
  facet_wrap(~label, scales = "free_y", ncol = 5)+
  scale_fill_manual(values = c("darkorchid", "mediumseagreen"))+
  xlab(NULL)+
  theme(axis.text = element_text(size = 12, colour = "black"))
p

results <- poms.fet %>%
  group_by(label) %>%
  summarise(
    p_value = wilcox.test(proportion ~ gestation_bin)$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)

#before running glmer lets just see what the proportions look like


fet.df.grouped <- fet.df %>% group_by(labels, bin.time) %>% summarise(mean = mean(prop))
ggplot(fet.df.grouped, aes(x = bin.time, y = mean, fill = labels))+geom_col()+theme_bw()

fet.df.grouped <- fet.df %>% group_by(labels, bin.time, sample_id) %>% summarise(mean = mean(prop))
ggplot(fet.df.grouped, aes(x = sample_id, y = mean, fill = labels))+geom_col()+theme_bw()+
  facet_wrap(~bin.time, scales = "free_x", ncol = 4)


ncells <- length(unique(fet.df$labels))

#res.prop <- glmer(I(c(prop)) ~ (1|labels)+(1|sample_id)+(1|bin.time)+(1|batch)+(1|sample_id:labels)+(1|bin.time:labels),
                  #data = fet.df,
                  #family = poisson, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))

#my version
res.prop <- glmer(number ~ (1|labels)+(1|sample_id)+(1|batch)+(1|bin.time)+(1|sample_id:labels)+(1|bin.time:labels), 
                  offset = log(total),
                  data = fet.df,
                  family = poisson, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))

summary(res.prop)

# standard errors of standard deviations (square root of the variance parameters)
devfun = update(res.prop, devFunOnly=T)
pars = getME(res.prop, c("theta","fixef"))
hess = hessian(devfun, unlist(pars))
sdse.prop = data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))

# posterior means and their standard deviations
res.prop.ranef = ranef(res.prop)


sdse.prop$effects <- rownames(sdse.prop)
ggplot(sdse.prop, aes(x = sd, y = effects))+geom_point()+
  geom_errorbar(aes(xmin = sd-se, xmax = sd+se), width = 0.1)+theme_bw()+
  geom_vline(xintercept = 0, linetype = "dashed")

getCondVal <- function(res.prop.ranef, id, ncells, nfactors=2){
  tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
  if(length(grep(":",tmp$grp))==0){
    cnam = matrix(as.character(tmp$term),ncells)[1,]
    rnam = matrix(as.character(tmp$grp),ncells)[,1] 
  }else if(nfactors==2){
    cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
  }else if(nfactors==3){
    cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
    cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
    cnam = paste(cnam1,cnam2,sep=":")
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
  }
  condval = matrix(tmp$condval,ncells)
  condsd  = matrix(tmp$condsd, ncells)
  rownames(condval)=rownames(condsd)=rnam
  colnames(condval)=colnames(condsd)=cnam
  lfsr = pnorm(condval,0,condsd)
  lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
  list(condval=condval, lfsr=lfsr)
}


res1 <- getCondVal(res.prop.ranef,"bin.time:labels",ncells)
res2 <- getCondVal(res.prop.ranef,"sample_id:labels",ncells)
postmean <- as.data.frame(cbind(res1[[1]]))  # condval matrices
lfsr <- as.data.frame(cbind(res1[[2]]))

postmean$labels <- rownames(postmean)
postmean <- pivot_longer(postmean, names_to = "gestation_weeks", values_to = "log2FC", 1:4)
lfsr$labels <- rownames(lfsr)
lfsr <- pivot_longer(lfsr, names_to = "gestation_weeks", values_to = "LTSR", 1:4)

lfsr <- lfsr[, 3]



abundance.data <- cbind(postmean, lfsr)

abundance.data$LTSR <- 1-abundance.data$LTSR

abundance.data <- mutate(abundance.data,
                         "LTSR>0.9" = LTSR > 0.9)

abundance.data$labels <- factor(abundance.data$labels, levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                                                    "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                                                    "Mono/DC", "pDC", "HSC"))


p <- ggplot(abundance.data, aes(x = gestation_weeks, y = labels, fill = log2FC, colour = LTSR>0.9, 
                                size = LTSR))+geom_point(shape = 21)+
  theme_bw()+scale_fill_distiller(palette = "RdBu")+scale_colour_manual(values = c("black", "red"))+
  theme(panel.grid = element_blank())+
   scale_size_continuous(
    breaks = c(0.5, 0.7, 0.9),  # Set your LTSR thresholds (these are your breaks)
    range = c(0.1, 10), limits = c(0,1))+
  theme(legend.position = "bottom")+theme(axis.text = element_text(size = 12, color = "black"))+
  ylab(NULL)+xlab("Gestational weeks")+scale_y_discrete(limits = rev)
p

#get average proportions to add as 2nd plot
fet.df.grouped <- fet.df %>% group_by(labels) %>% summarise(mean = mean(prop))

fet.df.grouped$labels <- factor(fet.df.grouped$labels, 
                                levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                                            "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                                            "Mono/DC", "pDC", "HSC"))
p2 <- ggplot(fet.df.grouped, aes(x = mean, y = labels))+geom_col(fill = "black")+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,65))+theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  xlab("Average proportion \n of PBMCs")+ylab(NULL)+
  theme(axis.text.y = element_blank())+scale_y_discrete(limits = rev)
p2

combo <- plot_grid(p,p2,rel_widths = c(3,1), align = "h")
combo

ggsave(plot = combo,filename = "glmer.differential.abundance.fetal.dotplot.png", dpi = 500,
       height = 6, width = 12)


#regular propotions boxplot
p <- Proportion_Plot(poms_seurat, group_by_var = "basic_labels", split.by = "patient_nice")+scale_fill_manual(values = umap_cols)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p

ggsave(plot = p, filename = "proportions.barplots.png", dpi = 500,
       height = 6, width = 12)


#DIFFERENTIAL EXPRESSION ANALYSIS===========================================================================================================================
#i want to run DEA on fetal vs maternal to start with, on pseudobulked samples
#attempting to run DEA on all cell subsets
#i need poms to be a seurat object too
poms.seurat <- as.Seurat(poms, counts = "counts", data = "logcounts")

poms.seurat$batch <- factor(poms.seurat$poolID, levels = c("sc12", "sc13", "sc14", "sc15", "sc16", "sc17", "sc18", "sc19"),
                        labels = c("sc12_13", "sc12_13", "sc14_15", "sc14_15", "sc16_17", "sc16_17", "sc18_19", "sc18_19"))

# Pseudobulk using Seurat's built-in function
agg <- AggregateExpression(poms.seurat, group.by = c("patient_nice", "basic_labels", "fetmat", "batch", "gestation_week"), slot = "counts", return.seurat = T)

table(agg$basic_labels)

DefaultAssay(agg) <- "originalexp"

agg@assays$originalexp$counts <- round(GetAssayData(agg, slot = "counts"))

celltypes <- unique(agg$basic_labels)

#runnig deseq2 outside of seurat to allow for covariates

# Store DE results for each cell type
de_results <- map(celltypes, function(ct) {
  message("Processing: ", ct)
  
  # Subset the pseudobulked Seurat object
  subset_obj <- subset(agg, subset = basic_labels == ct)
  
  # Extract the counts matrix
  counts <- GetAssayData(subset_obj, assay = "originalexp", slot = "counts")
  
  # Get sample metadata
  meta <- subset_obj@meta.data
  
  # Check condition balance (important!)
  subtype_counts <- table(meta$fetmat)
  if (length(subtype_counts) < 2 || any(subtype_counts < 2)) {
    message(paste("Skipping", ct, "due to insufficient cells per condition"))
    return(NULL)
  }
  
  # Prepare DESeq2 input
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ batch + gestation_week + fetmat  # Include covariates
  )
  
  # Filter out lowly expressed genes (optional but recommended)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results: fetal vs maternal
  res <- results(dds, contrast = c("fetmat", "Fetal", "Maternal"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$celltype <- ct
  
  return(res_df)
})

# Combine all celltype DE results into one dataframe
de_results_df <- bind_rows(de_results)

# Combine all results
names(de_results) <- celltypes[!map_lgl(de_results, is.null)]
#if some celltypes don't have degs due to lack of numbers then this will mess order above up
#thefore need to manually alter them
#names(de_results)[9] <- "DC"

de_results_df <- bind_rows(de_results, .id = "celltype")

de_results_filtered <- map(de_results, function(df) {
  if (is.null(df)) return(NULL)
  dplyr::filter(df, padj < 0.05)
})
de_results_filtered_df <- bind_rows(de_results_filtered, .id = "celltype")


write.csv(de_results_df, file = "poms_fetal_vs_maternal_unfiltered.csv")
write.csv(de_results_filtered_df, file = "poms_fetal_vs_maternal_significant.csv")

de_results_filtered_df <- read.csv("poms_fetal_vs_maternal_significant.csv")


#plot number of degs
deg_numbers <- de_results_filtered_df %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  group_by(celltype) %>%
  summarise(n_DEGs = n())

umap_cols <- brewer.pal(name = "Paired", n = 12)
umap_cols <- c(umap_cols,"grey50", "cyan4", "black")

#or my own palette

umap_cols <- c("steelblue2", "blue", "cyan4", "darkorchid", "darkorchid4", "mediumvioletred", "hotpink",
               "mediumseagreen", "forestgreen", "green2", "red", "indianred2", "darkred", "orange", "grey50")

deg_numbers$celltype <- factor(deg_numbers$celltype, levels = c("Naive-CD4+T", "Memory-CD4+T", "Treg", "Naive-CD8+T", "Memory-CD8+T", 
                                                                "MAIT", "NK", "Transitional-B", "Naive-B", "Memory-B", "CD14+-mono", "CD16+-mono",
                                                                "Mono/DC", "pDC", "HSC"))


p <- ggplot(deg_numbers, aes(y = celltype, x = n_DEGs, fill = celltype))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = umap_cols)+
  ggtitle("Number of DEGs Fetal vs Maternal (log2FC > |1| & FDR < 0.05)")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)
p

ggsave(plot = p, filename = "DEG_numbers.png", dpi = 500,
       height = 5, width = 7)

poms.seurat$basic_labels <- factor(poms.seurat$basic_labels, levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                                                                        "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                                                                        "Mono/DC", "pDC", "HSC"))

p <- DimPlot(poms.seurat, group.by = "basic_labels", label = T, repel = T, label.box = T, label.color = "white")+
  scale_colour_manual(values = umap_cols)+
  scale_fill_manual(values = umap_cols)+
  theme(legend.position = "none")
p

ggsave(plot = p, filename = "UMAP.annotated.png", dpi = 500,
       height = 5, width = 7)
#only degs above or below FC threshold - and then split by direction

#separated by up and down
deg_numbers <- de_results_filtered_df %>%
  mutate(Direction = case_when(
    log2FoldChange > 1  ~ "Up",
    log2FoldChange < -1 ~ "Down",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Direction)) %>%
  group_by(celltype, Direction) %>%
  summarise(n_DEGs = n(), .groups = "drop")
deg_numbers


umap_cols <- brewer.pal(name = "Paired", n = 12)
umap_cols <- c(umap_cols,"grey50", "black", "cyan4")


p <- ggplot(deg_numbers, aes(y = celltype, x = n_DEGs, fill = Direction))+
  geom_col(position = position_dodge(), colour = "black")+theme_bw()+
  scale_fill_manual(values = c("darkorchid2", "mediumseagreen"))+
  ggtitle("Number of DEGs Fetal vs Maternal (Fold change > |2|)")+
  ylab(NULL)+xlab("Number of DEGs")+
  theme(axis.text = element_text(size = 12, colour = "black"))+
  scale_y_discrete(limits = rev)
p


deg_numbers
#add deg number as a metadata column for visualisation
poms.seurat$degs <- factor(poms.seurat$basic_labels, levels = c("Naive_CD4+T", "Memory_CD4+T", "Treg", "Naive_CD8+T", "Memory_CD8+T", 
                                                                "MAIT", "NK", "Transitional_B", "Naive_B", "Memory_B", "CD14+_mono", "CD16+_mono",
                                                                "Mono/DC", "pDC", "HSC"),
                           labels = c(2140,308,596,1468,567,81,339,34,72,29,544,24,399,16,86))
table(poms.seurat$degs)

poms.seurat$degs <- as.numeric(as.character(poms.seurat$degs))
poms.seurat$degs_log <- log10(poms.seurat$degs)


p1 <- FeaturePlot(poms.seurat, features = "degs")+scale_colour_viridis(option = "F", direction = -1)+ggtitle("Number of DEGs")
p2 <- FeaturePlot(poms.seurat, features = "degs_log")+scale_colour_viridis(option = "F", direction = -1)+ggtitle("Number of DEGs (log10)")

combo <- p1+p2+plot_layout(ncol = 2)

ggsave(plot = combo, filename = "UMAPs.nDEGs.png", dpi = 500,
       height = 5, width = 12)



#next steps (or at least plan for when full data arrives) - run DEA on binned timepoints for fetal. compare one timepoint to all others, systematically
#then look at overlapping genes to see if there is a set of unique genes at each timepoint/bin

#CELLCHAT==============================================================================================
library(CellChat)
#cell-cell interactions
#first using cellchat
Idents(poms.seurat) <- "patient_nice"
poms.seurat <- subset(poms.seurat, idents = "Fetal_9", invert = T)
Idents(poms.seurat) <- "fetmat"
table(poms.seurat@active.ident)
fet <- subset(poms.seurat, idents = "Fetal")
mat <- subset(poms.seurat, idents = "Maternal")

fet@assays

fet.chat <- createCellChat(fet, group.by = "basic_labels", assay = "originalexp")
mat.chat <- createCellChat(mat, group.by = "basic_labels", assay = "originalexp")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#use of cellchatdb except for non-protein signalling
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
fet.chat@DB <- CellChatDB.use
mat.chat@DB <- CellChatDB.use
#HC.chat@DB <- CellChatDB.use

fet.chat <- subsetData(fet.chat) # This step is necessary even if using the whole database
mat.chat <- subsetData(mat.chat) # This step is necessary even if using the whole database
#HC.chat <- subsetData(HC.chat)


future::plan("multisession", workers = 4) # do parallel
fet.chat <- identifyOverExpressedGenes(fet.chat)
fet.chat <- identifyOverExpressedInteractions(fet.chat)
mat.chat <- identifyOverExpressedGenes(mat.chat)
mat.chat <- identifyOverExpressedInteractions(mat.chat)
#HC.chat <- identifyOverExpressedGenes(HC.chat)
#HC.chat <- identifyOverExpressedInteractions(HC.chat)

fet.chat <- subsetData(fet.chat)
mat.chat <- subsetData(mat.chat)

#fet.chat <- projectData(fet.chat, PPI.human)
fet.chat <- computeCommunProb(fet.chat, type = "triMean", population.size = T)
mat.chat <- computeCommunProb(mat.chat, type = "triMean", population.size = T)
#HC.chat <- computeCommunProb(HC.chat, type = "triMean")


fet.chat <- filterCommunication(fet.chat, min.cells = 10)
mat.chat <- filterCommunication(mat.chat, min.cells = 10)
#HCE.chat <- filterCommunication(HC.chat, min.cells = 10)

fet.chat <- computeCommunProbPathway(fet.chat)
mat.chat <- computeCommunProbPathway(mat.chat)
#HC.chat <- computeCommunProbPathway(HC.chat)

fet.chat <- aggregateNet(fet.chat)
mat.chat <- aggregateNet(mat.chat)
#HC.chat <- aggregateNet(HC.chat)

options(future.seed = TRUE)

fet.chat <- netAnalysis_computeCentrality(fet.chat, slot.name = "netP")
mat.chat <- netAnalysis_computeCentrality(mat.chat, slot.name = "netP")
#HC.chat <- netAnalysis_computeCentrality(HC.chat, slot.name = "netP")

#saving each cellchat object so i dont have to load again
saveRDS(fet.chat, file = "fet.cellchat.rds")
saveRDS(mat.chat, file = "mat.cellchat.rds")
#saveRDS(HC.chat, file = "HC.cellchat.rds")

fet.chat <- readRDS(file = "fet.cellchat.rds")
mat.chat <- readRDS(file = "mat.cellchat.rds")


object.list <- list(fet = fet.chat, mat = mat.chat)
poms.chat <- mergeCellChat(object.list, add.names = names(object.list))


netVisual_diffInteraction(poms.chat, comparison = c(1,2), vertex.size.max = 10, edge.width.max = 1,
                          label.edge = F)

p1 <- compareInteractions(poms.chat, show.legend = F, group = c("fet", "mat"))+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+#scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 16, colour = "black"))
p2 <- compareInteractions(poms.chat, show.legend = F, group = c("fet", "mat"), measure = "weight")+
  scale_fill_manual(values = c("mediumseagreen", "darkorchid"))+#scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 16, colour = "black"))

p <- p1+p2

p

ggsave(plot = p, filename = "fet_mat_interactions.number.comparison.png",
       dpi = 500,
       height = 5, width = 6)

netVisual_heatmap(poms.chat, comparison = c(1,2))
netVisual_heatmap(poms.chat, measure = "weight")

rankNet(poms.chat, mode = "comparison", comparison = c(1,2), measure = "weight", stacked = F, do.stat = T)
pathways <- rankNet(poms.chat, mode = "comparison", comparison = c(1,2), measure = "weight", stacked = T, do.stat = T)
pathways
pathways <- pathways$data
pathways <- filter(pathways, pvalues < 0.05)

diff.pathways <- unique(pathways$name)

p <- ggplot(pathways, aes(x = contribution.scaled, y = name, fill = group))+
  geom_col(position = "fill", width = 0.5)+
  theme_bw()+scale_x_continuous(expand = c(0,0))+geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30")+
  scale_fill_manual(values = c("darkorchid", "mediumseagreen"))+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  ylab("Pathway")+xlab("Relative signalling strength")

p

ggsave(plot = p, filename = "fet_matrelative.pathways.png", dpi = 500, height = 8, width = 8)

a <- netAnalysis_signalingRole_scatter(fet.chat)
a <- a$data
b <- netAnalysis_signalingRole_scatter(mat.chat)
b <- b$data

colnames(a) <- paste(colnames(a), "_fet")
colnames(b) <- paste(colnames(b), "_mat")

c <- cbind(a,b)

c$outgoing.diff <- c$`x _fet`-c$`x _mat`
c$incoming.diff <- c$`y _fet`-c$`y _mat`
c$count.diff <- c$`Count _fet`-c$`Count _mat`

p3 <- ggplot(c , aes(x = outgoing.diff, y = incoming.diff, fill = `labels _fet`, label = `labels _mat`, colour = `labels _fet`))+
  geom_text_repel()+
  geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(size = 5)+theme_bw()+scale_fill_manual(values = umap_cols)+
  scale_colour_manual(values = umap_cols)+theme(legend.position = "none")

p3

ggsave(plot = p3, filename = "ingoing.outgoing.change.png", dpi = 500,
       height = 5, width = 7)


p1 <- netAnalysis_signalingRole_scatter(fet.chat)+ylim(0,0.3)+xlim(0,0.3)
p2 <- netAnalysis_signalingRole_scatter(mat.chat)+ylim(0,0.3)+xlim(0,0.3)
combo <- p1+p2+p3+plot_layout(ncol = 3, widths = c(0.5,0.5,1))

ggsave(plot = combo, filename = "celltype.interaction.strength.changes.png", dpi = 500,
       height = 5, width = 14)


netVisual_chord_cell(fet.chat, signaling = "CCL", title.name = "Fetal CCL signalling")
netVisual_chord_cell(mat.chat, signaling = "CCL", title.name = "Meternal CCL signalling")

p1 <- netAnalysis_contribution(fet.chat, signaling = "CCL", title = "CCL ligand-receptor contribution (fetal)")
p2 <- netAnalysis_contribution(mat.chat, signaling = "CCL", title = "CCL ligand-receptor contribution (maternal)")

combo <- p1+p2+plot_layout(ncol = 1)

ggsave(plot = combo, filename = "CCL.LR.contribution.png", dpi = 500,
       height = 4, width = 6)

netVisual_chord_gene(fet.chat, signaling = "CCL", title.name = "Fetal CCL signaling at ligand-receptor level")

interactions$interaction <- paste0(interactions$source, "_", interactions$target, "_", interactions$interaction_name)


