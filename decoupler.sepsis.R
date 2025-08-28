#running decoupler on sepsis data to find out transcription factor activity
library(decoupleR)
library(Seurat)
setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/RIPs_Vincent_project/samples")

#load in sepsis data and pull out only sepsis and HC samples
sepsis_final <- readRDS(file = "sepsis_final.rds")
Idents(sepsis_final) <- "condition2"
table(sepsis_final@active.ident)
sepsis_only <- subset(sepsis_final, idents = c("HC", "Sepsis"))

#get the prior knowledge data of TFs and their targets
net <- decoupleR::get_collectri(organism = "human", split_complexes = F)

#get normalised log counts
#split by cell type as this is memory instensive (also makes sense to do it on one cluster at a time)
Idents(qiu_data) <- "labels"
mac_mono <- subset(qiu_data, idents = c("Mac/mono"))

mat <- GetAssayData(mac_mono, slot = "data")
rownames(mat)


acts <- run_ulm(mat = mat, network = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)



