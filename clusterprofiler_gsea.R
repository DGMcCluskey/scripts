#pathway enrichment analysis using clusterprofiler
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/RIPs_Vincent_project/samples")
setwd("C:/Users/dan94/OneDrive - University College London/UCL_Senior_Research_Fellow/Tim_data/")

degs <- read.csv("deseq2_naive_cd4_mat_vs_fet_all_genes.csv")

names(degs)[1] <- "gene"


#if loading in from analysis script
degs <- filter(mono_ns_vs_hc_degs, rownames(mono_ns_vs_hc_degs) %in% ns_only$`setdiff(go.obj@listA, go.obj@intersection)`)

degs$gene <- rownames(degs)

#gene set enrichment (gse) is run on all genes
gene_list <- degs$avg_log2FC
names(gene_list) <- degs$gene
gene_list <- sort(gene_list, decreasing = T)

gene_list_entrez <- bitr(names(gene_list), 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

gene_list_kegg <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list_kegg) <- gene_list_entrez$ENTREZID

# Remove any NAs and duplicates
gene_list_kegg <- gene_list_kegg[!is.na(names(gene_list_kegg))]
gene_list_kegg <- gene_list_kegg[!duplicated(names(gene_list_kegg))]

gse <- gseGO(geneList=gene_list_kegg, OrgDb = org.Hs.eg.db)
enrich_list <- names(gene_list_kegg)
enrich <- enrichGO(enrich_list, OrgDb = org.Hs.eg.db)

dotplot(enrich)
enrich_results <- enrich@result
#rename to go/kegg/reactome to save separate analysis, then continue with plotting code
go <- gse
go
#if using go rather than kegg
gse <- clusterProfiler::simplify(gse, cutoff=0.7, by="p.adjust",
                 select_fun=min)

results <- gse@result
dotplot(gse)+facet_wrap(.~.sign)

# Get top 10 positive and top 10 negative enrichment scores
top_positive <- results[order(results$NES, decreasing = TRUE), ][1:10, ]
#fewer than 10 positive so filter out the negative ones that will be selected
top_positive <- filter(top_positive, NES > 0)
top_negative <- results[order(results$NES, decreasing = FALSE), ][1:10, ]

# Combine them
top_pathways <- rbind(top_positive, top_negative)


p <- ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES), fill = -log10(p.adjust)))+
  geom_col(colour = "black")+
  theme_bw()+scale_fill_viridis_c(option = "D", direction = -1)+
  ylab(NULL)+xlab("Normalised enrichment score")+
  theme(axis.text = element_text(size = 14, colour = "black"))
p

ggsave(plot = p, filename = "gsea_go_top_pathways_naive_cd4_mat_vs_fet.png", dpi = 500,
       height = 6, width = 12)

#running comparecluster, which allows comparison of multiple gene lists at once
non_survival <- read.csv("deseq2_nonclassical_non_survival_vs_HC_all_genes.csv")
names(non_survival)[1] <- "gene"
non_survival_list <- non_survival$stat
names(non_survival_list) <- non_survival$gene
non_survival_list <- sort(non_survival_list, decreasing = T)

non_survival_list_entrez <- bitr(names(non_survival_list), 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

non_survival_list_input <- non_survival_list[non_survival_list_entrez$SYMBOL]
names(non_survival_list_input) <- non_survival_list_entrez$ENTREZID

icu <- read.csv("deseq2_nonclassical_icu_vs_HC_all_genes.csv")
names(icu)[1] <- "gene"
icu_list <- icu$stat
names(icu_list) <- icu$gene
icu_list <- sort(icu_list, decreasing = T)

icu_list_entrez <- bitr(names(icu_list), 
                                 fromType = "SYMBOL", 
                                 toType = "ENTREZID", 
                                 OrgDb = org.Hs.eg.db)

icu_list_input <- icu_list[icu_list_entrez$SYMBOL]
names(icu_list_input) <- icu_list_entrez$ENTREZID

no_icu <- read.csv("deseq2_nonclassical_no_icu_vs_HC_all_genes.csv")
names(no_icu)[1] <- "gene"
no_icu_list <- no_icu$stat
names(no_icu_list) <- no_icu$gene
no_icu_list <- sort(no_icu_list, decreasing = T)

no_icu_list_entrez <- bitr(names(no_icu_list), 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)

no_icu_list_input <- no_icu_list[no_icu_list_entrez$SYMBOL]
names(no_icu_list_input) <- no_icu_list_entrez$ENTREZID

gene_lists <- list("non_survival" = non_survival_list_input, 
                  "icu"= icu_list_input, 
                  "no_icu" = no_icu_list_input)

gene_lists <- lapply(gene_lists, function(x) sort(x, decreasing = TRUE))

cc_results <- compareCluster(
  geneClusters = gene_lists,
  fun = "gseGO",
  OrgDb = org.Hs.eg.db
)

cc_results <- clusterProfiler::simplify(cc_results, cutoff=0.7, by="p.adjust",
                select_fun=min)


dotplot(cc_results)
comparison_results <- cc_results@compareClusterResult

top_pathways_per_cluster <- comparison_results %>%
  group_by(Cluster) %>%
  arrange(desc(enrichmentScore)) %>%  # or desc(NES) if using GSEA
  slice(c(1:10, if(n() >= 20) (n()-9):n() else 11:n())) %>%
  ungroup()

# Check how many pathways per cluster you got
table(top_pathways_per_cluster$Cluster)

top_pathways_per_cluster$Cluster <- factor(top_pathways_per_cluster$Cluster, levels = c("no_icu", "icu", "non_survival"))
p <- ggplot(top_pathways_per_cluster, aes(x = Cluster, y = reorder(Description, NES), fill = NES, size = -log10(p.adjust)))+
  geom_point(shape = 21)+
  theme_bw()+scale_fill_distiller(palette = "RdBu")+
  ylab(NULL)+xlab(NULL)+
  theme(axis.text = element_text(size = 14, colour = "black"))
p

ggsave(plot = p, filename = "GSEA_top_pathways_by_condition.png", dpi = 500,
       height = 8, width = 14, path = "../plots/")

#running clusterprofiler on monocyte DEGs from non_survival vs HC and survival vs HC
#i overlapped genes previously and want to look at three gene sets
#common genes, ns only and s only
mono_ns_vs_hc_degs$gene <- rownames(mono_ns_vs_hc_degs)
#ns_input <- filter(mono_ns_vs_hc_degs, mono_ns_vs_hc_degs$gene %in% ns_only$`setdiff(go.obj@listA, go.obj@intersection)`)

mono_s_vs_hc_degs$gene <- rownames(mono_s_vs_hc_degs)
#s_input <- filter(mono_s_vs_hc_degs, mono_s_vs_hc_degs$gene %in% s_only$`setdiff(go.obj@listB, go.obj@intersection)`)

mono_ns_vs_hc_degs <- read.csv("../plots/nonclassical_ns_vs_hc_degs.csv")
mono_s_vs_hc_degs <- read.csv("../plots/nonclassical_s_vs_hc_degs.csv")

mono_ns_vs_hc_degs$gene <- mono_ns_vs_hc_degs$X
mono_s_vs_hc_degs$gene <- mono_s_vs_hc_degs$X


shared <- intersect(mono_ns_vs_hc_degs$gene, mono_s_vs_hc_degs$gene)
shared_ns <- mono_ns_vs_hc_degs %>% 
  filter(gene %in% shared)
shared_s <- mono_s_vs_hc_degs %>% 
  filter(gene %in% shared)
shared_merged <- merge(shared_ns, shared_s, by = "gene", suffixes = c("_ns", "_s"))
shared_merged$avg_log2FC_combined <- (shared_merged$avg_log2FC_ns + shared_merged$avg_log2FC_s) / 2

ns_genes <- mono_ns_vs_hc_degs$avg_log2FC
names(ns_genes) <- mono_ns_vs_hc_degs$gene
ns_entrez <- bitr(names(ns_genes), 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)

ns_genes <- ns_genes[ns_entrez$SYMBOL]
names(ns_genes) <- ns_entrez$ENTREZID


s_genes <- mono_s_vs_hc_degs$avg_log2FC
names(s_genes) <- mono_s_vs_hc_degs$gene
s_entrez <- bitr(names(s_genes), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

s_genes <- s_genes[s_entrez$SYMBOL]
names(s_genes) <- s_entrez$ENTREZID

shared_genes <- shared_merged$avg_log2FC_combined
names(shared_genes) <- shared_merged$gene
shared_entrez <- bitr(names(shared_genes), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

shared_genes <- shared_genes[shared_entrez$SYMBOL]
names(shared_genes) <- shared_entrez$ENTREZID


gene_list <- list(non_survival = ns_genes, survival = s_genes)
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
  slice_head(n = 15) %>%
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
    cluster_combo == "non_survival_survival" ~ "shared",
    cluster_combo == "survival" ~ "survival_only",
    cluster_combo == "non_survival" ~ "non_survival_only",
    TRUE ~ "other"
  ))

# Merge back into final_pathways
final_pathways <- final_pathways %>%
  left_join(cluster_counts, by = "Description")



final_pathways$group <- factor(final_pathways$group, levels = c("survival_only", "shared", "non_survival_only"))
p <- ggplot(final_pathways, aes(x = Cluster, y = reorder(Description, NES), fill = NES, size = -log10(p.adjust)))+
  geom_point(shape = 21)+scale_size_continuous(range = c(3,8))+
  theme_bw()+scale_fill_distiller(palette = "RdBu")+
  ylab(NULL)+xlab(NULL)+
  theme(axis.text = element_text(size = 14, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(group ~ ., space = "free", scales = "free")+
  theme(axis.text.y = element_text(size = 14, colour = "black"))
p

ggsave(plot = p, filename = "GSEA_GO_nonclassical_monocytes_qiu_survival_non_survival.png",
       dpi = 500, height = 12, width = 10, path = "../plots/")

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

