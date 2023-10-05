library(MAST)


################################################################################
# CREATE RDS OBJECT FOR INDIVIDUAL TUMORS

gbm.seurat <- readRDS("DepMap.Scored_scGBM_2023.02.21.RDS")
Idents(gbm.seurat) <- gbm.seurat$NPvNon
gbm.np <- subset(gbm.seurat, idents = "Neoplastic")

saveRDS(
  gbm.np, 
  "Output/Rdata/scRNA/scGBM_NP_2023.10.02.rds"
)

# isolate gbm tumors as list object
tum <- as.character(unique(gbm.np$orig.ident))

Idents(gbm.np) <- gbm.np$orig.ident
gbm.xPT.list <- list()
for (i in 1:length(tum)) {
  np <- subset(gbm.np, idents = tum[i])
  gbm.xPT.list[[length(gbm.xPT.list) + 1]] <- assign(tum[i], np)
  rm(np)
}
gc()

names(gbm.xPT.list) <- tum

saveRDS(
  gbm.xPT.list,
  "Output/Rdata/scRNA/scGBM_NP_xPT_list_2023.09.19.rds"
)


################################################################################
# ADD CLUSTER DATA TO SCRNA-SEQ DATA

names(gbm.xPT.list)
# [1] "GBM21" "GBM41" "GBM47" "GBM49" "GBM51" "GBM53"

clust.data.1 <- read.csv(
  "Output/Rdata/jaccard/03_GBM21_cluster.data_2023.09.18.csv", 
  row.names = 1
)
clust.data.2 <- read.csv(
  "Output/Rdata/jaccard/03_GBM41_cluster.data_2023.09.18.csv", 
  row.names = 1
)
clust.data.3 <- read.csv(
  "Output/Rdata/jaccard/03_GBM47_cluster.data_2023.09.18.csv", 
  row.names = 1
)
clust.data.4 <- read.csv(
  "Output/Rdata/jaccard/03_GBM49_cluster.data_2023.09.18.csv", 
  row.names = 1
)
clust.data.5 <- read.csv(
  "Output/Rdata/jaccard/03_GBM51_cluster.data_2023.09.18.csv", 
  row.names = 1
)
clust.data.6 <- read.csv(
  "Output/Rdata/jaccard/03_GBM53_cluster.data_2023.09.18.csv", 
  row.names = 1
)

clust.data <- rbind(
  clust.data.1, 
  clust.data.2, 
  clust.data.3, 
  clust.data.4, 
  clust.data.5, 
  clust.data.6
)

for (i in 1:length(gbm.xPT.list)) {
  gbm.xPT.list[[i]] <- AddMetaData(
    gbm.xPT.list[[i]], 
    clust.data$cluster[which(clust.data$tumor == names(gbm.xPT.list)[i])], 
    col.name = "jaccard.cluster"
  )
}

# get essential genes that had scaled scRNA-seq data
# jaccard.input <- read.csv(
#   "Output/Rdata/jaccard/00_binary.expr.mtx_2023.09.12.csv",
#   row.names = 1
# )
ess.genes <- rownames(jaccard.input)
ess <- read.csv("Output/Rdata/04_gbm.killing_2023.09.07.csv")

# ################################################################################
# ## PERFORM DIFFERENTIAL EXPRESSION BETWEEN CLUSTERS, ESSENTIAL GENES ONLY
# for (i in 1:length(tnbc.xPT.list)) {
#   
#   Idents(tnbc.xPT.list[[i]]) <- tnbc.xPT.list[[i]]@meta.data$jaccard.cluster
#   
#   deg <- FindAllMarkers(
#     tnbc.xPT.list[[i]],
#     assay = "RNA", 
#     test.use = "MAST", 
#     features = ess.genes,
#     only.pos = FALSE
#   )
#   
#   write.csv(
#     deg, 
#     paste0("Output/Rdata/scRNA/01_", names(tnbc.xPT.list)[i], "_DEGxJAC_ESSENTIALS_2023.09.21.csv")
#   )
#   
# }
# 
# Idents(tnbc.xPT.list[[8]]) <- tnbc.xPT.list[[8]]@meta.data$jaccard.cluster
# 
# markers <- deg %>%
#   group_by(cluster) %>%
#   top_n(n = 5, wt = avg_log2FC)
# 
# # # doesn't work?
# # avg <- AverageExpression(tnbc.xPT.list[[8]], return.seurat = TRUE)
# 
# DoHeatmap(
#   tnbc.xPT.list[[8]], 
#   assay = "RNA", 
#   group.by = "jaccard.cluster", 
#   features = markers$gene, 
# )
# 
# EnhancedVolcano(
#   deg, 
#   lab = rownames(deg), 
#   x = "avg_log2FC", 
#   y = "p_val_adj", 
#   pCutoff = 5e-2, 
#   FCcutoff = 0.25
# )


################################################################################
#### PERFORM DIFFERENTIAL EXPRESSION BETWEEN CLUSTERS, ALL GENES
for (i in 1:length(gbm.xPT.list)) {

  Idents(gbm.xPT.list[[i]]) <- gbm.xPT.list[[i]]@meta.data$jaccard.cluster

  DEG <- FindAllMarkers(
    gbm.xPT.list[[i]],
    assay = "RNA",
    logfc.threshold = 0.25,
    test.use = "MAST",
    only.pos = TRUE
  )

  write.csv(
    DEG,
    paste0("Output/Rdata/scRNA/01_", names(gbm.xPT.list)[i], "_DEGxJAC_ALL_2023.09.21.csv")
  )

}

################################################################################
# MASSAGING DEG DATA FOR ALL TUMORS

deg.1 <- read.csv(
  "Output/Rdata/scRNA/01_GBM21_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)
deg.2 <- read.csv(
  "Output/Rdata/scRNA/01_GBM41_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)
deg.3 <- read.csv(
  "Output/Rdata/scRNA/01_GBM47_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)
deg.4 <- read.csv(
  "Output/Rdata/scRNA/01_GBM49_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)
deg.5 <- read.csv(
  "Output/Rdata/scRNA/01_GBM51_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)
deg.6 <- read.csv(
  "Output/Rdata/scRNA/01_GBM53_DEGxJAC_ALL_2023.09.21.csv", 
  row.names = 1
)

deg.master <- rbind(
  deg.1, 
  deg.2, 
  deg.3, 
  deg.4, 
  deg.5, 
  deg.6
)

# 37 clusters to start
length(unique(deg.master$cluster)) # 33 clusters (4 don't have DEG)

# filter deg for only essential genes
deg.ess <- deg.master[which(deg.master$gene %in% ess$Gene), ]
deg.ess.filt <- deg.ess[which(deg.ess$p_val < 0.05), ]
# 123 genes remaining / 162
length(unique(deg.ess.filt$cluster)) # down to 22 clusters

################################################################################
# FINAL MASSAGE AND JACCARD SCORING
set.seed(707)

df <- deg.ess.filt[, -c(1:5)]

pivot_df <- df %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = cluster, values_from = present, values_fill = 0) %>%
  as.data.frame() 
rownames(pivot_df) <- pivot_df$gene
pivot_df$gene <- NULL
for (col in 1:length(colnames(pivot_df))) {
  pivot_df[, col] <- as.integer(pivot_df[, col])
}

dir.create("Output/Rdata/Analysis")
write.csv(
  pivot_df, 
  "Output/Rdata/Analysis/00_Jaccard.Input_2_2023.09.25.csv"
)

mtx2 <- simil(t(pivot_df), method = "jaccard") %>% 
  as.matrix()

range(mtx2, na.rm = TRUE)
# [1] 0.0000000 1
mean(mtx2, na.rm = TRUE)
# [1] 0.08467153

mtx2[is.na(mtx2)] <- 0 # should be 1, but to preserve scale on heatmap



write.csv(
  mtx2, 
  "Output/Rdata/Analysis/00_Jaccard.Matrix_2_2023.09.25.csv"
)


generate.colors.1 <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::turbo(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}
generate.colors <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::mako(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

dir.create("Output/Figures/Analysis")
pdf("Output/Figures/Analysis/00_Elbow_Essential.Clusters_transcriptome_2023.09.22.pdf")
fviz_nbclust(mtx2, kmeans, method = "wss", k.max = 12) + 
  theme_minimal() + 
  ggtitle("Elbow Method")
dev.off()

# gap statistic calculation takes too long when too many samples
gap_stat <- clusGap(mtx2, FUN = kmeans, nstart = 30, K.max = 12, B = 50)
gap <- fviz_gap_stat(gap_stat) +
  theme_minimal() +
  ggtitle("Gap Statistic Method")
pdf("Output/Figures/Analysis/00_Gap.Stat_Clusters_transcriptome_2023.09.25.pdf")
gap
dev.off()

clust.rows <- mtx2 %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method, maybe try average linkage? 

ann.col <- data.frame(Cluster = cutree(clust.rows, k = 2))
ann.col$Cluster <- as.character(ann.col$Cluster)

# ann.row <- meta.data[which(rownames(meta.data) %in% rownames(mtx2)), ]
ann <- data.frame(
  clust.ID = rownames(mtx2)
)

split <- function(x) {
  ugh <- strsplit(x, "_")
  return(ugh[[1]][2])
}

ann$tumor.ID <- sapply(ann$clust.ID, split)
rownames(ann) <- ann$clust.ID
ann$clust.ID <- NULL

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]], 
  # state = generate.colors(ann.row)[[2]], 
  # phase = generate.colors(ann.row)[[3]], 
  Cluster = generate.colors.1(ann.col)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/Analysis/01_Essential.Clusters_transcriptome_Heatmap_2023.09.25.pdf"), height = 7, width = 8.5)
pheatmap(
  mtx2,
  color = viridis_pal(direction = -1, option = "magma")(100), #cm.colors(100)
  annotation_row = ann,
  annotation_col = ann.col, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows, 
  cluster_cols = clust.rows, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE
)
dev.off()

rownames(mtx2)
# [1] "1_CID4465"   "3_CID4465"   "4_CID4465"   "5_CID4465"   "1_CID4495"   "2_CID4495"   "3_CID4495"   "4_CID4495"   "5_CID4495"  
# [10] "6_CID4495"   "7_CID4495"   "1_CID44971"  "2_CID44971"  "3_CID44971"  "4_CID44971"  "6_CID44971"  "1_CID44991"  "2_CID44991" 
# [19] "4_CID44991"  "6_CID44991"  "7_CID44991"  "8_CID44991"  "9_CID44991"  "10_CID44991" "11_CID44991" "12_CID44991" "13_CID44991"
# [28] "1_CID4513"   "2_CID4513"   "3_CID4513"   "4_CID4513"   "5_CID4513"   "6_CID4513"   "7_CID4513"   "8_CID4513"   "1_CID4515"  
# [37] "2_CID4515"   "3_CID4515"   "4_CID4515"   "6_CID4515"   "7_CID4515"   "8_CID4515"   "9_CID4515"   "1_CID4523"   "2_CID4523"  
# [46] "3_CID4523"   "4_CID4523"   "5_CID4523"   "6_CID4523"   "7_CID4523"   "8_CID4523"   "9_CID4523"   "10_CID4523"  "11_CID4523" 
# [55] "12_CID4523"  "1_CID3963"   "2_CID3963"   "3_CID3963"   "4_CID3963"   "6_CID3963"   "7_CID3963"  

# need to assign new cluster IDs (1-7) to cells with specific cluster labels
glimpse(tnbc.np.obj@meta.data)
glimpse(tnbc.xPT.list[[1]]@meta.data)

tnbc.np.obj <- AddMetaData(
  tnbc.np.obj,
  clust.data$cluster,
  col.name = "jaccard.cluster"
)

Idents(tnbc.np.obj) <- tnbc.np.obj@meta.data$jaccard.cluster
xClust <- subset(tnbc.np.obj, idents = rownames(mtx2))

scANN <- as.data.frame(cutree(clust.rows, k = 7))
colnames(scANN) <- "cluster"

# cannot use AddMetaData because rownames(scANN) are not cellIDs
xClust@meta.data <- xClust@meta.data %>%
  mutate(
    ess_cluster = 
      case_when(
        jaccard.cluster %in% rownames(scANN) ~ scANN$cluster[match(jaccard.cluster, rownames(scANN))]
      )
  )


saveRDS(
  tnbc.np.obj, 
  "/data/mrd/thesis/scTNBC_NP_DepMap.Clust.Ann_2023.09.21.rds"
)
saveRDS(
  xClust, 
  "Output/Rdata/scRNA/scTNBC_DEG.Clust.Only_2023.09.25.rds"
)

################################################################################
# FINAL DIFFERENTIAL EXPRESSION

# dir.create("Output/Rdata/Analysis")

Idents(xClust) <- xClust$ess_cluster
fin.deg <- FindAllMarkers(
  xClust,
  assay = "RNA", 
  logfc.threshold = 0.25, 
  test.use = "MAST", 
  only.pos = TRUE
)

write.csv(
  fin.deg, 
  "Output/Rdata/Analysis/01_Essential.Clusters_transcriptome_DEGs_2023.09.22.csv"
)

################################################################################
# ENRICHMENT ANALYSIS

library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(reactome.db)

my.theme <- theme(
  panel.background = element_rect(fill = NA), 
  # legend.position = "none", 
  legend.text = element_text(size = 14, colour = "black"), 
  # title = element_blank(), 
  plot.title = element_blank(), 
  axis.line = element_line(linewidth = 1, colour = "black"), 
  axis.text = element_text(size = 14, colour = "black"), 
  axis.title = element_text(size = 14, colour = "black")
)

fin.deg.filt <- fin.deg[which(fin.deg$p_val_adj < 0.05), ]
# fin.deg.filt.2 <- fin.deg.filt %>%
#   group_by(cluster) %>%
#   # Use slice_max to select the top 200 genes within each cluster
#   slice_max(order_by = avg_log2FC, n = 250) %>%
#   ungroup()

f1 <- fin.deg.filt[which(fin.deg.filt$cluster == "1"), ] %>% 
  arrange(desc(avg_log2FC)) 
f2 <- fin.deg.filt[which(fin.deg.filt$cluster == "2"), ] %>% 
  arrange(desc(avg_log2FC)) 
f3 <- fin.deg.filt[which(fin.deg.filt$cluster == "3"), ] %>% 
  arrange(desc(avg_log2FC)) 
f4 <- fin.deg.filt[which(fin.deg.filt$cluster == "4"), ] %>% 
  arrange(desc(avg_log2FC)) 
f5 <- fin.deg.filt[which(fin.deg.filt$cluster == "5"), ] %>% 
  arrange(desc(avg_log2FC)) 
f6 <- fin.deg.filt[which(fin.deg.filt$cluster == "6"), ] %>% 
  arrange(desc(avg_log2FC)) 
f7 <- fin.deg.filt[which(fin.deg.filt$cluster == "7"), ] %>% 
  arrange(desc(avg_log2FC)) 

clust.list <- list(
  mVC1 = f1$gene, 
  mVC2 = f2$gene, 
  mVC3 = f3$gene, 
  mVC4 = f4$gene, 
  mVC5 = f5$gene,
  mVC6 = f6$gene,
  mVC7 = f7$gene
)

msigdb <- msigdbr(species = "Homo sapiens", category = "H")

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list, 
  fun = "enricher", 
  TERM2GENE = msigdb_ref, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  qvalueCutoff = 0.05
)

pdf("Output/Figures/Poster_2023.09.29/06_Essential.Clusters_transcriptome_ClusterProfiler_2023.09.22.pdf", height = 10.5, width = 12)
dotplot(
  comp, 
  x = "Cluster", 
  color = "p.adjust", 
  showCategory = 5, 
  label_format = 50
) + 
  theme_classic() + 
  my.theme
dev.off()


library(dittoSeq)

pdf("Output/Figures/Analysis/03_Essential.Clusters_Barplot_GM_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust, 
  var = "gene.module.ann", 
  group.by = "ess_cluster", 
  scale = "percent", 
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
  xlab = NULL, 
  ylab = "Percent of Cells \n", 
  main = NULL, 
  color.panel = viridis_pal(option = "turbo")(7)
)
dev.off()

pdf("Output/Figures/Analysis/04_Essential.Clusters_Barplot_Tumor_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust, 
  var = "orig.ident", 
  group.by = "ess_cluster", 
  scale = "percent", 
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
  xlab = NULL, 
  ylab = "Percent of Cells \n", 
  main = NULL, 
  color.panel = viridis_pal(option = "turbo")(8)
)
dev.off()

pdf("Output/Figures/Analysis/05_Essential.Clusters_Barplot_Phase_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust, 
  var = "cell.cycle.phase", 
  group.by = "ess_cluster", 
  scale = "percent", 
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
  xlab = NULL, 
  ylab = "Percent of Cells \n", 
  main = NULL, 
  color.panel = viridis_pal(option = "turbo")(8)
)
dev.off()




essentials.STRAT <- fin.deg.filt[which(fin.deg.filt$gene %in% ess$Gene), ] 
head(essentials.STRAT)
# p_val avg_log2FC pct.1 pct.2     p_val_adj cluster   gene
# RAC1   2.614881e-280  0.3353545 0.688 0.894 6.977549e-276       1   RAC1
# PPDPF  1.820259e-273  0.5664593 0.860 0.948 4.857179e-269       1  PPDPF
# COX8A  8.377955e-271  0.2890868 0.633 0.878 2.235574e-266       1  COX8A
# TCEB2  4.188058e-270  0.3470466 0.678 0.891 1.117541e-265       1  TCEB2
# EIF4G2 1.661550e-265  0.2611260 0.407 0.738 4.433681e-261       1 EIF4G2
# KTN1   3.552924e-263  0.3038206 0.455 0.738 9.480621e-259       1   KTN1


# identify cell state-defining genes

# how many genes per state were DEG? 
summary(fin.deg.filt$cluster)
# 1    2    3    4    5    6    7 
# 411  628 1353  461  622 1312 1062

# # extract scaled data
# xclust.scaled <- t(xClust@assays$RNA@scale.data)
# # transpose so cells = columns, genes = rows
# xclust.scaled.t <- as.data.frame(t(xclust.scaled))

Idents(xClust) <- xClust$ess_cluster
DefaultAssay(xClust) <- "RNA"

# take top 200 genes for every cluster
def.genes <- fin.deg.filt %>%
  group_by(cluster) %>%
  # Use slice_max to select the top 200 genes within each cluster
  slice_max(order_by = avg_log2FC, n = 250) %>%
  ungroup()

write.csv(
  def.genes, 
  "Output/Rdata/Analysis/03_top.200.cluster_2023.09.26.csv"
)

# avg.clust <- AverageExpression(xClust, group.by = "ident", slot = "data", return.seurat = TRUE)

# pdf("Output/Figures/Analysis/06_Cluster.Markers_Heatmap_2023.09.26.pdf", height = 7, width = 7)
# DoHeatmap(
#   avg.clust$RNA,
#   features = def.genes$gene,
#   # group.by = "ess_cluster",
#   label = TRUE, 
#   draw.lines = FALSE
# ) +
#   scale_fill_distiller(palette = "RdBu") + #  scale_fill_viridis(option = "mako") +
#   theme(
#     panel.background = element_rect(fill = NA),
#     legend.position = "none",
#     legend.text = element_text(size = 14, colour = "black"),
#     # title = element_blank(),
#     plot.title = element_blank(),
#     # axis.line = element_line(linewidth = 1, colour = "black"),
#     axis.text = element_blank(),
#     axis.title = element_text(size = 14, colour = "black")
#   )
# dev.off()

# def.scaled <- xclust.scaled.t[which(rownames(xclust.scaled.t) %in% def.genes$gene), ]

# def.meta <- as.data.frame(
#   cbind(
#     as.character(xClust$orig.ident),
#     as.character(xClust$jaccard.cluster),
#     as.character(xClust$ess_cluster)
#   )
# )
# colnames(def.meta) <- c("tumor", "jaccard.cluster", "VC")
# rownames(def.meta) <- rownames(xClust@meta.data)
# 
# # sort xpression data by vulnerability cluster
# def.meta <- def.meta %>%
#   arrange(VC)

# order_vec <- match(rownames(def.meta), colnames(def.scaled))
# sorted_df <- def.scaled[, order_vec]
# 
# 
# ann.colors <- list(
#   tumor = generate.colors(def.meta)[[1]], 
#   jaccard.clsuter = generate.colors.1(def.meta)[[2]], 
#   VC = generate.colors.2(def.meta)[[3]]
# )
#

avg.clust <- AverageExpression(xClust, group.by = "ident", slot = "scale.data")

order_vec <- match(unique(def.genes$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% def.genes$gene), ]
# sort expression by gene

gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% ess$Gene)]
)

library(ComplexHeatmap)
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

an.col <- data.frame(
  row.names = c(1:7), 
  mVC = c("mVC1", "mVC2", "mVC3", "mVC4", "mVC5", "mVC6", "mVC7")
)
an.colors <- list(mVC = generate.colors.1(an.col)[[1]])

ph <- pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col, 
  annotation_colors = an.colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  use_raster = TRUE
)
dev.off() 


pdf("Output/Figures/Analysis/06_Cluster.Markers_Heatmap_2023.09.26.pdf", height = 7, width = 7)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()



