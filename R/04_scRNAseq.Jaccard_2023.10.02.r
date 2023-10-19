library(tidyverse)
library(Seurat)
library(viridis)
library(pheatmap)
library(proxy)
library(corrplot)
library(factoextra)
library(FactoMineR)
library(NbClust)
library(cowplot)
library(cluster)
library(MAST)

################################################################################
# GENERATE COLORS FOR HEATMAPS

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

################################################################################
# CREATE RDS OBJECT FOR INDIVIDUAL TUMORS

obj <- readRDS("data/Suter_scGBM_2023.10.04.RDS")

# isolate neoplastic cells only
Idents(obj) <- obj$NPvNon
obj.np <- subset(obj, idents = "Neoplastic")
tum <- as.character(unique(obj.np$orig.ident))

# make obj list that contains np cells for each tumor
Idents(obj.np) <- obj.np$orig.ident
gbm.xPT.list <- list()
for (i in 1:length(tum)) {
  np <- subset(obj.np, idents = tum[i])
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
    col.name = "VC"
  )
}

# get essential genes that had scaled scRNA-seq data
jaccard.input <- read.csv(
  "Output/Rdata/jaccard/00_binary.expr.mtx_2023.09.12.csv",
  row.names = 1
)
ess.genes <- rownames(jaccard.input)
ess <- read.csv("Output/Rdata/04_gbm.killing_2023.09.07.csv")

################################################################################
# PERFORM DIFFERENTIAL EXPRESSION BETWEEN CLUSTERS, ALL GENES

for (i in 1:length(gbm.xPT.list)) {

  Idents(gbm.xPT.list[[i]]) <- gbm.xPT.list[[i]]@meta.data$VC

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
# number of clusters remaining afterr DEG: (good bc no clusters dropped from
# low cell counts)
length(unique(deg.master$cluster)) # 75 clusters
# how many unique genes?
length(unique(deg.master$gene)) # 9233

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

write.csv(
  pivot_df,
  "Output/Rdata/jaccard/04_jaccard.input_2_2023.09.25.csv"
)

#### CALCULATE JACCARD MATRIX --------------------------------------------------

mtx2 <- simil(t(pivot_df), method = "jaccard") %>%
  as.matrix()

range(mtx2, na.rm = TRUE)
# [1] 0.0000000 1
mean(mtx2, na.rm = TRUE)
# [1] 0.08467153

mtx2[is.na(mtx2)] <- 0 # should be 1, but to preserve scale on heatmap

write.csv(
  mtx2,
  "Output/Rdata/jaccard/00_Jaccard.Matrix_2_2023.09.25.csv"
)

################################################################################
# CLUSTERING AND PLOTTING

#### IDENTIFY APPROPRIATE K ----------------------------------------------------

pdf("Output/Figures/00_Elbow_Essential.Clusters_2023.09.22.pdf")
fviz_nbclust(mtx2, kmeans, method = "wss", k.max = 12) +
  theme_minimal() +
  ggtitle("Elbow Method")
dev.off()

# gap statistic calculation takes too long when too many samples
gap_stat <- clusGap(mtx2, FUN = kmeans, nstart = 30, K.max = 12, B = 50)
gap <- fviz_gap_stat(gap_stat) +
  theme_minimal() +
  ggtitle("Gap Statistic Method")
pdf("Output/Figures/00_Gap.Stat_Clusters_2023.09.25.pdf")
gap
dev.off()

#### CLUSTERING AND HEATMAP ----------------------------------------------------

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
pdf(paste0("Output/Figures/Analysis/01_Essential.Clusters_Heatmap_2023.09.25.pdf"), height = 7, width = 8.5)
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

length(rownames(mtx2))


################################################################################
# SINGLE CELL ANNOTATION WITH VC/MVC

# need to assign new cluster IDs (1-7) to cells with specific cluster labels
glimpse(obj.np@meta.data)
glimpse(gbm.xPT.list[[1]]@meta.data)

# add VC cluster data to obj containing np cells only
obj.np <- AddMetaData(
  obj.np,
  clust.data$cluster,
  col.name = "VC"
)

# which clusters are missing?
unique(obj.np$VC)[which(!unique(obj.np$VC) %in% rownames(mtx2))]
# "2_CID44971"  "5_CID44971"  "3_CID44991"  "11_CID44991" "12_CID44991" "17_CID44991" "8_CID4515"
missing.vc <- c("2_CID44971", "5_CID44971", "3_CID44991", "11_CID44991", "12_CID44991", "17_CID44991", "8_CID4515")
for (i in 1:length(missing.vc)) {
  print(table(obj.np$VC[which(obj.np$VC == missing.vc[i])]))
}

# get only single cells that have mVC data (bc of 7 dropped clusters)
Idents(obj.np) <- obj.np@meta.data$VC
xClust <- subset(obj.np, idents = rownames(mtx2))

# add mVC metadata to Seurat object
scANN <- as.data.frame(cutree(clust.rows, k = 2))
colnames(scANN) <- "cluster"

# cannot use AddMetaData because rownames(scANN) are not cellIDs
xClust@meta.data <- xClust@meta.data %>%
  mutate(
    mVC.pre =
      case_when(
        VC %in% rownames(scANN) ~ scANN$cluster[match(VC, rownames(scANN))]
      )
  )

xClust@meta.data <- xClust@meta.data %>%
  mutate(
    mVC = case_when(
      xClust$mVC.pre == 1 ~ "mVC1",
      xClust$mVC.pre == 2 ~ "mVC2",
      xClust$mVC.pre == 3 ~ "mVC3",
      xClust$mVC.pre == 4 ~ "mVC4",
      xClust$mVC.pre == 5 ~ "mVC5",
      xClust$mVC.pre == 6 ~ "mVC6",
      xClust$mVC.pre == 7 ~ "mVC7",
    )
  )
xClust$mVC.pre <- NULL

saveRDS(
  xClust,
  "Output/Rdata/scRNA/scGBM_xClust.Only_2023.09.25.rds"
)


#### ASSIGN METADATA TO OBJ AND OBJ.NP AND SAVE BOTH ---------------------------

# assign VC cluster-level data to obj
obj <- AddMetaData(
  obj,
  obj.np$VC,
  col.name = "VC"
)
obj$VC[is.na(obj$VC)] <- "Non-Neoplastic"
table(obj$VC)

# assign mVC cluster-level data to obj
obj <- AddMetaData(
  obj,
  xClust$mVC,
  col.name = "pre.mVC"
)
dropped.cells <- rownames(obj.np@meta.data)[which(!rownames(obj.np@meta.data) %in% rownames(xClust@meta.data))]
# 1834 dropped cells
obj@meta.data <- obj@meta.data %>%
  mutate(
    mVC = case_when(
      rownames(obj@meta.data) %in% dropped.cells ~ "uncat",
      is.na(obj$pre.mVC) ~ "Non-Neoplastic",
      obj$pre.mVC == 1 ~ "mVC1",
      obj$pre.mVC == 2 ~ "mVC2",
      obj$pre.mVC == 3 ~ "mVC3",
      obj$pre.mVC == 4 ~ "mVC4",
      obj$pre.mVC == 5 ~ "mVC5",
      obj$pre.mVC == 6 ~ "mVC6",
      obj$pre.mVC == 7 ~ "mVC7",
    )
  )
obj$pre.mVC <- NULL
table(obj$mVC)

# save RDS
saveRDS(
  obj,
  "Output/Rdata/scRNA/scGBM_mVC.Ann_2023.09.28.rds"
)

################################################################################
# FINAL DIFFERENTIAL EXPRESSION

Idents(xClust) <- xClust$mVC
fin.deg <- FindAllMarkers(
  xClust,
  assay = "RNA",
  logfc.threshold = 0.25,
  test.use = "MAST",
  only.pos = TRUE
)

write.csv(
  fin.deg,
  "Output/Rdata/05_mVC_DEGs_2023.10.17.csv"
)

# fin.deg <- read.csv(
#   "Output/Rdata/05_mVC_DEGs_2023.10.17.csv",
#   row.names = 1
# )
table(fin.deg$cluster)


