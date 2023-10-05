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

# dir.create("Output/Rdata/jaccard")
# dir.create("Output/Figures/jaccard")

# obj <- readRDS("/home/mrd78/DEPMAP/DepMap.Scored_scGBM_2023.02.21.RDS")

ess <- read.csv("Output/Rdata/04_gbm.killing_2023.09.07.csv")
# # ess.ref <- read.csv("Output/Rdata/scRNA/04_REF_gbm.killing_2023.09.07.csv")
# ess.xPT.list <- readRDS("Output/Rdata/scRNA/xPT/00_REF_xPT_list_2023.09.11.rds")

# cells = rows, genes = columns
scaled.expr <- read.csv("Output/Rdata/scRNA/scaled.rna.expr_gbm_suter_2023.09.07.csv")
scaled.expr <- as.data.frame(t(scaled.expr))

################################################################################
# EXTRACT SCALED EXPRESSION OF ESSENTIAL GENES

scaled.ess <- scaled.expr[which(rownames(scaled.expr) %in% ess$Gene), ] # 162/168 ess genes
rm(scaled.expr)
gc()

write.csv(
  scaled.ess, 
    "Output/Rdata/jaccard/00_scaled.np.expr_ess.only_2023.09.12.csv"
)

################################################################################
# CALCULATE BINARY EXPRESSION MATRIX

# generate matrix of binary essential gene expression (1, 0)
# scaled expression > 0

jaccard.input <- scaled.ess

binarize <- function(val) {
  ifelse(val > 0, 1, 0)
}

jaccard.input[] <- lapply(jaccard.input, binarize)

write.csv(
    jaccard.input, 
    "Output/Rdata/jaccard/00_binary.expr.mtx_2023.09.12.csv"
)

################################################################################
# CALCULATE JACCARD SIMILARITY MATRIX

# # jaccard calculation function
# jaccard <- function(x, y) {
#     i.len <- length(intersect(x, y)) 
#     u.len <- length(union(x, y))
#    return(round(i.len / u.len, 3))
# } DOESN'T WORK BC THE DATA IS BINARY

mtx <- simil(t(jaccard.input), method = "jaccard") %>%
  as.matrix()

range(mtx, na.rm = TRUE)
# [1] 0 1

write.csv(
  mtx, 
  "Output/Rdata/jaccard/01_jaccard.matrix_2023.09.12.csv"
)

################################################################################
# VISUALIZATION

# annotation data for heatmap
meta.data <- as.data.frame(
    cbind(
        as.character(gbm.np$orig.ident), 
        as.character(gbm.np$CellStateClassification4),
        as.character(gbm.np$Phase)
    )
)
colnames(meta.data) <- c("tumor", "state", "phase")
rownames(meta.data) <- rownames(gbm.np@meta.data)

write.csv(
  meta.data,
  "Output/Rdata/jaccard/02_jaccard.meta.data_2023.09.15.csv"
)

set.seed(707)

# pdf("Output/Figures/jaccard/01_jaccard.ComplexHeatmap_clustered_2023.09.12.pdf", height = 7, width = 8.2)
# Heatmap(
#   mtx,
#   use_raster = TRUE,
#   cluster_rows = TRUE,
#   cluster_columns = TRUE,
#   show_row_names = FALSE,
#   show_column_names = FALSE
# )
# dev.off()

################################################################################
# EXTRACT JACCARD MATRIX BY TUMOR

split.cells <- strsplit(rownames(mtx), "_")
tumorID <- unique(sapply(split.cells, function(x) x[1])) # isolate tumor ids

tumor.mtx <- list()

for (tumor_id in tumorID) {
  tumor_rows <- grep(tumor_id, rownames(mtx))
  tumor_cols <- grep(tumor_id, colnames(mtx))
  tumor.mtx[[tumor_id]] <- mtx[tumor_rows, tumor_cols]
}
rm(split.cells)

################################################################################
# IDENTIFY OPTIMAL K FOR EACH TUMOR
set.seed(707)
names(tumor.mtx)
# [1] "GBM21" "GBM41" "GBM47" "GBM49" "GBM51" "GBM53"

# determine optimal number of clusters for kmeans clustering for every tumor

for (i in 1:length(tumor.mtx)) {
  
  jac <- tumor.mtx[[i]]
  jac[is.na(jac)] <- 1
  
  # elbow <- fviz_nbclust(jac, kmeans, method = "wss", k.max = 25) + 
  #   theme_minimal() + 
  #   ggtitle("Elbow Method")
  
  # gap statistic calculation takes too long when too many samples
  gap_stat <- clusGap(jac, FUN = kmeans, K.max = 20, B = 50)
  gap <- fviz_gap_stat(gap_stat) +
    theme_minimal() +
    ggtitle("Gap Statistic Method")

  
  # silhouette <- fviz_nbclust(jac, kmeans, method = "silhouette", k.max = 12) +
  #   theme_minimal() +
  #   ggtitle("Silhouette Method")
  
  pdf(paste0("Output/Figures/QC/optimal.clusters_GAP_", names(tumor.mtx)[i], "_2023.09.20.pdf"), height = 7, width = 7)
  print(gap)
  dev.off()
  
}

# # gap statistic calculation takes too long when too many samples
# gap_stat <- clusGap(jac, FUN = kmeans, K.max = 25, B = 50)
# gap <- fviz_gap_stat(gap_stat) +
#   theme_minimal() +
#   ggtitle("Gap Statistic Method")
# pdf("Output/Figures/QC/optimal.clusters_GAP_GBM21_2023.09.20.pdf")
# gap
# dev.off()

generate.colors <- function(df) { # tumor
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::mako(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}
generate.colors.1 <- function(df) { # cluster
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::turbo(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}
generate.colors.2 <- function(df) { # state
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::rocket(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}
generate.colors.3 <- function(df) { # phase
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis::plasma(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

################################################################################
# CLUSTER INDIVIDUAL TUMORS 

#### TUMOR 1 ##############################

# [1] "GBM21" "GBM41" "GBM47" "GBM49" "GBM51" "GBM53"

# compute clustering via Ward's method
clust.rows.1 <- tumor.mtx[[1]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col <- data.frame(Cluster = cutree(clust.rows.1, k = 12))
ann.col$Cluster <- as.character(ann.col$Cluster)

ann.row <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[1]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row)[[1]], 
  state = generate.colors.2(ann.row)[[2]], 
  phase = generate.colors.3(ann.row)[[3]], 
  Cluster = generate.colors.1(ann.col)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[1], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[1]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[1]])), ],
  annotation_col = ann.col, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows.1, 
  cluster_cols = clust.rows.1, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
) 
dev.off()

# generate clustering metadata data frame
clust.data.1 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[1]])), ],
  cluster = paste0(cutree(clust.rows.1, k = 12), "_", names(tumor.mtx)[1])
)

table(clust.data.1$cluster) # #cells/cluster
# 1_GBM21 2_GBM21 3_GBM21 4_GBM21 5_GBM21 6_GBM21 
# 395     137     173      73     436     186 

write.csv(
  clust.data.1, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[1], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[1], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.1 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[1], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.1 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[1], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.1 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[1], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.1 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

#### TUMOR 2 ##############################

# compute clustering via Ward's method
clust.rows.2 <- tumor.mtx[[2]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col.2 <- data.frame(Cluster = cutree(clust.rows.2, k = 6))
ann.col.2$Cluster <- as.character(ann.col.2$Cluster)

ann.row.2 <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[2]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row.2)[[1]], 
  state = generate.colors.2(ann.row.2)[[2]], 
  phase = generate.colors.3(ann.row.2)[[3]], 
  Cluster = generate.colors.1(ann.col.2)[[1]]
)

library(ComplexHeatmap)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[2], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[2]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = ann.row.2,
  annotation_col = ann.col.2, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows.2, 
  cluster_cols = clust.rows.2, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
)
dev.off()

# generate clustering metadata data frame
clust.data.2 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[2]])), ],
  cluster = paste0(cutree(clust.rows.2, k = 6), "_", names(tumor.mtx)[2])
)

table(clust.data.2$cluster) # #cells/cluster
# 1_GBM41 2_GBM41 3_GBM41 4_GBM41 5_GBM41 6_GBM41 
# 2211    1180     762    1115    2177     227 

write.csv(
  clust.data.2, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[2], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[2], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.2 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[2], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.2 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[2], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.2 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[2], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.2 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

#### TUMOR 3 ##############################

# compute clustering via Ward's method
clust.rows.3 <- tumor.mtx[[3]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col.3 <- data.frame(Cluster = cutree(clust.rows.3, k = 6))
ann.col.3$Cluster <- as.character(ann.col.3$Cluster)

ann.row.3 <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[3]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row.3)[[1]], 
  state = generate.colors.2(ann.row.3)[[2]], 
  phase = generate.colors.3(ann.row.3)[[3]], 
  Cluster = generate.colors.1(ann.col.3)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[3], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[3]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = ann.row.3,
  annotation_col = ann.col.3, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows.3, 
  cluster_cols = clust.rows.3, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
)
dev.off()

# generate clustering metadata data frame
clust.data.3 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[3]])), ],
  cluster = paste0(cutree(clust.rows.3, k = 6), "_", names(tumor.mtx)[3])
)

table(clust.data.3$cluster) # #cells/cluster
# 1_GBM47 2_GBM47 3_GBM47 4_GBM47 5_GBM47 6_GBM47 
#   1299     571     984    1313    1493    1138

write.csv(
  clust.data.3, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[3], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[3], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.3 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[3], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.3 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[3], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.3 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[3], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.3 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

#### TUMOR 4 ##############################

# compute clustering via Ward's method
clust.rows.4 <- tumor.mtx[[4]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col.4 <- data.frame(Cluster = cutree(clust.rows.4, k = 6))
ann.col.4$Cluster <- as.character(ann.col.4$Cluster)

ann.row.4 <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[4]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row.4)[[1]], 
  state = generate.colors.2(ann.row.4)[[2]], 
  phase = generate.colors.3(ann.row.4)[[3]], 
  Cluster = generate.colors.1(ann.col.4)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[4], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[4]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = ann.row.4,
  annotation_col = ann.col.4,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows.4, 
  cluster_cols = clust.rows.4, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
)
dev.off()

# generate clustering metadata data frame
clust.data.4 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[4]])), ],
  cluster = paste0(cutree(clust.rows.4, k = 6), "_", names(tumor.mtx)[4])
)

table(clust.data.4$cluster) # #cells/cluster
# 1_GBM49 2_GBM49 3_GBM49 4_GBM49 5_GBM49 6_GBM49 
#   2433     417    1479    2964     349    1130 

write.csv(
  clust.data.4, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[4], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[4], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.4 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[4], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.4 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[4], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.4 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[4], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.4 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

#### TUMOR 5 ##############################

# compute clustering via Ward's method
clust.rows.5 <- tumor.mtx[[5]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col.5 <- data.frame(Cluster = cutree(clust.rows.5, k = 6))
ann.col.5$Cluster <- as.character(ann.col.5$Cluster)

ann.row.5 <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[5]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row.5)[[1]], 
  state = generate.colors.2(ann.row.5)[[2]], 
  phase = generate.colors.3(ann.row.5)[[3]], 
  Cluster = generate.colors.1(ann.col.5)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[5], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[5]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = ann.row.5,
  annotation_col = ann.col.5, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows.5, 
  cluster_cols = clust.rows.5, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
)
dev.off()

# generate clustering metadata data frame
clust.data.5 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[5]])), ],
  cluster = paste0(cutree(clust.rows.5, k = 6), "_", names(tumor.mtx)[5])
)

table(clust.data.5$cluster) 
# 1_GBM51 2_GBM51 3_GBM51 4_GBM51 5_GBM51 6_GBM51 
#   1202     998     956    1154     660     118

write.csv(
  clust.data.5, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[5], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[5], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.5 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[5], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.5 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[5], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.5 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[5], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.5 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

#### TUMOR 6 ##############################

# compute clustering via Ward's method
clust.rows.6 <- tumor.mtx[[6]] %>%
  dist(method = "euclidean") %>% # compute dissimilarity matrix
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

ann.col.6 <- data.frame(Cluster = cutree(clust.rows.6, k = 7))
ann.col.6$Cluster <- as.character(ann.col.6$Cluster)

ann.row.6 <- meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[6]])), ]

ann.colors.col <- list(
  tumor = generate.colors(ann.row.6)[[1]], 
  state = generate.colors.2(ann.row.6)[[2]], 
  phase = generate.colors.3(ann.row.6)[[3]], 
  Cluster = generate.colors.1(ann.col.6)[[1]]
)

# generate heatmap object
pdf(paste0("Output/Figures/jaccard/04_", names(tumor.mtx)[6], "_jaccaard.pheatmap_clustered_2023.09.18.pdf"), height = 7, width = 8.2)
pheatmap(
  tumor.mtx[[6]],
  color = viridis_pal(direction = -1, option = "magma")(100), 
  annotation_row = ann.row.6,
  annotation_col = ann.col.6, 
  annotation_colors = ann.colors.col, 
  cluster_rows = clust.rows.6, 
  cluster_cols = clust.rows.6, 
  show_rownames = FALSE,
  show_colnames = FALSE, 
  use_raster = TRUE, 
)
dev.off()

# generate clustering metadata data frame
clust.data.6 <- cbind(
  meta.data[which(rownames(meta.data) %in% rownames(tumor.mtx[[6]])), ],
  cluster = paste0(cutree(clust.rows.6, k = 7), "_", names(tumor.mtx)[6])
)

table(clust.data.6$cluster)
# 1_GBM53 2_GBM53 3_GBM53 4_GBM53 5_GBM53 6_GBM53 7_GBM53 
# 333     471      87     695     416     310     117 

write.csv(
  clust.data.6, 
  paste0("Output/Rdata/jaccard/03_", names(tumor.mtx)[6], "_cluster.data_2023.09.18.csv")
)

#### SUMMARY STATISTICS - CATEGORICAL

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[6], "_summary.stats_STATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.6 %>%
  dplyr::count(state) %>%
  ggplot(aes(x = factor(state), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[6], "_summary.stats_PHASE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.6 %>%
  dplyr::count(phase) %>%
  ggplot(aes(x = factor(phase), y = n, fill = factor(phase))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[6], "_summary.stats_CLUSTER_2023.09.18.pdf"), height = 7, width = 7)
clust.data.6 %>%
  dplyr::count(cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Count", fill = "") + 
  geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()

pdf(paste0("Output/Figures/jaccard/05_", names(tumor.mtx)[6], "_summary.stats_CLUSTERxSTATE_2023.09.18.pdf"), height = 7, width = 7)
clust.data.6 %>%
  dplyr::count(state, cluster) %>%
  ggplot(aes(x = factor(cluster), y = n, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "Cluster", y = "Count", fill = "State") + 
  # geom_text(stat = "identity", aes(label = n), vjust = -0.5, size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)))
dev.off()
