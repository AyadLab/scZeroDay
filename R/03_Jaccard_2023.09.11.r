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

gbm.seurat <- readRDS("data/Suter_scGBM_2023.10.04.RDS")

ess <- read.csv("Output/Rdata/04_gbm.killing_2023.09.07.csv")
# # ess.ref <- read.csv("Output/Rdata/scRNA/04_REF_gbm.killing_2023.09.07.csv")
# ess.xPT.list <- readRDS("Output/Rdata/scRNA/xPT/00_REF_xPT_list_2023.09.11.rds")

################################################################################
# EXTRACT SCALED EXPRESSION OF ESSENTIAL GENES

Idents(gbm.seurat) <- gbm.seurat$NPvNon
gbm.np <- subset(gbm.seurat, idents = "Neoplastic")

# extract scaled expression data
scaled <- as.data.frame(gbm.np@assays$RNA@scale.data)

# subset for only essential genes
scaled.ess <- scaled[which(rownames(scaled) %in% ess$Gene), ] # 162/168 ess genes
rm(scaled)
rm(gbm.np)
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

saveRDS(
  tumor.mtx,
  "Output/Rdata/jaccard/03_jaccard.mtx_xPT_2023.010.10.RDS"
)
