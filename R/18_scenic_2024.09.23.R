library(tidyverse)
library(Seurat)
library(pheatmap)
library(RColorBrewer)

dir.create("Output/Rdata/18_scenic")
dir.create("Output/Figures/18_scenic")

# # ################################################################################
# # # EXPORT DATA AS RAW COUNTS FOR pySCENIC
#
# suter <- readRDS("data/Suter_scGBM_2023.10.04.RDS")
#
# scenic.input <- suter@assays$RNA@counts
# scenic.input <- as.matrix(t(scenic.input))
# # need to confirm columns are genes, rows are cells
# head(scenic.input)[1:5, 1:5]
#
# meta <- suter@meta.data
# write.table(
#   meta,
#   "/data/mrd/cognitive.seeds_zero/Output/Rdata/07_suter_2024.05.06/suter_meta.data.csv",
#   sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE
# )
#
#
# # avg.clust <- AggregateExpression(
# #   suter,
# #   assays = "RNA",
# #   features = rownames(suter[["RNA"]]$counts),
# #   group.by = "mVC",
# #   normalization.method = "CLR",
# #   return.seurat = TRUE
# # )
# #
# # pseudo.counts <- avg.clust[["RNA"]]$counts
# # pseudo.counts <- as.matrix(pseudo.counts)
# # # genes as rows, cells as columns
# # write.table(
# #   pseudo.counts,
# #   "Output/Rdata/01_pseudobulk.np.cell.counts.tsv",
# #   sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE
# # )
#
# # genes as columns, cells as rows
# write.table(
#   scenic.input,
#   "/data/mrd/scenic/gbm/expr_mat.tsv",
#   sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE
# )
#
#
# # dir.create("Output/Figures/scenic")

################################################################################
# READ IN DATA

obj <- readRDS("data/Suter_scGBM.NP_2023.10.04.RDS")

auc <- read.csv("Output/Rdata/18_scenic/auc_mtx.csv", header = TRUE, row.names = 1)
glimpse(auc)
head(auc)[1:5]

################################################################################
# DATA CLEAN

# clean up auc column names

cols <- as.data.frame(colnames(auc))

for (i in 1:length(rownames(cols))) {
  cols$newCol[i] <- unlist(strsplit(
    cols$`colnames(auc)`[i],
    split = "...",
    fixed = TRUE
  ))[1]
}

colnames(auc) <- cols$newCol

################################################################################
# IDK

# pdf("Output/Figures/scenic/scTNBC_DimPlot.pdf", height = 7, width = 8.5)
# DimPlot(
#   obj,
#   cols.highlight = viridis_pal(option = "turbo", direction = -1)(6)
# ) +
#   theme(
#     panel.background = element_rect(fill = NA),
#     legend.position = "right",
#     plot.title = element_blank(),
#     axis.line = element_blank(),
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank()
#   )
# dev.off()

################################################################################
# AUC DATA ANNOTATION AND SUBSET

ann <- data.frame(
  Tumor = obj$orig.ident,
  Cell.Type = obj$CellType,
  Seurat.Clust = obj$seurat_clusters,
  CS = obj$ClassificationNPvsNon,
  mVC = obj$mVC
)

auc.ann <- auc %>%
  mutate(
    mVC = case_when(
      rownames(auc) %in% rownames(ann) ~ paste0(ann$mVC)
    )
  )

sub.auc.ann <- subset(auc.ann, auc.ann$mVC != "Non-Neoplastic")

aucXmvc <- sub.auc.ann %>%
  group_by(mVC) %>%
  summarise_all(.funs = mean)
aucXmvc <- as.data.frame(aucXmvc)
rownames(aucXmvc) <- aucXmvc$mVC
aucXmvc$mVC <- NULL

################################################################################


annot <- data.frame(
  row.names = rownames(auc.ann),
  mVC = auc.ann$mVC
)

my.breaks <- c(seq(min(auc), 0, length.out=ceiling(100/2) + 1),
               seq(max(auc)/100, max(auc), length.out=floor(100/2)))

pdf("Output/Figures/18_scenic/00_Suter_RegulonActivity_2024.11.05.pdf", height = 8, width = 8)
ComplexHeatmap::pheatmap(
  t(as.matrix(auc)),
  annotation_col = annot,
  scale = "row",
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  show_rownames = FALSE,
  show_colnames = FALSE,
  # breaks = my.breaks,
  use_raster = TRUE
)
dev.off()

################################################################################
# SCALE AND PLOT

new.auc <- as.data.frame(t(aucXmvc))
pdf("Output/Figures/18_scenic/01_Suter_RegulonActivity.mVC_2024.11.06.pdf", height = 8, width = 8)
pheatmap(
  new.auc,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  show_rownames = FALSE,
  scale = "row"
)
dev.off()

# scale_aucXmvc <- as.data.frame(t(scale(aucXmvc, center = TRUE, scale = TRUE)))
#
# my.breaks <- c(seq(min(scale_aucXmvc), 0, length.out=ceiling(100/2) + 1),
#                seq(max(scale_aucXmvc)/100, max(scale_aucXmvc), length.out=floor(100/2)))
#
# pdf("Output/Figures/18_scenic/01_Suter_RegulonActivity.mVC_2024.01.31.pdf", height = 7, width = 7)
# pheatmap(
#   scale_aucXmvc,
#   color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
#   show_rownames = FALSE,
#   breaks = my.breaks
# )
# dev.off()


################################################################################
# IDENTIFY TOP REGULONS

ann.top <- new.auc
ann.top$regulons <- rownames(ann.top)

top.regulons <- reshape2::melt(ann.top)
glimpse(top.regulons)
colnames(top.regulons) <- c("regulon", "mVC", "relative.activity")

# TOP 25 REGULONS PER CELL STATE
top.25 <- top.regulons %>%
  group_by(mVC) %>%
  slice_max(n = 25, order_by = relative.activity)
glimpse(top.25)

scale.25 <- subset(new.auc, rownames(new.auc) %in% unique(top.25$regulon))

pdf("Output/Figures/18_scenic/02_Suter_RegulonActivity.mVC_top25_2024.11.06.pdf", height = 8, width = 8)
pheatmap(
  scale.25,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45",
  scale = "row"
)
dev.off()

# TOP 10 REGULONS PER CELL STATE
top.10 <- top.regulons %>%
  group_by(mVC) %>%
  slice_max(n = 10, order_by = relative.activity)
glimpse(top.10)

scale.10 <- subset(new.auc, rownames(new.auc) %in% unique(top.10$regulon))

pdf("Output/Figures/18_scenic/03_Suter_RegulonActivity.mVC_top10_2024.11.06.pdf", height = 8, width = 8)
pheatmap(
  scale.10,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45",
  scale = "row"
)
dev.off()

################################################################################
# if you scale the matrix first, then look for the top regulons,
# you get more distinct results
scale_aucXmvc <- as.data.frame(t(scale(aucXmvc, center = TRUE, scale = TRUE)))

ann.top <- scale_aucXmvc
ann.top$regulons <- rownames(ann.top)

top.regulons <- reshape2::melt(ann.top)
glimpse(top.regulons)
colnames(top.regulons) <- c("regulon", "mVC", "relative.activity")

# TOP 25 REGULONS PER CELL STATE
top.25 <- top.regulons %>%
  group_by(mVC) %>%
  slice_max(n = 25, order_by = relative.activity)
glimpse(top.25)

scale.10 <- subset(scale_aucXmvc, rownames(scale_aucXmvc) %in% unique(top.25$regulon))

my.breaks.1 <- c(seq(min(scale.10), 0, length.out=ceiling(100/2) + 1),
                 seq(max(scale.10)/100, max(scale.10), length.out=floor(100/2)))

pdf("Output/Figures/18_scenic/02_Suter_RegulonActivity.mVC_top25_2024.01.31.pdf", height = 14, width = 14)
pheatmap(
  scale.10,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45", fontsize_row = 8, fontsize_col = 15,
  breaks = my.breaks.1
)
dev.off()

# TOP 10 REGULONS PER CELL STATE
top.10 <- top.regulons %>%
  group_by(mVC) %>%
  slice_max(n = 10, order_by = relative.activity)
glimpse(top.10)

scale.20 <- subset(scale_aucXmvc, rownames(scale_aucXmvc) %in% unique(top.10$regulon))

my.breaks.1 <- c(seq(min(scale.20), 0, length.out=ceiling(100/2) + 1),
                 seq(max(scale.20)/100, max(scale.20), length.out=floor(100/2)))

pdf("Output/Figures/18_scenic/03_Suter_RegulonActivity.mVC_top10_2024.01.31.pdf", height = 7, width = 7)
pheatmap(
  scale.20,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45", fontsize_row = 8, fontsize_col = 15,
  breaks = my.breaks.1
)
dev.off()

################################################################################
# LIMMA FOR SIGNIFICANTLY DIFFERENT REGULONS

# organize data matrix
mtx <- sub.auc.ann[,-c(length(colnames(sub.auc.ann)))]
dat.mtx <- as.matrix(t(mtx))

# create factor for groups
group <- factor(sub.auc.ann$mVC)

# create design matrix
design <- model.matrix(~0 + group)  # ~0 removes the intercept, so each group is a separate column
colnames(design) <- levels(group)

fit <- lmFit(dat.mtx, design)

# Define contrasts to compare the groups, for example:
contrast_matrix <- makeContrasts(
  VS1_vs_VS2 = mVC1 - mVC2,
  VS1_vs_VS3 = mVC1 - mVC3,
  VS2_vs_VS3 = mVC2 - mVC3,
  levels = design
)

# Apply the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes moderation
fit2 <- eBayes(fit2)
saveRDS(fit2, "Output/Rdata/18_scenic/limma_results_2024.11.13.RDS")

# Get the top table of results for each contrast
results_VS1.VS2 <- topTable(fit2, coef = "VS1_vs_VS2", number = Inf, adjust.method = "bonferroni", p.value = 0.05)
results_VS1.VS3 <- topTable(fit2, coef = "VS1_vs_VS3", number = Inf, adjust.method = "bonferroni", p.value = 0.05)
results_VS2.VS3 <- topTable(fit2, coef = "VS2_vs_VS3", number = Inf, adjust.method = "bonferroni", p.value = 0.05)

# Get the top table of results for each contrast
res.list <- list()
for (i in 1:length(colnames(contrast_matrix))) {
  comp <- colnames(contrast_matrix)[i]
  res <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH", p.value = 0.05)
  res.list[[i]] <- res
}

names(res.list) <- colnames(contrast_matrix)
saveRDS(res.list, "Output/Rdata/18_scenic/scenic_diff.regulons_2024.11.13.RDS")

# merge results across comparisons
df_list <- lapply(names(res.list), function(name) {
  df <- res.list[[name]]
  var1 <- strsplit(name, "_")[[1]][1]
  var2 <- strsplit(name, "_")[[1]][3]

  if (length(df) != 0) {
    df <- df %>%
      mutate(
        Gene = rownames(df),
        Comparison = name, # comparison annotation
        State = case_when( # state annotation
          logFC > 0 ~ var1,
          logFC <= 0 ~ var2
        )
      )
  }
  else {
    df <- data.frame()
  }
  return(df)
})
glimpse(df_list)

combined <- bind_rows(df_list)
glimpse(combined)
length(unique(combined$Gene))
write.csv(
  combined,
  "Output/Rdata/18_scenic/scenic_diff.regulons_annotated_2024.11.13.csv",
  row.names = TRUE
)

top <- combined %>%
  group_by(State) %>%
  slice_max(n = 250, order_by = logFC)
clust.list <- split(x = combined.grp$Gene, f = combined.grp$State)
glimpse(clust.list)

plot <- subset(new.auc, rownames(new.auc) %in% top$Gene)

pdf("Output/Figures/18_scenic/05_topLIMMA_250_2024.11.13.pdf", height = 6, width = 12)
pheatmap(
  t(plot),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # angle_col = "45",
  scale = "column", show_colnames = FALSE, # fontsize_col = 6
)
dev.off()


top <- combined %>%
  group_by(State) %>%
  slice_max(n = 25, order_by = logFC)
clust.list <- split(x = combined.grp$Gene, f = combined.grp$State)
glimpse(clust.list)

plot <- subset(new.auc, rownames(new.auc) %in% top$Gene)

pdf("Output/Figures/18_scenic/05_topLIMMA_25_2024.11.13.pdf", height = 4, width = 12)
pheatmap(
  t(plot),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45",
  scale = "column", show_colnames = TRUE, fontsize_col = 6
)
dev.off()


top <- combined %>%
  group_by(State) %>%
  slice_max(n = 50, order_by = logFC)
clust.list <- split(x = combined.grp$Gene, f = combined.grp$State)
glimpse(clust.list)

plot <- subset(new.auc, rownames(new.auc) %in% top$Gene)

pdf("Output/Figures/18_scenic/05_topLIMMA_50_2024.11.13.pdf", height = 4, width = 12)
pheatmap(
  t(plot),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45",
  scale = "column", show_colnames = TRUE, fontsize_col = 4
)
dev.off()






# VS1-specific regulons
sig.1.2 <- results_VS1.VS2[which(results_VS1.VS2$logFC > 0), ]
sig.1.3 <- results_VS1.VS3[which(results_VS1.VS3$logFC > 0), ]
sig.1.2_25 <- sig.1.2 %>%
  slice_max(n = 25, order_by = logFC)
sig.1.3_25 <- sig.1.3 %>%
  slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.1.2_25)), unique(rownames(sig.1.3_25))))

# VS2-specific drugs
sig.2.1 <- results_VS1.VS2[which(results_VS1.VS2$logFC < 0), ]
sig.2.3 <- results_VS2.VS3[which(results_VS2.VS3$logFC > 0), ]
sig.2.1_25 <- sig.2.1 %>%
  slice_min(n = 25, order_by = logFC)
sig.2.3_25 <- sig.2.3 %>%
  slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.2.1_25)), unique(rownames(sig.2.3_25))))

# VS3-specific drugs
sig.3.1 <- results_VS1.VS3[which(results_VS1.VS3$logFC < 0), ]
sig.3.2 <- results_VS2.VS3[which(results_VS2.VS3$logFC < 0), ]
sig.3.1_25 <- sig.3.1 %>%
  slice_min(n = 25, order_by = logFC)
sig.3.2_25 <- sig.3.2 %>%
  slice_min(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.3.1_25)), unique(rownames(sig.3.2_25))))

top.reg <- c(rownames(sig.1.2_25), rownames(sig.1.3_25),
               rownames(sig.2.1_25), rownames(sig.2.3_25),
               rownames(sig.3.1_25), rownames(sig.3.2_25))
top.reg <- unique(top.reg) ####


top <- subset(new.auc, rownames(new.auc) %in% top.reg)

pdf("Output/Figures/18_scenic/04_Suter_RegulonActivity.mVC_topLIMMA_2024.11.11.pdf", height = 6, width = 12)
pheatmap(
  t(top),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  angle_col = "45",
  scale = "column", fontsize_col = 6
)
dev.off()
