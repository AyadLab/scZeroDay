library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(viridis)
library(EnhancedVolcano)
library(MAST)
library(SummarizedExperiment)
library(clustree)
library(dittoSeq)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(msigdbr)
library(clusterProfiler)

dir.create("Output/Rdata/02_scRNA.neftel")
dir.create("Output/Figures/02_scRNA.neftel")
dir.create("Output/Figures/02_scRNA.neftel/scUMAP")

mrd.theme <- readRDS("Output/mrd.fig.theme.RDS")

umap.theme <- theme(
  panel.background = element_rect(fill = NA),
  title = element_blank(),
  plot.title = element_blank(),
  legend.position = "none",
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank()
)
saveRDS(umap.theme, "Output/mrd.umap.theme.RDS")

# # color palettes
# # RED/WHITE/BLUE (CONTINUOUS, DIVERGING, HEATMAPS)
# paletteer::scale_colour_paletteer_c("grDevices::RdBu")
# paletteer::scale_color_paletteer_c("grDevices::RdBu")
# paletteer::scale_fill_paletteer_c("grDevices::RdBu")
# paletteer::paletteer_c("grDevices::RdBu")
# # VIRIDIS INFERNO BLACK/RED/YELLOW (CONTINUOUS, SEQUENTIAL, HEATMAPS)
# paletteer::scale_colour_paletteer_c("viridis::inferno")
# paletteer::scale_color_paletteer_c("viridis::inferno")
# paletteer::scale_fill_paletteer_c("viridis::inferno")
# paletteer::paletteer_c("viridis::inferno")
# # VIRIDIS TURBO (DISCRETE, QUALITATIVE)
# viridis::viridis_pal(option = "turbo")
# viridis::scale_color_viridis(option = "turbo")
# viridis::scale_fill_viridis(option = "turbo")
## RCOLOR BREWER RDBU (DISCRETE, SIMILAR COLORS)
# colorRampPalette(brewer.pal(n = 11, name = "RdBu"))


################################################################################
# READ DATA

# Neftel et al. scRNA-seq data (pre-processed)
obj <- readRDS("data/neftel_ss_2023.11.10.RDS")

# remove pediatric samples
Idents(obj) <- obj$GBMType
obj <- subset(obj, idents = "Adult")

# read in list of essential genes
killing <- read.table(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv",
  row.names = 1, header = TRUE, sep = ","
)

################################################################################
# ANNOTATE NEFTEL DATA OBJECT####

glimpse(obj@meta.data)

obj@meta.data <- obj@meta.data %>%
  mutate(
    deconv = case_when(
      CellAssignment == "Macrophage" ~ CellAssignment,
      CellAssignment == "Malignant" ~ cell.state,
      CellAssignment == "Oligodendrocyte" ~ CellAssignment,
      CellAssignment == "T-cell" ~ CellAssignment
    )
  )

obj@meta.data <- obj@meta.data %>%
  mutate(
    CellAssignment = case_when(
      CellAssignment == "Macrophage" ~ CellAssignment,
      CellAssignment == "Malignant" ~ "Neoplastic",
      CellAssignment == "Oligodendrocyte" ~ CellAssignment,
      CellAssignment == "T-cell" ~ CellAssignment
    )
  )

# correct improper tumor ID
obj@meta.data <- obj@meta.data %>%
  mutate(
    tumor.id = case_when(
      Sample == "BT85" ~ "MGH85",
      Sample != "BT85" ~ Sample
    )
  )

################################################################################
# EXPLORE, FILTER, FIX ####

glimpse(obj@meta.data)
table(obj$deconv)
table(obj$CellAssignment)
table(obj$neoplastic.state)
table(obj$cell.state)
table(obj$tumor.id)

# cleaning up and checking cell.state
obj@meta.data <- obj@meta.data %>%
  mutate(
    cell.state = case_when(
      cell.state == "Unassigned" ~ deconv,
      cell.state != "Unassigned" ~ cell.state
    )
  )
table(obj$cell.state)

obj@meta.data <- obj@meta.data %>%
  mutate(
    cell.state = case_when(
      cell.state == "AC" ~ cell.state,
      cell.state == "MES1" ~ cell.state,
      cell.state == "MES2" ~ cell.state,
      cell.state == "OPC" ~ cell.state,
      cell.state == "NPC1" ~ cell.state,
      cell.state == "NPC2" ~ cell.state,
      cell.state == "Unassigned" ~ cell.state,
      cell.state == "Macrophage" ~ "Non-Neoplastic",
      cell.state == "Oligodendrocyte" ~ "Non-Neoplastic",
      cell.state == "T-cell" ~ "Non-Neoplastic"
    )
  )
table(obj$cell.state)

# remove all "Unassigned" neoplastic cells
Idents(obj) <- obj@meta.data$cell.state
obj <- subset(obj, idents = "Unassigned", invert = TRUE)

# # confirm
Idents(obj) <- obj$neoplastic.state
pdf("Output/Figures/02_scRNA.neftel/00_cell.state.proportions_2024.05.01.pdf", height = 4, width = 8)
dittoBarPlot(
  subset(obj, idents = "Neoplastic"),
  var = "cell.state",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(8)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/00_cell.type.proportions_2024.05.01.pdf", height = 4, width = 8)
dittoBarPlot(
  obj,
  var = "CellAssignment",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(8)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

################################################################################
# PROCESS ####

dir.create("Output/Figures/02_scRNA.neftel/scPCA")

# re-scale expression data on only neoplastic cells
DefaultAssay(obj) <- "integrated"
obj <- FindVariableFeatures(obj, selection.method = "vst")
# this will have only been done on variable features, not all genes
obj <- ScaleData(
  obj,
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)

# principal component analysis
obj <- RunPCA(obj, verbose = FALSE, seed.use = 42)

set.seed(123)
pdf("Output/Figures/02_scRNA.neftel/scPCA/00_Elbow_2023.06.07.pdf", height = 8, width = 8)
ElbowPlot(obj, ndims = 50)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/01_PCA1v2_2023.06.07.pdf", height = 8, width = 8)
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/02_PCA.DimPlot_2023.06.07.pdf", height = 8, width = 8)
DimPlot(obj, reduction = "pca")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/03_PCA1_2023.06.07.pdf", height = 8, width = 8)
DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/04_PCA1.12_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/05_PCA13.24_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/06_PCA25.36_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj, dims = 25:36, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scPCA/07_PCA37.48_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj, dims = 37:48, cells = 500, balanced = TRUE)
dev.off()

# Number of dimensions chosen : 13

# uniform manifold approximation
obj <- RunUMAP(
  obj,
  dims = 1:13,
  n.neighbors = 25,
  n.epochs = 350,
  reduction = "pca",
  seed.use = 42
)

dir.create("Output/Figures/02_scRNA.neftel/scClust")

# clustering
obj <- FindNeighbors(obj, dims = 1:13)
obj.clust <- FindClusters(obj, resolution = seq(0, 1, 0.05), random.seed = 42)

pdf("Output/Figures/02_scRNA.neftel/scClust/01_clusters_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scClust/02_clusters.stability_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Cluster resolution chosen : 0.4

rm(obj.clust)
gc()

obj <- FindClusters(obj, resolution = 0.4, random.seed = 42)
Idents(obj) <- obj@meta.data$seurat_clusters
glimpse(obj@meta.data)
# cluster data for obj w/o peds samples

# plotting

pdf("Output/Figures/02_scRNA.neftel/scUMAP/03_clusters.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "integrated_snn_res.0.4",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/04_Neftel.States.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "cell.state",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/05_CellCycle.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "cell.cycle.phase",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/06_Tumor.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "tumor.id",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/07_Cell.Types.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "CellAssignment",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

################################################################################
# MAIN MARKERS FOR CELL TYPES ####

dir.create("Output/Figures/02_scRNA.neftel/scUMAP/Features")

DefaultAssay(obj) <- "RNA"
obj <- ScaleData(
  obj,
  features = rownames(obj@assays$RNA$data),
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)
# scale.data now on RNA assay now done on all genes

pdf("Output/Figures/02_scRNA.neftel/scUMAP/Features/01_EGFR_2023.11.03.pdf", height = 3.5, width = 3.5)
FeaturePlot(
  obj,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("EGFR"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/Features/02_CD68_2023.11.03.pdf", height = 3.5, width = 3.5)
FeaturePlot(
  obj,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("CD68"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/Features/03_CD3E_2023.11.03.pdf", height = 3.5, width = 3.5)
FeaturePlot(
  obj,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("CD3E"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/Features/04_MBP_2023.11.03.pdf", height = 3.5, width = 3.5)
FeaturePlot(
  obj,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("MBP"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scUMAP/Features/05_MKI67_2023.11.03.pdf", height = 3.5, width = 3.5)
FeaturePlot(
  obj,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("MKI67"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()


################################################################################
# CELL TYPE ANALYSIS ####

# need to re-do PCA and UMAP calculations after having removed pediatric tumors

Idents(obj) <- obj$CellAssignment
nin <- WhichCells(obj, idents = "Neoplastic")
pdf("Output/Figures/02_new_scRNA.neftel/scUMAP/01_cell.type.umap_2023.06.07.pdf", height = 7, width = 7)
DimPlot(
  obj,
  reduction = "umap",
  cells.highlight = nin,
  group.by = "CellAssignment",
  cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  my.theme.umap
dev.off()

pdf("Output/Figures/02_new_scRNA.neftel/scUMAP/01_cell.type.umap_col2_2023.06.07.pdf", height = 7, width = 7)
DimPlot(
  obj,
  reduction = "umap",
  group.by = "CellAssignment",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE, pt.size = 1,
  cols = DiscretePalette(4, palette = "parade")
) +
  labs(title = "") +
  my.theme.umap
dev.off()

pdf("Output/Figures/02_new_scRNA.neftel/01_cell.type.ditto_2023.09.27.pdf", height = 3.5, width = 8.75)
dittoBarPlot(
  obj,
  var = "CellAssignment",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "mako", direction = -1)(4)
) +
  my.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/02_new_scRNA.neftel/01_cell.type.ditto_col2_2023.09.27.pdf", height = 3.5, width = 8.75)
dittoBarPlot(
  obj,
  var = "CellAssignment",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = DiscretePalette(4, palette = "parade")
) +
  my.theme +
  theme(legend.position = "right")
dev.off()

Idents(obj) <- obj$neoplastic.state
pdf("Output/Figures/02_new_scRNA.neftel/scUMAP/02_np.non_col2_2024.02.12.pdf", height = 7, width = 7)
DimPlot(
  obj,
  reduction = "umap",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE, pt.size = 1,
  cols = DiscretePalette(2, palette = "parade")
) +
  my.theme.umap
dev.off()


################################################################################
# ESSENTIAL GENES IN NEFTEL DATASET ####

# violin plot, module score of essential genes, version 1 ----

obj <- AddModuleScore(obj, features = list(killing$Gene), name = "ess.MS_")
# CARS1, EPRS1, SPOUT1, RTF2, SGO1, POLR1H, RARS1, UTP25, NOPCHAP1, NDC1, POLR1G -> missing in scRNA-seq

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_cell.type_2024.02.01.pdf", height = 4, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "CellAssignment",
  cols = viridis_pal(option = "turbo")(4)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# with statistics -- commented out would have been for t-test, need anova
# comp <- list(
#   c("Macrophage", "Neoplastic"),
#   c("Neoplastic", "Oligodendrocyte"),
#   c("Neoplastic", "T-cell")
# )
# ast <- list(
#   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
#   symbols = c("****", "***", "**", "*", "ns")
# )
test.df <- data.frame(
  MS = obj$ess.MS_1,
  CellTypes = obj$CellAssignment
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_cell.type_statistics_2024.02.01.pdf", height = 4, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "CellAssignment",
  cols = viridis_pal(option = "turbo")(4)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.5) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# anova and tukey's for cell type plot

test.res <- aov(MS ~ CellTypes, data = test.df)
print(test.res)

# Call:
#   aov(formula = MS ~ CellTypes, data = test.df)
#
# Terms:
#                 CellTypes Residuals
# Sum of Squares   21.89585 210.64049
# Deg. of Freedom         3      5685
#
# Residual standard error: 0.1924889
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = MS ~ CellTypes, data = test.df)
#
# $CellTypes
# diff          lwr          upr     p adj
# Neoplastic-Macrophage       0.16372051  0.141170423  0.186270598 0.0000000
# Oligodendrocyte-Macrophage -0.04778090 -0.088072061 -0.007489734 0.0124090
# T-cell-Macrophage           0.01245695 -0.046845441  0.071759334 0.9492832
# Oligodendrocyte-Neoplastic -0.21150141 -0.246364875 -0.176637942 0.0000000
# T-cell-Neoplastic          -0.15126356 -0.207020686 -0.095506442 0.0000000
# T-cell-Oligodendrocyte      0.06023784 -0.004752371  0.125228059 0.0806702

# by neftel state ----

test.df <- data.frame(
  MS = obj$ess.MS_1,
  State = obj$cell.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/02_scRNA.neftel/08_essential.genes_module.score_cell.state_statistics_2024.02.01.pdf", height = 4, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "cell.state",
  cols = viridis_pal(option = "turbo")(7)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.5) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/09_essential.genes_module.score_cell.state_statistics_big_2024.02.01.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "cell.state",
  cols = viridis_pal(option = "turbo")(7)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.5) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# anova and tukey's for cell type plot

test.res <- aov(MS ~ State, data = test.df)
print(test.res)

# Call:
#   aov(formula = MS ~ State, data = test.df)
#
# Terms:
#   State Residuals
# Sum of Squares   23.63511 208.90123
# Deg. of Freedom         6      5682
#
# Residual standard error: 0.1917432
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = MS ~ State, data = test.df)
#
# $State
#                             diff          lwr           upr     p adj
# MES1-AC              0.030899128  0.003801622  0.0579966343 0.0136692
# MES2-AC              0.047451867  0.017874167  0.0770295672 0.0000468
# Non-Neoplastic-AC   -0.162621789 -0.186416723 -0.1388268546 0.0000000
# NPC1-AC              0.013277899 -0.011104469  0.0376602658 0.6784091
# NPC2-AC             -0.019367977 -0.044697294  0.0059613405 0.2663221
# OPC-AC               0.038687172  0.011625263  0.0657490810 0.0005025
# MES2-MES1            0.016552739 -0.018844602  0.0519500799 0.8131234
# Non-Neoplastic-MES1 -0.193520917 -0.224250875 -0.1627909593 0.0000000
# NPC1-MES1           -0.017621230 -0.048808266  0.0135658067 0.6386246
# NPC2-MES1           -0.050267105 -0.082199937 -0.0183342736 0.0000719
# OPC-MES1             0.007788044 -0.025535797  0.0411118840 0.9932186
# Non-Neoplastic-MES2 -0.210073656 -0.243011413 -0.1771358997 0.0000000
# NPC1-MES2           -0.034173969 -0.067538572 -0.0008093657 0.0405692
# NPC2-MES2           -0.066819844 -0.100882598 -0.0327570900 0.0000002
# OPC-MES2            -0.008764695 -0.044134793  0.0266054026 0.9907115
# NPC1-Non-Neoplastic  0.175899687  0.147535003  0.2042643713 0.0000000
# NPC2-Non-Neoplastic  0.143253812  0.114071115  0.1724365086 0.0000000
# OPC-Non-Neoplastic   0.201308961  0.170610388  0.2320075338 0.0000000
# NPC2-NPC1           -0.032645875 -0.062309502 -0.0029822491 0.0201661
# OPC-NPC1             0.025409273 -0.005746839  0.0565653856 0.1963128
# OPC-NPC2             0.058055149  0.026152519  0.0899577789 0.0000017

# violin plot, module score of essential genes, version 2 ----
# aggregate tme cell types into "non-neoplastic" to use t-test, simplify fig

comp.2 <- list(
  c("Neoplastic", "Non-Neoplastic")
)
ast <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  symbols = c("****", "***", "**", "*", "ns")
)
test.df <- data.frame(
  MS = obj$ess.MS_1,
  CellTypes = obj$neoplastic.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_npVnon_2024.02.05.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "neoplastic.state",
  cols = viridis_pal(option = "turbo")(4)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.1) +
  stat_compare_means(
    comparisons = comp.2, method = "t.test",
    symnum.args = ast, p.adjust.methods = "bonferroni",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_npVnon_down.size_2024.02.05.pdf", height = 4, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "neoplastic.state",
  cols = viridis_pal(option = "turbo")(4)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.1) +
  stat_compare_means(
    comparisons = comp.2, method = "t.test",
    symnum.args = ast, p.adjust.methods = "bonferroni",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_npVnon_tall.narrow_2024.02.05.pdf", height = 8, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "neoplastic.state",
  cols = viridis_pal(option = "turbo")(4)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.1) +
  stat_compare_means(
    comparisons = comp.2, method = "t.test",
    symnum.args = ast, p.adjust.methods = "bonferroni",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# double-check statistics ----
test.df <- data.frame(
  MS = obj$ess.MS_1,
  NP = obj$neoplastic.state
)
hist(test.df$MS)
dev.off()

test.res <- t.test(MS ~ NP, data = test.df) # two-sided by default
print(test.res$p.value)
# 7.912825e-145

################################################################################

# heatmap of expression of essential genes ----

expr.seurat <- AggregateExpression(
  obj,
  assays = "RNA",
  features = killing$Gene,
  group.by = "CellAssignment",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/02_scRNA.neftel/04_essential.genes_heatmap_cell.type_2024.02.01.pdf", height = 24, width = 16)
DoHeatmap(
  expr.seurat,
  features = killing$Gene,
  group.by = "CellAssignment",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(4)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10)
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/05_essential.genes_heatmap_cell.type_no.gene.names_2024.02.01.pdf", height = 8, width = 4)
DoHeatmap(
  expr.seurat,
  features = killing$Gene,
  group.by = "CellAssignment",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(4)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
dev.off()

# heatmap of expression of essential genes ----
# by Neftel cell states

expr.seurat <- AggregateExpression(
  obj,
  assays = "RNA",
  features = killing$Gene,
  group.by = "cell.state",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/02_scRNA.neftel/06_essential.genes_heatmap_state_2024.02.01.pdf", height = 24, width = 16)
DoHeatmap(
  expr.seurat,
  features = killing$Gene,
  group.by = "cell.state",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(7)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10)
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/07_essential.genes_heatmap_state_no.gene.names_2024.02.01.pdf", height = 8, width = 4)
DoHeatmap(
  expr.seurat,
  features = killing$Gene,
  group.by = "cell.state",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(7)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
dev.off()

################################################################################################################################################################
# REPEAT FOR ONLY NEOPLASTIC CELLS (PROCESSING, ETC.)

dir.create("Output/Figures/02_scRNA.neftel/scNP.only")

Idents(obj) <- obj$neoplastic.state
obj.np <- subset(obj, idents = "Neoplastic")

# re-scale expression data on only neoplastic cells
DefaultAssay(obj.np) <- "integrated"
obj.np <- FindVariableFeatures(obj.np, selection.method = "vst")
# this will have only been done on variable features, not all genes
obj.np <- ScaleData(
  obj.np,
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)

# principal component analysis
obj.np <- RunPCA(obj.np, verbose = FALSE, seed.use = 42)

set.seed(123)
pdf("Output/Figures/02_scRNA.neftel/scNP.only/00_Elbow_2023.06.07.pdf", height = 8, width = 8)
ElbowPlot(obj.np, ndims = 50)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/01_PCA1v2_2023.06.07.pdf", height = 8, width = 8)
VizDimLoadings(obj.np, dims = 1:2, reduction = "pca")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/02_PCA.DimPlot_2023.06.07.pdf", height = 8, width = 8)
DimPlot(obj.np, reduction = "pca", group.by = "cell.state")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/03_PCA1_2023.06.07.pdf", height = 8, width = 8)
DimHeatmap(obj.np, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/04_PCA1.12_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/05_PCA13.24_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/06_PCA25.36_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 25:36, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/07_PCA37.48_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 37:48, cells = 500, balanced = TRUE)
dev.off()

# Number of dimensions chosen : 15

# uniform manifold approximation
obj.np <- RunUMAP(
  obj.np,
  dims = 1:15,
  n.neighbors = 25,
  n.epochs = 350,
  reduction = "pca",
  seed.use = 42
)

# clustering
obj.np <- FindNeighbors(obj.np, dims = 1:15)
obj.clust <- FindClusters(obj.np, resolution = seq(0, 1, 0.05), random.seed = 42)

pdf("Output/Figures/02_scRNA.neftel/scNP.only/01_clusters_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.")
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/02_clusters.stability_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Cluster resolution chosen : 0.3

rm(obj.clust)
gc()

obj.np <- FindClusters(obj.np, resolution = 0.3, random.seed = 42)
Idents(obj.np) <- obj.np@meta.data$seurat_clusters
glimpse(obj.np@meta.data)
# cluster data for obj w/o peds samples

# plotting

pdf("Output/Figures/02_scRNA.neftel/scNP.only/03_clusters.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "integrated_snn_res.0.3",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/04_Neftel.States.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "cell.state",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/05_CellCycle.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "cell.cycle.phase",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/06_Tumor.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "tumor.id",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

################################################################################
# SCALE DATA ON RNA ASSAY ####

DefaultAssay(obj.np) <- "RNA"
obj.np <- ScaleData(
  obj.np,
  features = rownames(obj.np@assays$RNA$data),
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)
# scale.data now on RNA assay now done on all genes

################################################################################
# ESSENTIAL GENES IN NEFTEL DATASET ####

Idents(obj.np) <- obj.np$cell.state
obj.np <- AddModuleScore(obj.np, features = list(killing$Gene), name = "ess.MS_")
# CARS1, EPRS1, SPOUT1, RTF2, SGO1, POLR1H, RARS1, UTP25, NOPCHAP1, NDC1, POLR1G -> missing in scRNA-seq

# violin plot, module score of essential genes, version 1 ----

pdf("Output/Figures/02_scRNA.neftel/scNP.only/03_essential.genes_module.score_cell.state_2024.02.01.pdf", height = 4, width = 4)
VlnPlot(
  obj.np,
  features = "ess.MS_1",
  group.by = "cell.state",
  cols = viridis_pal(option = "turbo")(6)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# with statistics -- commented out would have been for t-test, need anova

test.df <- data.frame(
  MS = obj.np$ess.MS_1,
  State = obj.np$cell.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/03_essential.genes_module.score_cell.type_statistics_2024.02.01.pdf", height = 4, width = 4)
VlnPlot(
  obj.np,
  features = "ess.MS_1",
  group.by = "cell.state",
  cols = viridis_pal(option = "turbo")(6)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], 1.5) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# anova and tukey's for cell type plot

test.res <- aov(MS ~ State, data = test.df)
print(test.res)

# Call:
#   aov(formula = MS ~ State, data = test.df)
#
# Terms:
#   State Residuals
# Sum of Squares    2.20678 185.21578
# Deg. of Freedom         5      4859
#
# Residual standard error: 0.1952385
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = MS ~ State, data = test.df)
#
# $State
#                   diff          lwr          upr     p adj
# MES1-AC    0.029853096  0.003183410  0.056522783 0.0178768
# MES2-AC    0.045431050  0.016320327  0.074541772 0.0001284
# NPC1-AC    0.009191432 -0.014805983  0.033188846 0.8847936
# NPC2-AC   -0.023293190 -0.048222604  0.001636224 0.0828335
# OPC-AC     0.038241933  0.011607282  0.064876585 0.0006161
# MES2-MES1  0.015577954 -0.019260528  0.050416436 0.7989306
# NPC1-MES1 -0.020661664 -0.051356315  0.010032986 0.3904692
# NPC2-MES1 -0.053146286 -0.084574957 -0.021717615 0.0000218
# OPC-MES1   0.008388837 -0.024408881  0.041186555 0.9783898
# NPC1-MES2 -0.036239618 -0.069077455 -0.003401781 0.0206425
# NPC2-MES2 -0.068724239 -0.102249205 -0.035199274 0.0000001
# OPC-MES2  -0.007189116 -0.042000786  0.027622553 0.9918068
# NPC2-NPC1 -0.032484622 -0.061679914 -0.003289329 0.0190085
# OPC-NPC1   0.029050502 -0.001613713  0.059714716 0.0752110
# OPC-NPC2   0.061535123  0.030136177  0.092934069 0.0000004

################################################################################

# heatmap of expression of essential genes ----
# by Neftel cell states

expr.np <- AggregateExpression(
  obj.np,
  assays = "RNA",
  features = killing$Gene,
  group.by = "cell.state",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/02_scRNA.neftel/scNP.only/06_essential.genes_heatmap_state_2024.02.01.pdf", height = 24, width = 16)
DoHeatmap(
  expr.np,
  features = killing$Gene,
  group.by = "cell.state",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(6)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10)
  )
dev.off()

pdf("Output/Figures/02_scRNA.neftel/scNP.only/07_essential.genes_heatmap_state_no.gene.names_2024.02.01.pdf", height = 8, width = 4)
DoHeatmap(
  expr.np,
  features = killing$Gene,
  group.by = "cell.state",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(6)
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
dev.off()


################################################################################
# SAVE ####

saveRDS(
  obj,
  "Output/Rdata/02_scRNA.neftel/01_neftel.noPed_2024.05.02.RDS"
)

saveRDS(
  obj.np,
  "Output/Rdata/02_scRNA.neftel/02_neftel.noPed.NP_2024.05.02.RDS"
)
