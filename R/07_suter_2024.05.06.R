library(tidyverse)
library(Seurat)
library(singscore)
library(dittoSeq)
library(viridis)
library(GSEABase)
library(paletteer)
library(ggpubr)

dir.create("Output/Rdata/07_suter_2024.05.06")
dir.create("Output/Figures/07_suter_2024.05.06")

mrd.theme <- readRDS("Output/mrd.fig.theme.RDS")
umap.theme <- readRDS("Output/mrd.umap.theme.RDS")

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
# FIGURE SET-UP

generate.colors <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- DiscretePalette(n_levels, palette = "parade")
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

generate.colors_states <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- viridis_pal(option = "turbo")(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

################################################################################
# LOAD DATA

suter <- readRDS(
  "data/Suter_scGBM_2023.10.04.RDS"
)
glimpse(suter@meta.data)

ess <- read.csv(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv"
)

################################################################################
# BASIC PLOTS

Idents(suter) <- suter$NPvNon
nin <- WhichCells(suter, idents = "Neoplastic")

DefaultAssay(suter) <- "integrated"
pdf("Output/Figures/07_suter_2024.05.06/01_cell.type.umap_2023.06.07.pdf", height = 4, width = 4)
DimPlot(
  suter,
  reduction = "umap",
  group.by = "CellType",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 2.5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/02_cell.type.ditto_2023.09.27.pdf", height = 4, width = 8)
dittoBarPlot(
  suter,
  var = "CellType",
  group.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(length(unique(suter$CellType)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

################################################################################

pdf("Output/Figures/07_suter_2024.05.06/03_cell.state.umap_2023.06.07.pdf", height = 4, width = 4)
DimPlot(
  suter,
  reduction = "umap",
  group.by = "ClassificationNPvsNon",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 2.5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/04_cell.state.ditto_2023.09.27.pdf", height = 4, width = 8)
dittoBarPlot(
  subset(suter, idents = "Neoplastic"),
  var = "ClassificationNPvsNon",
  group.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(length(unique(suter$ClassificationNPvsNon)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

################################################################################
# FEATURRE PLOTS FOR CELL TYPES

dir.create("Output/Figures/07_suter_2024.05.06/Features")

pdf("Output/Figures/07_suter_2024.05.06/Features/01_EGFR_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
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

pdf("Output/Figures/07_suter_2024.05.06/Features/02_CD68_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
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

pdf("Output/Figures/07_suter_2024.05.06/Features/03_CD3E_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
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

pdf("Output/Figures/07_suter_2024.05.06/Features/04_MBP_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
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

pdf("Output/Figures/07_suter_2024.05.06/Features/05_MKI67_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
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

pdf("Output/Figures/07_suter_2024.05.06/Features/06_GFAP_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  suter,
  reduction = "umap",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("GFAP"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

################################################################################
# MODULE SCORE PLOTS

DefaultAssay(suter) <- "RNA"
suter <- AddModuleScore(suter, features = list(ess$Gene), name = "ess.MS_")
# CARS1, EPRS1, NOPCHAP1, POLR1G, POLR1H, RARS1

# Cell Type ----

test.df <- data.frame(
  MS = suter$ess.MS_1,
  CellTypes = suter$CellType
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/05_Vln.GDS.MS.Type_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  suter,
  features = "ess.MS_1",
  group.by = "CellType", pt.size = 0.01,
  cols = viridis_pal(option = "turbo", direction = -1)(length(unique(suter$CellType)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(suter$ess.MS_1)[1], range(suter$ess.MS_1)[2]+0.05) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

test.res <- aov(MS ~ CellTypes, data = test.df)
print(test.res)

# Terms:
#   CellTypes Residuals
# Sum of Squares    6.50296 118.01099
# Deg. of Freedom         6     42832
#
# Residual standard error: 0.05249006
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $CellTypes
# diff           lwr           upr     p adj
# Myeloid-Neoplastic           -0.023484991 -0.0255494816 -0.0214205013 0.0000000
# Oligodendrocytes-Neoplastic  -0.048147921 -0.0519456263 -0.0443502161 0.0000000
# Astrocytes-Neoplastic        -0.013628512 -0.0184333057 -0.0088237187 0.0000000
# T-Cells-Neoplastic           -0.020315086 -0.0277515014 -0.0128786708 0.0000000
# Fibroblasts-Neoplastic       -0.003828763 -0.0123659131  0.0047083870 0.8417416
# Endothelial-Neoplastic       -0.009864028 -0.0192024647 -0.0005255913 0.0304679
# Oligodendrocytes-Myeloid     -0.024662930 -0.0288096413 -0.0205162182 0.0000000
# Astrocytes-Myeloid            0.009856479  0.0047713348  0.0149416236 0.0000002
# T-Cells-Myeloid               0.003169905 -0.0044506541  0.0107904647 0.8840999
# Fibroblasts-Myeloid           0.019656228  0.0109582065  0.0283542502 0.0000000
# Endothelial-Myeloid           0.013620963  0.0041352346  0.0231066923 0.0004593
# Astrocytes-Oligodendrocytes   0.034519409  0.0285178143  0.0405210036 0.0000000
# T-Cells-Oligodendrocytes      0.027832835  0.0195724858  0.0360931843 0.0000000
# Fibroblasts-Oligodendrocytes  0.044319158  0.0350554640  0.0535828522 0.0000000
# Endothelial-Oligodendrocytes  0.038283893  0.0282769208  0.0482908657 0.0000000
# T-Cells-Astrocytes           -0.006686574 -0.0154556586  0.0020825108 0.2697850
# Fibroblasts-Astrocytes        0.009799749  0.0000796923  0.0195198060 0.0466208
# Endothelial-Astrocytes        0.003764484 -0.0066663822  0.0141953506 0.9385195
# Fibroblasts-T-Cells           0.016486323  0.0052304918  0.0277421543 0.0003159
# Endothelial-T-Cells           0.010451058 -0.0014240084  0.0223261247 0.1274384
# Endothelial-Fibroblasts      -0.006035265 -0.0186289087  0.0065583789 0.7951008

# Neftel Cell State ----

test.df <- data.frame(
  MS = suter$ess.MS_1,
  CellState = suter$ClassificationNPvsNon
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/06_Vln.GDS.MS.Nef.State_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  suter,
  features = "ess.MS_1",
  group.by = "ClassificationNPvsNon", pt.size = 0.01,
  cols = viridis_pal(option = "turbo", direction = -1)(length(unique(suter$ClassificationNPvsNon)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(suter$ess.MS_1)[1], range(suter$ess.MS_1)[2]+0.05) +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

test.res <- aov(MS ~ CellState, data = test.df)
print(test.res)

# Terms:
#   CellState Residuals
# Sum of Squares   14.17098 110.34297
# Deg. of Freedom         6     42832
#
# Residual standard error: 0.05075609
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $CellState
# diff           lwr          upr     p adj
# NPC2-OPC             0.027630827  0.0240806334  0.031181021 0.0000000
# NPC1-OPC             0.026276130  0.0227693329  0.029782926 0.0000000
# AC-OPC              -0.006624847 -0.0095914916 -0.003658203 0.0000000
# MES2-OPC             0.032859340  0.0295915660  0.036127113 0.0000000
# MES1-OPC            -0.004793660 -0.0077938440 -0.001793476 0.0000506
# Non-Neoplastic-OPC  -0.017324407 -0.0201997409 -0.014449073 0.0000000
# NPC1-NPC2           -0.001354698 -0.0048987044  0.002189309 0.9199995
# AC-NPC2             -0.034255675 -0.0372662125 -0.031245137 0.0000000
# MES2-NPC2            0.005228512  0.0019208388  0.008536186 0.0000646
# MES1-NPC2           -0.032424487 -0.0354680813 -0.029380894 0.0000000
# Non-Neoplastic-NPC2 -0.044955234 -0.0478758344 -0.042034634 0.0000000
# AC-NPC1             -0.032900977 -0.0358602137 -0.029941740 0.0000000
# MES2-NPC1            0.006583210  0.0033221599  0.009844260 0.0000001
# MES1-NPC1           -0.031069790 -0.0340626491 -0.028076930 0.0000000
# Non-Neoplastic-NPC1 -0.043600536 -0.0464682271 -0.040732845 0.0000000
# MES2-AC              0.039484187  0.0368125200  0.042155854 0.0000000
# MES1-AC              0.001831187 -0.0005055977  0.004167972 0.2387907
# Non-Neoplastic-AC   -0.010699559 -0.0128737274 -0.008525392 0.0000000
# MES1-MES2           -0.037653000 -0.0403618612 -0.034944138 0.0000000
# Non-Neoplastic-MES2 -0.050183746 -0.0527536439 -0.047613849 0.0000000
# Non-Neoplastic-MES1 -0.012530747 -0.0147504610 -0.010311032 0.0000000

################################################################################################################################################################
# NEOPLASTIC-ONLY CLUSTERING / FIGURES

dir.create("Output/Figures/07_suter_2024.05.06/scNP.only")

Idents(suter) <- suter$NPvNon
obj.np <- subset(suter, idents = "Neoplastic")

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
pdf("Output/Figures/07_suter_2024.05.06/scNP.only/00_Elbow_2023.06.07.pdf", height = 8, width = 8)
ElbowPlot(obj.np, ndims = 50)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/01_PCA1v2_2023.06.07.pdf", height = 8, width = 8)
VizDimLoadings(obj.np, dims = 1:2, reduction = "pca")
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/02_PCA.DimPlot_2023.06.07.pdf", height = 8, width = 8)
DimPlot(obj.np, reduction = "pca", group.by = "ClassificationNPvsNon")
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/03_PCA1_2023.06.07.pdf", height = 8, width = 8)
DimHeatmap(obj.np, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/04_PCA1.12_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/05_PCA13.24_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/06_PCA25.36_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 25:36, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/07_PCA37.48_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 37:48, cells = 500, balanced = TRUE)
dev.off()

# Number of dimensions chosen : 18

# uniform manifold approximation
obj.np <- RunUMAP(
  obj.np,
  dims = 1:18,
  n.neighbors = 25,
  n.epochs = 350,
  reduction = "pca",
  seed.use = 42
)

# clustering
obj.np <- FindNeighbors(obj.np, dims = 1:18)
obj.clust <- FindClusters(obj.np, resolution = seq(0, 1, 0.1), random.seed = 42)

library(clustree)
pdf("Output/Figures/07_suter_2024.05.06/scNP.only/01_clusters_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.")
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/02_clusters.stability_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Cluster resolution chosen : 0.1

rm(obj.clust)
gc()

obj.np <- FindClusters(obj.np, resolution = 0.1, random.seed = 42)
Idents(obj.np) <- obj.np@meta.data$seurat_clusters
glimpse(obj.np@meta.data)

# plotting

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/03_clusters.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "integrated_snn_res.0.1",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/04_Neftel.States.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "ClassificationNPvsNon",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/05_CellCycle.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "CellCycleIdents",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/06_Tumor.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "orig.ident",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

# SCALE DATA ON RNA ASSAY ####

DefaultAssay(obj.np) <- "RNA"
obj.np <- ScaleData(
  obj.np,
  features = rownames(obj.np@assays$RNA$data),
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)

################################################################################################################################################################
# CELL STATE ASSIGNMENTS

markers <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)
names(markers)

# set up ranked data based on scRNA-seq signatures
obj.data <- GetAssayData(obj.np, assay = "RNA", layer = "scale.data")
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
rm(obj.data)
rm(obj.data.df)
gc()

# score each cell for each cell state
for (i in 1:length(markers)) {
  # capture unique gene signature
  sig <- markers[[i]]$gene
  # capture name of cell state ,
  state <- names(markers)[i]
  # singscore scoring
  set <- GeneSet()
  set@geneIds <- as.character(sig)
  scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
  # rename $Sig to appropriate cell state
  scored$Sig <- as.character(state)
  assign(paste0("scored.", state), scored)
}

# make list of singscores
singscores.list <- c(
  "scored.VS1",
  "scored.VS2",
  "scored.VS3"
)

# append singscores to seurat object
for (n in singscores.list) {
  state <- unlist(strsplit(as.character(n), split = ".", fixed = TRUE))[2]
  scored <- eval(parse(text = n))
  obj.np <- AddMetaData(
    obj.np,
    metadata = scored[, "TotalScore"],
    col.name = paste0(state, "singscore")
  )
}

# add state identities to each cell
singscore.df <- data.frame(
  row.names = rownames(obj.np@meta.data),
  mVC1 = obj.np@meta.data$VS1singscore,
  mVC2 = obj.np@meta.data$VS2singscore,
  mVC3 = obj.np@meta.data$VS3singscore
)

assignments <- data.frame(
  row.names = rownames(obj.np@meta.data)
)
assignments$cell.state <- "Ambiguous"

for (i in 1:length(rownames(singscore.df))) {
  state <- which.max(singscore.df[i, ])
  assignments$cell.state[i] <- names(state)
}

# add cell states to seurat object
obj.np <- AddMetaData(obj.np, metadata = assignments, col.name = "mVC")
table(obj.np@meta.data$mVC)


suter <- AddMetaData(suter, metadata = assignments, col.name = "mVC")
suter@meta.data <- suter@meta.data %>%
  mutate(
    mVC = case_when(
      NPvNon == "Neoplastic" ~ mVC,
      NPvNon == "Non-Neoplastic" ~ "Non-Neoplastic"
    )
  )
table(suter@meta.data$mVC)

# dir.create("Output/Figures/07_suter_2024.05.06/singscore.disp")
# singscore qc, dispersion figures
pdf("Output/Figures/07_suter_2024.05.06/singscore.disp/p01_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS1,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/07_suter_2024.05.06/singscore.disp/p02_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS2,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/07_suter_2024.05.06/singscore.disp/p03_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS3,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/04_Essential.Clusters_Barplot_Tumor.NonNP_2023.12.08.pdf", height = 4, width = 8)
dittoBarPlot(
  suter,
  var = "mVC",
  cells.use = nin,
  group.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(suter$mVC)))
) +
  mrd.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/07_VS.umap_2024.05.06.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  group.by = "mVC",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/scNP.only/07_VS.umap_2024.05.06.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  group.by = "mVC",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

################################################################################
# MODULE SCORES, VS, VIOLIN PLOT

test.df <- data.frame(
  MS = suter$ess.MS_1,
  State = suter$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/07_suter_2024.05.06/07_Vln.GDS.MS.VS_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  suter,
  features = "ess.MS_1",
  group.by = "mVC", pt.size = 0.01,
  cols = viridis_pal(option = "turbo")(length(unique(suter$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(suter$ess.MS_1)[1], range(suter$ess.MS_1)[2]+0.05) +
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

# Terms:
#   State Residuals
# Sum of Squares    8.78373 115.73022
# Deg. of Freedom         3     42835
#
# Residual standard error: 0.05197853
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $State
# diff          lwr         upr p adj
# mVC2-mVC1            0.03082992  0.028652501  0.03300734     0
# mVC3-mVC1            0.00928308  0.007652806  0.01091335     0
# Non-Neoplastic-mVC1 -0.01713007 -0.018820629 -0.01543952     0
# mVC3-mVC2           -0.02154684 -0.023787872 -0.01930581     0
# Non-Neoplastic-mVC2 -0.04795999 -0.050245253 -0.04567474     0
# Non-Neoplastic-mVC3 -0.02641315 -0.028184886 -0.02464142     0

################################################################################
# VS module scores

dir.create("Output/Figures/07_suter_2024.05.06/sigs.by.vs")

for (i in 1:length(markers)) {

  nm <- names(markers)[i]
  vs <- markers[[i]]
  suter <- AddModuleScore(suter, features = list(vs$gene), name = paste0(nm, "_100_MS_"))

  pdf(paste0("Output/Figures/07_suter_2024.05.06/sigs.by.vs/06_MS_Vln_", nm, "_2024.02.27.pdf"), height = 7, width = 7)
  print(
    VlnPlot(
      suter,
      features = paste0(nm, "_100_MS_1"),
      group.by = "mVC",
      cols = viridis_pal(option = "turbo")(length(unique(suter$mVC)))
    ) +
      geom_boxplot() +
      ylab("Module Score \n") +
      stat_compare_means(
        method = "anova"
      ) +
      mrd.theme +
      theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      )
  )
  dev.off()

}

pdf("Output/Figures/07_suter_2024.05.06/06_MSxVS_dotplot_2024.02.27.pdf", height = 4, width = 5.5)
DotPlot(
  subset(suter, idents = "Non-Neoplastic", invert = TRUE),
  features = c("VS1_100_MS_1", "VS2_100_MS_1", "VS3_100_MS_1"),
  group.by = "mVC",
  cols = c("white", "#C91F1F"),
  col.min = 0
) +
  RotatedAxis() +
  coord_flip() +
  mrd.theme +
  theme(
    legend.position = "right"
  )
dev.off()

################################################################################
# differential expression between assigned cell states

Idents(obj.np) <- obj.np@meta.data$mVC
sut.deg <- FindAllMarkers(
  obj.np,
  assay = "RNA",
  logfc.threshold = 0,
  test.use = "MAST",
  only.pos = TRUE
)

write.csv(
  sut.deg,
  "Output/Rdata/07_suter_2024.05.06/02_mVC_DEGs_2024.02.26.csv"
)

################################################################################

saveRDS(
  suter,
  "data/Suter_scGBM_2023.10.04.RDS"
)

saveRDS(
  obj.np,
  "data/Suter_scGBM.NP_2023.10.04.RDS"
)
