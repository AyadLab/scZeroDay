library(tidyverse)
library(Seurat)
library(singscore)
library(dittoSeq)
library(viridis)
library(GSEABase)
library(paletteer)
library(ggpubr)

dir.create("Output/Rdata/08_johnson_2024.05.06")
dir.create("Output/Figures/08_johnson_2024.05.06")

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

john <- readRDS(
  "data/johnson_2023.10.04.RDS"
)
glimpse(john@meta.data)
john@meta.data <- john@meta.data[, -c(70:81)]

ess <- read.csv(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv"
)

################################################################################
# FEATURE PLOTS FOR CELL TYPES

DefaultAssay(john) <- "RNA"

dir.create("Output/Figures/08_johnson_2024.05.06/Features")

pdf("Output/Figures/08_johnson_2024.05.06/Features/01_EGFR_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
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

pdf("Output/Figures/08_johnson_2024.05.06/Features/02_CD68_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
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

pdf("Output/Figures/08_johnson_2024.05.06/Features/03_CD3E_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
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

pdf("Output/Figures/08_johnson_2024.05.06/Features/04_MBP_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
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

pdf("Output/Figures/08_johnson_2024.05.06/Features/05_MKI67_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
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

pdf("Output/Figures/08_johnson_2024.05.06/Features/06_PDGFRB_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("PDGFRB"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/Features/07_PTPRC_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("PTPRC"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/Features/08_CLDN5_2023.11.03.pdf", height = 4, width = 4)
FeaturePlot(
  john,
  reduction = "umap_int",
  slot = "scale.data",
  min.cutoff = 0,
  features = c("CLDN5"),
  order = TRUE, raster = FALSE
) +
  umap.theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, colour = "black")
  )
dev.off()

################################################################################
# CELL TYPE ASSIGNMENTS

DefaultAssay(john) <- "RNA"

glimpse(john@meta.data)

table(john$CellType)
table(john$NPvNon)
table(john$singscorebin) # correct one (?)
table(john$merged_singscorebin)
table(john$cell_state)


john@meta.data <- john@meta.data %>%
  mutate(
    cell.type = case_when(
      cell_state == "Diff.-like" ~ "Neoplastic",
      cell_state == "Prolif. stem-like" ~ "Neoplastic",
      cell_state == "Stem-like" ~ "Neoplastic",
      cell_state == "B cell" ~ cell_state,
      cell_state == "Dendritic cell" ~ cell_state,
      cell_state == "Endothelial" ~ cell_state,
      cell_state == "Fibroblast" ~ cell_state,
      cell_state == "Granulocyte" ~ cell_state,
      cell_state == "Myeloid" ~ cell_state,
      cell_state == "Oligodendrocyte" ~ cell_state,
      cell_state == "Pericyte" ~ cell_state,
      cell_state == "T cell" ~ cell_state,
    )
  )
table(john$cell.type)

john@meta.data <- john@meta.data %>%
  mutate(
    nef.state = case_when(
      cell.type == "Neoplastic" ~ paste0(singscorebin, "_RM"),
      cell.type == "B cell" ~ "Non-Neoplastic",
      cell.type == "Dendritic cell" ~ "Non-Neoplastic",
      cell.type == "Endothelial" ~ "Non-Neoplastic",
      cell.type == "Fibroblast" ~ "Non-Neoplastic",
      cell.type == "Granulocyte" ~ "Non-Neoplastic",
      cell.type == "Myeloid" ~ "Non-Neoplastic",
      cell.type == "Oligodendrocyte" ~ "Non-Neoplastic",
      cell.type == "Pericyte" ~ "Non-Neoplastic",
      cell.type == "T cell" ~ "Non-Neoplastic",
    )
  )
table(john$nef.state)
john@meta.data[which(john$nef.state == "Non-Neoplastic_RM"), ]
john$nef.state[which(john$nef.state == "Non-Neoplastic_RM")] <- "AC_RM"
table(john$nef.state)
john@meta.data <- john@meta.data %>%
  mutate(
    nef.state = case_when(
      nef.state == "AC_RM" ~ "AC",
      nef.state == "MES1_RM" ~ "MES1",
      nef.state == "MES2_RM" ~ "MES2",
      nef.state == "NPC1_RM" ~ "NPC1",
      nef.state == "NPC2_RM" ~ "NPC2",
      nef.state == "OPC_RM" ~ "OPC",
      nef.state == "Non-Neoplastic" ~ "Non-Neoplastic"
    )
  )
table(john$nef.state)

john@meta.data <- john@meta.data %>%
  mutate(
    np.state = case_when(
      cell.type == "Neoplastic" ~ "Neoplastic",
      cell.type == "B cell" ~ "Non-Neoplastic",
      cell.type == "Dendritic cell" ~ "Non-Neoplastic",
      cell.type == "Endothelial" ~ "Non-Neoplastic",
      cell.type == "Fibroblast" ~ "Non-Neoplastic",
      cell.type == "Granulocyte" ~ "Non-Neoplastic",
      cell.type == "Myeloid" ~ "Non-Neoplastic",
      cell.type == "Oligodendrocyte" ~ "Non-Neoplastic",
      cell.type == "Pericyte" ~ "Non-Neoplastic",
      cell.type == "T cell" ~ "Non-Neoplastic",
    )
  )
table(john$np.state)

################################################################################
# BASIC PLOTS

pdf("Output/Figures/08_johnson_2024.05.06/01_cell.type.umap_2023.06.07.pdf", height = 4, width = 4)
DimPlot(
  john,
  reduction = "umap_int",
  group.by = "cell.type",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 2.5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") +
  umap.theme
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/02_cell.type.ditto_2023.09.27.pdf", height = 4, width = 8)
dittoBarPlot(
  john,
  var = "cell.type",
  group.by = "patient.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(length(unique(john$cell.type)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

################################################################################

Idents(john) <- john$nef.state

pdf("Output/Figures/08_johnson_2024.05.06/03_cell.state.umap_2023.06.07.pdf", height = 4, width = 4)
DimPlot(
  john,
  reduction = "umap_int",
  group.by = "nef.state",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 2.5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/04_cell.state.ditto_2023.09.27.pdf", height = 4, width = 8)
dittoBarPlot(
  subset(john, idents = "Non-Neoplastic", invert = TRUE),
  var = "nef.state",
  group.by = "patient.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo", direction = -1)(length(unique(john$nef.state)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

################################################################################
# MODULE SCORE PLOTS

DefaultAssay(john) <- "RNA"
john <- AddModuleScore(john, features = list(ess$Gene), name = "ess.MS_")
# CARS1, EPRS1, NOPCHAP1, POLR1G, POLR1H, RARS1

# Cell Type ----

test.df <- data.frame(
  MS = john$ess.MS_1,
  CellTypes = john$cell.type
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/05_Vln.GDS.MS.Type_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  john,
  features = "ess.MS_1",
  group.by = "cell.type", pt.size = 0.001,
  cols = viridis_pal(option = "turbo", direction = -1)(length(unique(john$cell.type)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(john$ess.MS_1)[1], range(john$ess.MS_1)[2]+0.05) +
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
# Sum of Squares   13.49437  64.92639
# Deg. of Freedom         9     28120
#
# Residual standard error: 0.04805106
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $CellTypes
# diff           lwr          upr     p adj
# Dendritic cell-B cell          -0.004232675 -0.0409733206  0.032507970 0.9999982
# Endothelial-B cell             -0.007510794 -0.0839189707  0.068897383 0.9999995
# Fibroblast-B cell               0.008084705 -0.0591288514  0.075298261 0.9999973
# Granulocyte-B cell             -0.009636460 -0.0449800090  0.025707089 0.9974815
# Myeloid-B cell                 -0.013279260 -0.0482080330  0.021649513 0.9721338
# Neoplastic-B cell               0.028627606 -0.0062647558  0.063519968 0.2194446
# Oligodendrocyte-B cell         -0.041090413 -0.0762021763 -0.005978649 0.0081360
# Pericyte-B cell                 0.004037652 -0.0535571783  0.061632483 1.0000000
# T cell-B cell                   0.001017493 -0.0376880338  0.039723019 1.0000000
# Endothelial-Dendritic cell     -0.003278119 -0.0722382320  0.065681995 1.0000000
# Fibroblast-Dendritic cell       0.012317380 -0.0462910655  0.070925825 0.9996858
# Granulocyte-Dendritic cell     -0.005403785 -0.0183054974  0.007497927 0.9481210
# Myeloid-Dendritic cell         -0.009046585 -0.0207643887  0.002671219 0.3004363
# Neoplastic-Dendritic cell       0.032860281  0.0212514628  0.044469100 0.0000000
# Oligodendrocyte-Dendritic cell -0.036857738 -0.0491102343 -0.024605241 0.0000000
# Pericyte-Dendritic cell         0.008270327 -0.0389996656  0.055540320 0.9999320
# T cell-Dendritic cell           0.005250168 -0.0151313316  0.025631667 0.9983822
# Fibroblast-Endothelial          0.015595498 -0.0734173340  0.104608331 0.9999312
# Granulocyte-Endothelial        -0.002125666 -0.0703516779  0.066100345 1.0000000
# Myeloid-Endothelial            -0.005768466 -0.0737805341  0.062243602 0.9999999
# Neoplastic-Endothelial          0.036138400 -0.0318549757  0.104131776 0.8062541
# Oligodendrocyte-Endothelial    -0.033579619 -0.1016858460  0.034526608 0.8673850
# Pericyte-Endothelial            0.011548446 -0.0704441556  0.093541047 0.9999894
# T cell-Endothelial              0.008528287 -0.0614984226  0.078554996 0.9999970
# Granulocyte-Fibroblast         -0.017721165 -0.0754640541  0.040021725 0.9938320
# Myeloid-Fibroblast             -0.021363964 -0.0788539115  0.036125983 0.9760714
# Neoplastic-Fibroblast           0.020542902 -0.0369249307  0.078010734 0.9816270
# Oligodendrocyte-Fibroblast     -0.049175117 -0.1067764261  0.008426191 0.1730699
# Pericyte-Fibroblast            -0.004047052 -0.0775469797  0.069452875 1.0000000
# T cell-Fibroblast              -0.007067212 -0.0669269864  0.052792563 0.9999977
# Myeloid-Granulocyte            -0.003642800 -0.0096924567  0.002406857 0.6657200
# Neoplastic-Granulocyte          0.038264066  0.0324283070  0.044099826 0.0000000
# Oligodendrocyte-Granulocyte    -0.031453953 -0.0384835396 -0.024424366 0.0000000
# Pericyte-Granulocyte            0.013674112 -0.0325183508  0.059866575 0.9952868
# T cell-Granulocyte              0.010653953 -0.0070857467  0.028393652 0.6691515
# Neoplastic-Myeloid              0.041906866  0.0396912770  0.044122455 0.0000000
# Oligodendrocyte-Myeloid        -0.027811153 -0.0323131347 -0.023309171 0.0000000
# Pericyte-Myeloid                0.017316912 -0.0285589683  0.063192792 0.9734274
# T cell-Myeloid                  0.014296753 -0.0026014629  0.031194968 0.1830821
# Oligodendrocyte-Neoplastic     -0.069718019 -0.0739281944 -0.065507844 0.0000000
# Pericyte-Neoplastic            -0.024589954 -0.0704381180  0.021258210 0.7976937
# T cell-Neoplastic              -0.027610114 -0.0444329380 -0.010787289 0.0000092
# Pericyte-Oligodendrocyte        0.045128065 -0.0008872926  0.091143423 0.0599404
# T cell-Oligodendrocyte          0.042107906  0.0248346190  0.059381192 0.0000000
# T cell-Pericyte                -0.003020159 -0.0518330161  0.045792697 1.0000000

# Neftel Cell State ----

test.df <- data.frame(
  MS = john$ess.MS_1,
  CellState = john$nef.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/06_Vln.GDS.MS.Nef.State_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  john,
  features = "ess.MS_1",
  group.by = "nef.state", pt.size = 0.01,
  cols = viridis_pal(option = "turbo", direction = -1)(length(unique(john$nef.state)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(john$ess.MS_1)[1], range(john$ess.MS_1)[2]+0.05) +
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
# Sum of Squares   15.95015  62.47060
# Deg. of Freedom         6     28123
#
# Residual standard error: 0.04713104
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $CellState
# diff           lwr           upr     p adj
# MES1-AC              0.0100296710  0.0062821456  0.0137771964 0.0000000
# MES2-AC              0.0434726089  0.0400628611  0.0468823567 0.0000000
# Non-Neoplastic-AC   -0.0362362573 -0.0382815307 -0.0341909840 0.0000000
# NPC1-AC              0.0180051426  0.0152286862  0.0207815990 0.0000000
# NPC2-AC              0.0103834815  0.0067792279  0.0139877350 0.0000000
# OPC-AC               0.0050977504  0.0002237225  0.0099717784 0.0334704
# MES2-MES1            0.0334429379  0.0287744753  0.0381114004 0.0000000
# Non-Neoplastic-MES1 -0.0462659284 -0.0500542400 -0.0424776167 0.0000000
# NPC1-MES1            0.0079754715  0.0037473665  0.0122035766 0.0000006
# NPC2-MES1            0.0003538105 -0.0044585493  0.0051661702 0.9999914
# OPC-MES1            -0.0049319206 -0.0107563786  0.0008925374 0.1603740
# Non-Neoplastic-MES2 -0.0797088662 -0.0831633907 -0.0762543418 0.0000000
# NPC1-MES2           -0.0254674663 -0.0293992987 -0.0215356339 0.0000000
# NPC2-MES2           -0.0330891274 -0.0376433825 -0.0285348723 0.0000000
# OPC-MES2            -0.0383748585 -0.0439879454 -0.0327617716 0.0000000
# NPC1-Non-Neoplastic  0.0542413999  0.0514101336  0.0570726662 0.0000000
# NPC2-Non-Neoplastic  0.0466197388  0.0429730962  0.0502663814 0.0000000
# OPC-Non-Neoplastic   0.0413340078  0.0364285509  0.0462394647 0.0000000
# NPC2-NPC1           -0.0076216611 -0.0117233156 -0.0035200066 0.0000009
# OPC-NPC1            -0.0129073921 -0.0181599224 -0.0076548619 0.0000000
# OPC-NPC2            -0.0052857311 -0.0110190554  0.0004475933 0.0937030

################################################################################################################################################################
# NEOPLASTIC-ONLY CLUSTERING / FIGURES

dir.create("Output/Figures/08_johnson_2024.05.06/scNP.only")

Idents(john) <- john$nef.state
obj.np <- subset(john, idents = "Non-Neoplastic", invert = TRUE)

# obj.np <- DietSeurat(
#   obj.np,
#   layers = c("counts", "data", "scale.data"),
#   assays = c("RNA"),
#   dimreducs = names(obj.np@reductions),
#   misc = TRUE
# )
#
# obj.np <- UpdateSeuratObject(obj.np)
# # RNA slot okay
# # integrated slot NOT
# # SCT slot NOT

# re-scale expression data on only neoplastic cells
DefaultAssay(obj.np) <- "RNA"
obj.np <- FindVariableFeatures(obj.np, selection.method = "vst")
obj.np <- PercentageFeatureSet(obj.np, pattern = "^MT-", col.name = "percent.mt")
# this will have only been done on variable features, not all genes
obj.np <- ScaleData(
  obj.np,
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)

# principal component analysis
obj.np <- RunPCA(obj.np, verbose = FALSE, seed.use = 42)

set.seed(123)
pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/00_Elbow_2023.06.07.pdf", height = 8, width = 8)
ElbowPlot(obj.np, ndims = 50)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/01_PCA1v2_2023.06.07.pdf", height = 8, width = 8)
VizDimLoadings(obj.np, dims = 1:2, reduction = "pca")
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/02_PCA.DimPlot_2023.06.07.pdf", height = 8, width = 8)
DimPlot(obj.np, reduction = "pca", group.by = "nef.state")
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/03_PCA1_2023.06.07.pdf", height = 8, width = 8)
DimHeatmap(obj.np, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/04_PCA1.12_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/05_PCA13.24_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/06_PCA25.36_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 25:36, cells = 500, balanced = TRUE)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/07_PCA37.48_2023.06.07.pdf", height = 16, width = 16)
DimHeatmap(obj.np, dims = 37:48, cells = 500, balanced = TRUE)
dev.off()

# Number of dimensions chosen : 19

# uniform manifold approximation
obj.np <- RunUMAP(
  obj.np,
  dims = 1:19,
  n.neighbors = 25,
  n.epochs = 350,
  reduction = "pca",
  seed.use = 42
)

# clustering
obj.np <- FindNeighbors(obj.np, dims = 1:19)
obj.clust <- FindClusters(obj.np, resolution = seq(0, 1, 0.1), random.seed = 42)

library(clustree)
pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/01_clusters_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "RNA_snn_res.")
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/02_clusters.stability_2023.06.07.pdf", height = 18, width = 16)
clustree(obj.clust, prefix = "RNA_snn_res.", node_colour = "sc3_stability")
dev.off()

# Cluster resolution chosen : 0.4

rm(obj.clust)
gc()

obj.np <- FindClusters(obj.np, resolution = 0.4, random.seed = 42)
Idents(obj.np) <- obj.np@meta.data$seurat_clusters
glimpse(obj.np@meta.data)

# plotting

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/03_clusters.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "RNA_snn_res.0.4",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/04_Neftel.States.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "nef.state",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

# pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/05_CellCycle.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
# DimPlot(
#   obj.np,
#   reduction = "umap",
#   # cells.highlight = nin,
#   group.by = "CellCycleIdents",
#   # cols.highlight = "#49C1ADFF",
#   label = TRUE, label.box = TRUE, repel = TRUE, label.size = 3,
#   order = FALSE, raster = FALSE
# ) +
#   labs(title = "") + # fade out all but neoplastic
#   umap.theme
# dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/06_Tumor.umap_2023.11.03.pdf", height = 3.5, width = 3.5)
DimPlot(
  obj.np,
  reduction = "umap",
  # cells.highlight = nin,
  group.by = "patient.ident",
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


john <- AddMetaData(john, metadata = assignments, col.name = "mVC")
john@meta.data <- john@meta.data %>%
  mutate(
    mVC = case_when(
      np.state == "Neoplastic" ~ mVC,
      np.state == "Non-Neoplastic" ~ "Non-Neoplastic"
    )
  )
table(john@meta.data$mVC)

dir.create("Output/Figures/08_johnson_2024.05.06/singscore.disp")
# singscore qc, dispersion figures
pdf("Output/Figures/08_johnson_2024.05.06/singscore.disp/p01_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS1,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/08_johnson_2024.05.06/singscore.disp/p02_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS2,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/08_johnson_2024.05.06/singscore.disp/p03_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.VS3,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/04_Essential.Clusters_Barplot_Tumor.NonNP_2023.12.08.pdf", height = 4, width = 8)
dittoBarPlot(
  john,
  var = "mVC",
  cells.use = nin,
  group.by = "patient.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(john$mVC)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/scNP.only/07_VS.umap_2024.05.06.pdf", height = 3.5, width = 3.5)
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
  MS = john$ess.MS_1,
  State = john$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/08_johnson_2024.05.06/07_Vln.GDS.MS.VS_2024.02.23.pdf", height = 4, width = 4)
VlnPlot(
  john,
  features = "ess.MS_1",
  group.by = "mVC", pt.size = 0.01,
  cols = viridis_pal(option = "turbo")(length(unique(john$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(john$ess.MS_1)[1], range(john$ess.MS_1)[2]+0.05) +
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
# Sum of Squares  17.01297  61.40779
# Deg. of Freedom        3     28126
#
# Residual standard error: 0.04672591
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $State
#                            diff         lwr         upr p adj
# mVC2-mVC1            0.04377172  0.04125086  0.04629259     0
# mVC3-mVC1            0.02374116  0.02130971  0.02617260     0
# Non-Neoplastic-mVC1 -0.01875397 -0.02120403 -0.01630391     0
# mVC3-mVC2           -0.02003057 -0.02193115 -0.01812998     0
# Non-Neoplastic-mVC2 -0.06252569 -0.06445003 -0.06060135     0
# Non-Neoplastic-mVC3 -0.04249513 -0.04430074 -0.04068951     0

################################################################################
# VS module scores

dir.create("Output/Figures/08_johnson_2024.05.06/sigs.by.vs")

for (i in 1:length(markers)) {

  nm <- names(markers)[i]
  vs <- markers[[i]]
  john <- AddModuleScore(john, features = list(vs$gene), name = paste0(nm, "_100_MS_"))

  pdf(paste0("Output/Figures/08_johnson_2024.05.06/sigs.by.vs/06_MS_Vln_", nm, "_2024.02.27.pdf"), height = 7, width = 7)
  print(
    VlnPlot(
      john,
      features = paste0(nm, "_100_MS_1"),
      group.by = "mVC",
      cols = viridis_pal(option = "turbo")(length(unique(john$mVC)))
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

Idents(john) <- john$mVC
pdf("Output/Figures/08_johnson_2024.05.06/06_MSxVS_dotplot_2024.02.27.pdf", height = 4, width = 5.5)
DotPlot(
  subset(john, idents = "Non-Neoplastic", invert = TRUE),
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
jon.deg <- FindAllMarkers(
  obj.np,
  assay = "RNA",
  logfc.threshold = 0,
  test.use = "MAST",
  only.pos = TRUE
)

write.csv(
  jon.deg,
  "Output/Rdata/08_johnson_2024.05.06/02_mVC_DEGs_2024.02.26.csv"
)

################################################################################

saveRDS(
  john,
  "data/johnson_2024.10.04.RDS"
)

saveRDS(
  obj.np,
  "data/johnson.NP_2024.10.04.RDS"
)
