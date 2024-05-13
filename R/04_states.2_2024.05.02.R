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
library(singscore)
library(GSEABase)
library(dittoSeq)


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
# # DISCRETE PALETTE PARADE (DISCRETE, QUALITATIVE, ALTERNATE)
# Seurat::DiscretePalette(n_levels, palette = "parade")
## RCOLOR BREWER RDBU (DISCRETE, SIMILAR COLORS)
# colorRampPalette(brewer.pal(n = 11, name = "RdBu"))

dir.create("Output/Rdata/04_states.2_2024.05.02")
dir.create("Output/Figures/04_states.2_2024.05.02")

################################################################################
# GENERATE COLORS FOR HEATMAPS

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
# READ IN DATA

obj <- readRDS("Output/Rdata/02_scRNA.neftel/01_neftel.noPed_2024.05.02.RDS")
DefaultAssay(obj) <- "RNA"

obj.np <- readRDS("Output/Rdata/02_scRNA.neftel/02_neftel.noPed.NP_2024.05.02.RDS")
DefaultAssay(obj.np) <- "RNA"

ess <- read.csv(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv"
)

# list files from 03_SCRIPT output (signatures / tumor)
files <- list.files(
  path = "Output/Rdata/03_states.1_2024.05.02",
  pattern = "^08_filt.fin"
)
length(files) # 11 tumors / 20

# function to create pivoted and binarized data frame from listed files
pivt <- function(frame) {
  fin <- frame %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = k_clust, values_from = present, values_fill = 0) %>%
    as.data.frame()
  rownames(fin) <- fin$gene
  fin$gene <- NULL
  for (col in 1:length(colnames(fin))) {
    fin[, col] <- as.integer(fin[, col])
  }
  return(fin)
}

# for each file, add k_clust column, remove tumor, and bind together with data
# from other files
path = "Output/Rdata/03_states.1_2024.05.02/"
df <- data.frame()
for (i in 1:length(files)) {
  data <- read.csv(
    paste0(path, files[i]),
    row.names = 1
  )
  data$k_clust <- paste0(data$k_clust, "_", data$tumor)
  data$tumor <- NULL
  df <- rbind(df, data)
}
glimpse(df)
table(df$k_clust)

# binarize and pivot data frame for analysis, save
pivot <- pivt(frame = df)
glimpse(pivot)
saveRDS(pivot, "Output/Rdata/04_states.2_2024.05.02/01_pivot.mtx_2024.04.05.RDS")

# jaccard similarities between identified clusters across 11 tumors
# save; this is the true "jaccard matrix" for plotting
to.plot <- t(pivot) %>%
  simil(method = "Jaccard") %>%
  as.matrix()
glimpse(to.plot)
saveRDS(to.plot, "Output/Rdata/04_states.2_2024.05.02/02_to.plot_2024.04.05.RDS")

# 25 clusters from 11 tumors

################################################################################
# CLUSTERING AND PLOTTING

#### IDENTIFY APPROPRIATE K ----------------------------------------------------

dir.create("Output/Figures/04_states.2_2024.05.02/optimal.k")

calc <- t(pivot)
diag(calc) <- 0

set.seed(707)
pdf("Output/Figures/04_states.2_2024.05.02/optimal.k/01_Elbow_Essential.Clusters_2023.09.22.pdf")
fviz_nbclust(calc, kmeans, method = "wss", k.max = 24) +
  theme_minimal() +
  ggtitle("Elbow Method")
dev.off()

set.seed(707)
# gap statistic calculation takes too long when too many samples
gap_stat <- clusGap(calc, FUN = kmeans, K.max = 24, B = 50)
gap <- fviz_gap_stat(gap_stat) +
  theme_minimal() +
  ggtitle("Gap Statistic Method")
pdf("Output/Figures/04_states.2_2024.05.02/optimal.k/02_Gap.Stat_Clusters_2023.09.25.pdf")
gap
dev.off()

#### CLUSTERING AND HEATMAP ----------------------------------------------------

# confirm ward's method by calculating and finding highest agglomerative
# coefficient
ac <- function(x) {
  agnes(pivot, method = x)$ac
}
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
meth <- sapply(m, ac)
meth
# average    single  complete      ward
# 0.9608725 0.9410459 0.9710364 0.9956652

set.seed(707)
clust.rows <- t(pivot) %>%
  dist(method = "binary") %>% # compute dissimilarity matrix (opposite of Jaccard)
  hclust(method = "ward.D2") # hierarchical clustering using Ward's method

# first real cluster split at k = 3 -> 5, 7
ann.col <- data.frame(Cluster = cutree(clust.rows, k = 3))
ann.col$Cluster <- as.character(ann.col$Cluster)

# ann.row <- meta.data[which(rownames(meta.data) %in% rownames(mtx2)), ]
ann <- data.frame(
  clust.ID = rownames(to.plot)
)

splt <- function(x) {
  ugh <- strsplit(x, "_")
  return(ugh[[1]][4])
}

ann$tumor.ID <- sapply(ann$clust.ID, splt)
rownames(ann) <- ann$clust.ID
ann$clust.ID <- NULL

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)

diag(to.plot) <- NA
# generate heatmap object
pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_2023.09.25.pdf"), height = 8, width = 10)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100),
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()


pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_small_2023.09.25.pdf"), height = 4, width = 5)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100),
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()

################################################################################

pdf("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Corrplot_k3_2024.05.08.pdf", height = 8, width = 8)
corrplot(
  to.plot,
  tl.pos = "n",
  order = "hclust",
  hclust.method = "ward.D2",
  method = "shade",
  type = "full", addrect = 3, rect.col = "red",
  col = viridis_pal(option = "inferno", direction = -1)(100),
  is.corr = FALSE,
  diag = TRUE,
  na.label = "square", na.label.col = "grey"
)
dev.off()

################################################################################
#### DIFFERENT K ####

ann.col <- data.frame(Cluster = cutree(clust.rows, k = 2))
ann.col$Cluster <- as.character(ann.col$Cluster)

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)

pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_k2_small_2023.09.25.pdf"), height = 8, width = 10)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100), #cm.colors(100)
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()

#### DIFFERENT K ####

ann.col <- data.frame(Cluster = cutree(clust.rows, k = 4))
ann.col$Cluster <- as.character(ann.col$Cluster)

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)

pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_k4_small_2023.09.25.pdf"), height = 8, width = 10)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100), #cm.colors(100)
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()

#### DIFFERENT K ####

ann.col <- data.frame(Cluster = cutree(clust.rows, k = 5))
ann.col$Cluster <- as.character(ann.col$Cluster)

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)

pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_k5_small_2023.09.25.pdf"), height = 8, width = 10)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100), #cm.colors(100)
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()

#### DIFFERENT K ####

ann.col <- data.frame(Cluster = cutree(clust.rows, k = 7))
ann.col$Cluster <- as.character(ann.col$Cluster)

ann.colors.col <- list(
  tumor.ID = generate.colors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)

pdf(paste0("Output/Figures/04_states.2_2024.05.02/01_Essential.Clusters_Heatmap_k7_small_2023.09.25.pdf"), height = 8, width = 10)
pheatmap(
  to.plot,
  color = viridis_pal(option = "inferno", direction = -1)(100), #cm.colors(100)
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE
)
dev.off()

################################################################################
# CAPTURE VULNERABILITY STATE SIGNATURES

# capture signatures from each cluster
# calculate for k = 2 through k = 8
ann <- as.data.frame(cutree(clust.rows, k = seq(2, 8, 1)))
colnames(ann) <- paste0("k_", colnames(ann))
glimpse(ann)
ann$k_clust <- rownames(ann)
glimpse(ann)

# add "tumor" column to data frame for later merging
glimpse(df)
for (i in 1:length(df$k_clust)) {
  df$tumor[i] <- strsplit(df$k_clust[i], "_")[[1]][4]
}
glimpse(df)

# merge df and ann frames using reduce and dplyr
# creates new -- df with all genes // final signature with tumor and layer 1
# cluster information
lst <- list(df, ann)
new <- lst %>%
  purrr::reduce(full_join, by = "k_clust")
glimpse(new)
new$merge <- paste(new$k_clust, new$gene, sep = "_")

# split new dataframe by k of choice to create 3 signatures
signatures <- split.data.frame(new, new$k_3) # SELECT K OF CHOICE HERE#
glimpse(signatures)
# examine number of genes in each signature and tumors that contribute to each
# signature
for (i in 1:length(signatures)) {
  print(length(rownames(signatures[[i]])))
  print(unique(signatures[[i]]$tumor))
}
# [1] 1772
# [1] "MGH101" "MGH104" "MGH105" "MGH115" "MGH122" "MGH129" "MGH143" "MGH66"
# [1] 1286
# [1] "MGH104" "MGH106" "MGH122" "MGH124" "MGH143" "MGH66"
# [1] 226
# [1] "MGH125" "MGH66"

for (i in 1:length(signatures)) {
  print(unique(signatures[[i]]$k_clust))
}

# create data frame to aggregate non-unique genes w/in each signature by
# avg_log2FC
# need data frame with:
# gene # logFC # k_clust # signature.number -- FETCH FROM 03_SCRIPT OUTPUT

# list all signature files (will have to go through and check whether
# each is in the final sigantures object)
fle <- list.files(
  path = "Output/Rdata/03_states.1_2024.05.02",
  pattern = "06.5_signatures.N5"
)
length(fle)

for (i in 1:length(fle)) {
  data <- read.csv(
    paste0("Output/Rdata/03_states.1_2024.05.02/", fle[i]),
    row.names = 1
  )
  data$merge <- paste(data$k_clust, data$gene, sep = "_")
  for (tm in 1:length(data$merge)) {
    data$tumor[tm] <- strsplit(data$k_clust[tm], "_")[[1]][4]
  }

  for (j in 1:length(signatures)) {
    if (unique(data$tumor) %in% unique(signatures[[j]]$tumor)) {
      data2 <- data[which(data$merge %in% signatures[[j]]$merge), ]

      for (k in 1:length(rownames(signatures[[j]]))) {
        if (signatures[[j]]$merge[k] %in% data2$merge) {
          signatures[[j]]$avg_log2FC[k] <- data2$avg_log2FC[which(data2$merge == signatures[[j]]$merge[k])]
        }
      }
    }
  }
}

glimpse(signatures)
# make sure matches above
for (i in 1:length(signatures)) {
  print(length(rownames(signatures[[i]])))
  print(unique(signatures[[i]]$tumor))
}
# now have added avg_log2FC to the signatures object

saveRDS(
  signatures,
  "Output/Rdata/04_states.2_2024.05.02/04_VS.signatures.all_2024.04.05.RDS"
)

# create final meta-cluster signatures list by aggregating genes by log2FC
# take top 100 genes per signature
sigs <- list()
for (i in 1:length(signatures)) {
  sig <- signatures[[i]] %>%
    aggregate(avg_log2FC ~ gene, FUN = mean)

  nature <- sig %>%
    slice_max(order_by = avg_log2FC, n = 100) # top 100 genes / signature

  sigs[[i]] <- nature

}
glimpse(sigs)
names(sigs) <- c("VS1", "VS2", "VS3")

saveRDS(
  sigs,
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)

# create list of all 300 genes that are in final signature
# (without signature # info -- just vector of genes)
genes <- c()
for (i in 1:length(sigs)) {
  one <- sigs[[i]]$gene
  genes <- c(genes, one)
}
genes

write.csv(
  genes,
  "Output/Rdata/04_states.2_2024.05.02/02_all.VS.markers_noANN_2024.04.05.csv"
)

################################################################################
# SINGSCORE-BASED CELL STATE ASSIGNMENT, ALL CELLS

# add signatures as independent module scores in the object
for (i in 1:length(sigs)) {
  obj.np <- AddModuleScore(
    obj.np,
    features = list(sigs[[i]]$gene),
    name = paste0("VS_", i, "_"),
    slot = "data"
  )
}

# save whole scale.data on all genes separate
obj.data <- GetAssayData(obj.np, assay = "RNA", layer = "scale.data")
saveRDS(obj.data, "Output/Rdata/04_states.2_2024.05.02/05_neftel.noPed.NP_scale.data.ALL_2024.05.05.RDS")

# scale only on genes you are scoring?
# only include non-neoplastic cells when scoring?
obj.np <- ScaleData(
  obj.np,
  features = genes,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
)

# obj.data <- readRDS("Output/Rdata/04_states.2_2024.05.02/05_neftel.noPed.NP_scale.data.ALL_2024.05.05.RDS")
obj.data <- GetAssayData(obj.np, assay = "RNA", layer = "scale.data")
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
rm(obj.data)
rm(obj.data.df)
gc()

# score each cell for each cell state
for (i in 1:length(sigs)) {
  # capture unique gene signature
  sig <- sigs[[i]]$gene
  # capture name of cell state ,
  state <- paste0("Cluster_", i)
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
  "scored.Cluster_1",
  "scored.Cluster_2",
  "scored.Cluster_3"
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
  mVC1 = obj.np@meta.data$Cluster_1singscore,
  mVC2 = obj.np@meta.data$Cluster_2singscore,
  mVC3 = obj.np@meta.data$Cluster_3singscore
)

assignments <- data.frame(
  row.names = rownames(obj.np@meta.data)
)
assignments$cell.state <- "Ambiguous"

for (i in 1:length(rownames(singscore.df))) {
  state <- which.max(singscore.df[i, ])
  assignments$cell.state[i] <- names(state)
}

# add cell states to seurat object (np only)
obj.np <- AddMetaData(obj.np, metadata = assignments, col.name = "mVC")
table(obj.np@meta.data$mVC)

# add cell states to seurat object (all)
obj <- AddMetaData(obj, metadata = assignments, col.name = "mVC")
obj@meta.data <- obj@meta.data %>%
  mutate(
    mVC = case_when(
      neoplastic.state == "Neoplastic" ~ mVC,
      neoplastic.state == "Non-Neoplastic" ~ "Non-Neoplastic"
    )
  )
table(obj@meta.data$mVC)

# dir.create("Output/Figures/04_states.2_2024.05.02/dispersions")
# singscore qc, dispersion figures
pdf("Output/Figures/04_states.2_2024.05.02/dispersions/p01_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.Cluster_1,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/04_states.2_2024.05.02/dispersions/p02_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.Cluster_2,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()
pdf("Output/Figures/04_states.2_2024.05.02/dispersions/p03_singscore.dispersions.QC_2023.11.09.pdf", height = 8, width = 8)
plotDispersion(
  scored.Cluster_3,
  annot = obj.np@meta.data$mVC,
  isInteractive = FALSE
)
dev.off()

################################################################################
# SOME PLOTS WITH CELL STATE ANNOTATIONS DONE

pdf("Output/Figures/04_states.2_2024.05.02/03.5_VC.umap_2023.11.03.pdf", height = 4, width = 4)
DimPlot(
  obj.np,
  reduction = "umap",
  group.by = "mVC",
  # cols.highlight = "#49C1ADFF",
  label = TRUE, label.box = TRUE, repel = TRUE, label.size = 5,
  order = FALSE, raster = FALSE
) +
  labs(title = "") + # fade out all but neoplastic
  umap.theme
dev.off()

pdf("Output/Figures/04_states.2_2024.05.02/04_Essential.Clusters_Barplot_Tumor.NonNP_2023.12.08.pdf", height = 4, width = 8)
dittoBarPlot(
  obj.np,
  var = "mVC",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$mVC)))
) +
  mrd.theme
dev.off()

saveRDS(
  obj,
  "Output/Rdata/04_states.2_2024.05.02/00_neftel.noPed_VC.ANN_2024.04.05.RDS"
)
saveRDS(
  obj.np,
  "Output/Rdata/04_states.2_2024.05.02/00_neftel.noPed.NP_VC.ANN_2024.05.02.RDS"
)

################################################################################
# FINAL DIFFERENTIAL EXPRESSION

Idents(obj.np) <- obj.np$mVC
fin.deg <- FindAllMarkers(
  obj.np,
  assay = "RNA",
  logfc.threshold = 0,
  test.use = "MAST",
  only.pos = TRUE
)

write.csv(
  fin.deg,
  "Output/Rdata/04_states.2_2024.05.02/05_mVC_DEGs_2023.10.17.csv"
)
