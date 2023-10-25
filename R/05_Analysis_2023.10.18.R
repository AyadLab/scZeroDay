library(tidyverse)
library(Seurat)
library(viridis)
library(pheatmap)
library(proxy)
library(corrplot)
library(factoextra)
library(cowplot)
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(reactome.db)
library(dittoSeq)
library(ComplexHeatmap)
library(RColorBrewer)

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
# READ IN DATA

fin.deg <- read.csv(
  "Output/Rdata/05_mVC_DEGs_2023.10.17.csv",
  row.names = 1
)
fin.deg.filt <- fin.deg[which(fin.deg$p_val_adj < 0.05), ]


xClust <- readRDS(
  "Output/Rdata/scRNA/scGBM_xClust.Only_2023.09.25.rds"
)

obj <- readRDS(
  "Output/Rdata/scRNA/scGBM_mVC.Ann_2023.09.28.rds"
)

################################################################################
# EVALUATE mVC CONTENT

pdf("Output/Figures/03_Essential.Clusters_Barplot_CS_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust,
  var = "CellStateClassification4",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(7)
)
dev.off()

pdf("Output/Figures/04_Essential.Clusters_Barplot_Tumor_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust,
  var = "orig.ident",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(8)
)
dev.off()

pdf("Output/Figures/05_Essential.Clusters_Barplot_Phase_2023.09.25.pdf", height = 7, width = 7)
dittoBarPlot(
  xClust,
  var = "Phase",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(8)
)
dev.off()

pdf("Output/Figures/06_Essential.Clusters_Barplot_VC_2023.09.25.pdf", height = 7, width = 14)
dittoBarPlot(
  xClust,
  var = "VC",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(68) # paletteer_c("grDevices::Spectral", 68)
)
dev.off()

pdf("Output/Figures/07_Essential.Clusters_Barplot_VCxTumor_2023.09.25.pdf", height = 7, width = 14)
dittoBarPlot(
  xClust,
  var = "VC",
  group.by = "mVC",
  split.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(68) # paletteer_c("grDevices::Spectral", 68)
)
dev.off()

pdf("Output/Figures/08_Essential.Clusters_Barplot_mVCxTumor_2023.09.25.pdf", height = 7, width = 14)
dittoBarPlot(
  xClust,
  var = "VC",
  group.by = "orig.ident",
  split.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(68) # paletteer_c("grDevices::Spectral", 68)
)
dev.off()

################################################################################
# EXPLORE MVC GENE EXPRESSION PROFILES

pdf("Output/Figures/09_Essential.Clusters_FeaturePlot_EssentialGenes_2023.10.18.pdf", height = 861, width = 49)
FeaturePlot(
  xClust,
  features = ess$Gene,
  split.by = "mVC",
  order = TRUE
)
dev.off()

xClust <- AddModuleScore(xClust, features = list(ess$Gene), name = "ess.MS_")
pdf("Output/Figures/10_Essential.Clusters_Vln_EssentialGene.Module.Score_noNonNeo_2023.10.18.pdf", height = 7, width = 7)
VlnPlot(
  xClust,
  features = "ess.MS_1",
  group.by = "mVC"
) +
  geom_boxplot()
dev.off()

obj <- AddModuleScore(obj, features = list(ess$Gene), name = "ess.MS_")
pdf("Output/Figures/11_Essential.Clusters_Vln_EssentialGene.Module.Score_2023.10.18.pdf", height = 7, width = 7)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "mVC"
) +
  geom_boxplot()
dev.off()

################################################################################
# DIFFERENTIAL EXPRESSION HEATMAP

# how many genes per state were DEG?
table(fin.deg.filt$cluster)

# take top 200 genes for every cluster
def.genes <- fin.deg.filt %>%
  group_by(cluster) %>%
  # Use slice_max to select the top 200 genes within each cluster
  slice_max(order_by = avg_log2FC, n = 200) %>%
  ungroup()

write.csv(
  def.genes,
  "Output/Rdata/06_top.200.mVC_2023.09.26.csv"
)

avg.clust <- AverageExpression(xClust, group.by = "ident", slot = "scale.data")

order_vec <- match(unique(def.genes$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% def.genes$gene), ]
# sort expression by gene

gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% ess$Gene)]
)

# function to add specific gene labels to heatmap (i.e. just essential genes)
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
  row.names = c("mVC1", "mVC2", "mVC3", "mVC4", "mVC5", "mVC6", "mVC7"),
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


pdf("Output/Figures/12_Cluster.Markers_Heatmap_2023.09.26.pdf", height = 7, width = 7)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()

#### FOR ONLY ESSENTIAL GENES --------------------------------------------------

ess.deg <- fin.deg.filt[which(fin.deg.filt$gene %in% ess$Gene), ]
length(unique(ess.deg$gene)) # 66
range(ess.deg$avg_log2FC) # 0.2505186 0.8370831
ess.deg <- ess.deg %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 200) %>% # using 200 to capture all
  ungroup()

order_vec2 <- match(unique(ess.deg$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def2 <- avg.clust$RNA[order_vec2, ]

filt_avg.clust.ess <- sorted_avg.clust.def2[which(rownames(sorted_avg.clust.def2) %in% ess.deg$gene), ]

pdf("Output/Figures/13_Essential.Genes_Heatmap_Clustering_2023.10.18.pdf", height = 10.5, width = 7)
pheatmap::pheatmap(
  filt_avg.clust.ess,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = ess.deg$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  use_raster = TRUE
)
dev.off()

pdf("Output/Figures/13_Essential.Genes_Heatmap_noClustering_2023.10.18.pdf", height = 10.5, width = 7)
pheatmap::pheatmap(
  filt_avg.clust.ess,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = ess.deg$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  use_raster = TRUE
)
dev.off()

################################################################################
# ENRICHMENT ANALYSIS

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

pdf("Output/Figures/02_Essential.Clusters_ClusterProfiler_2023.09.22.pdf", height = 10.5, width = 12)
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

################################################################################
# mVC DIMPLOT

Idents(obj) <- obj$mVC
names <- c("mVC1", "mVC2", "mVC3", "mVC4", "mVC5", "mVC6", "mVC7")
mVC <- list()
for (i in 1:length(names)) {
  cells <- WhichCells(obj, idents = names[i])
  mVC[[length(mVC) + 1]] <- assign(names[i], cells)
  rm(cells)
}

pdf("Output/Figures/14_UMAP_mVC_2023.10.23.pdf", height = 7, width = 7.5)
DimPlot(
  obj,
  group.by = "mVC",
  cells.highlight = mVC,
  cols.highlight = viridis_pal(option = "mako")(7),
  order = TRUE, raster = FALSE
) +
  theme(
    panel.background = element_rect(fill = NA),
    legend.position = "right",
    plot.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
dev.off()
