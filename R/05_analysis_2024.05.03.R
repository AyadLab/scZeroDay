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
library(ggpubr)
library(paletteer)

dir.create("Output/Rdata/05_analysis_2024.05.03")
dir.create("Output/Figures/05_analysis_2024.05.03")

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

fin.deg <- read.csv(
  "Output/Rdata/04_states.2_2024.05.02/05_mVC_DEGs_2023.10.17.csv",
  row.names = 1
)
fin.deg.filt <- fin.deg[which(fin.deg$p_val_adj < 0.05), ]
table(fin.deg.filt$cluster)

obj.np <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/00_neftel.noPed.NP_VC.ANN_2024.05.02.RDS"
)
obj <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/00_neftel.noPed_VC.ANN_2024.04.05.RDS"
)

ess <- read.csv(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv"
)

sigs <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)

################################################################################
# EVALUATE mVC CONTENT

pdf("Output/Figures/05_analysis_2024.05.03/01_Essential.Clusters_Barplot_CS_2023.09.25.pdf", height = 8, width = 8)
dittoBarPlot(
  obj.np,
  var = "cell.state",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$cell.state)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/02_Essential.Clusters_Barplot_Tumor_2023.09.25.pdf", height = 8, width = 8)
dittoBarPlot(
  obj.np,
  var = "tumor.id",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$tumor.id)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/03_Essential.Clusters_Barplot_Tumor_2023.09.25.pdf", height = 8, width = 12)
dittoBarPlot(
  obj.np,
  var = "mVC",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$mVC))+1)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/03_Essential.Clusters_Barplot_Tumor_small_2023.09.25.pdf", height = 4, width = 6)
dittoBarPlot(
  obj.np,
  var = "mVC",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$mVC))+1)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/04_Essential.Clusters_Barplot_Phase_2023.09.25.pdf", height = 8, width = 8)
dittoBarPlot(
  obj.np,
  var = "cell.cycle.phase",
  group.by = "mVC",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$cell.cycle.phase)))
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()


################################################################################
# EXPLORE MVC GENE EXPRESSION PROFILES

pdf("Output/Figures/05_analysis_2024.05.03/05_Essential.Clusters_Vln_Module.Score_NP.only_2023.10.18.pdf", height = 8, width = 8)
VlnPlot(
  obj.np,
  features = "ess.MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj.np$mVC))+1)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/06_Essential.Clusters_Vln_Module.Score_2023.10.18.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# with statistics ----

pdf("Output/Figures/05_analysis_2024.05.03/09_Essential.Clusters_Vln_Module.Score_stats_tall_2023.10.18.pdf", height = 8, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], range(obj$ess.MS_1)[2]+0.2) +
  mrd.theme +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/09_Essential.Clusters_Vln_Module.Score_stats_big_2023.10.18.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], range(obj$ess.MS_1)[2]+0.2) +
  mrd.theme +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/09_Essential.Clusters_Vln_Module.Score_stats_small_2023.10.18.pdf", height = 4, width = 4)
VlnPlot(
  obj,
  features = "ess.MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(obj$ess.MS_1)[1], range(obj$ess.MS_1)[2]+0.2) +
  mrd.theme +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# stats
test.df <- data.frame(
  MS = obj$ess.MS_1,
  VS = obj$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

test.res <- aov(MS ~ VS, data = test.df)
print(test.res)

# Call:
#   aov(formula = MS ~ VS, data = test.df)
#
# Terms:
#   VS Residuals
# Sum of Squares   35.94836 196.58799
# Deg. of Freedom         3      5685
#
# Residual standard error: 0.1859573
# Estimated effects may be unbalanced

# tukey's hsd post-hoc
tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $VS
# diff        lwr         upr     p adj
# mVC2-mVC1            0.132519799  0.1122772  0.15276240 0.0000000
# mVC3-mVC1           -0.002020097 -0.0185414  0.01450121 0.9892597
# Non-Neoplastic-mVC1 -0.148017539 -0.1695303 -0.12650482 0.0000000
# mVC3-mVC2           -0.134539896 -0.1521863 -0.11689344 0.0000000
# Non-Neoplastic-mVC2 -0.280537338 -0.3029257 -0.25814893 0.0000000
# Non-Neoplastic-mVC3 -0.145997442 -0.1650875 -0.12690735 0.0000000


################################################################################
# MAKE VARIOUS DEG CUTOFFS (FOR USE DOWNSTREAM IF NECESSARY)

# how many genes per state were sig DEG?
table(fin.deg.filt$cluster)
fin.deg.filt$cluster <- factor(fin.deg.filt$cluster, levels = c("mVC1", "mVC2", "mVC3"))
# deg cutoffs: p < 0.05, log2FC > 0.5
deg.fc5 <- fin.deg.filt[which(fin.deg.filt$avg_log2FC > 0.5), ]
table(deg.fc5$cluster)

# take top 300 genes for every cluster
def.genes.300 <- deg.fc5 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 300) %>%
  ungroup()
write.csv(
  def.genes.300,
  "Output/Rdata/05_analysis_2024.05.03/01_VC.top.300_2023.09.26.csv"
)

def.genes.200 <- deg.fc5 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 200) %>%
  ungroup()
write.csv(
  def.genes.200,
  "Output/Rdata/05_analysis_2024.05.03/02_VC.top.200_2023.09.26.csv"
)

def.genes.100 <- deg.fc5 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  ungroup()
write.csv(
  def.genes.100,
  "Output/Rdata/05_analysis_2024.05.03/03_VC.top.100_2023.09.26.csv"
)

def.genes.50 <- deg.fc5 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  ungroup()
write.csv(
  def.genes.50,
  "Output/Rdata/05_analysis_2024.05.03/04_VC.top.50_2023.09.26.csv"
)

################################################################################
# DEG HEATMAPS - DOHEATMAP

avg.clust <- AggregateExpression(
  obj.np,
  assays = "RNA",
  features = def.genes.100$gene,
  group.by = "mVC",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/05_analysis_2024.05.03/07_DEG.100.Heatmap_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  avg.clust,
  features = def.genes.100$gene,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank() # element_text(size = 10)
  )
dev.off()

################################################################################
# PHEATMAP

# # NOTE: OBJ.NP SCALE.DATA SLOT NOW CONTAINS SCALED TOP 100 DEG DATA
# DID NOT DO 5/5/2024
# obj.np <- ScaleData(
#   obj.np,
#   assay = "RNA",
#   features = def.genes.100$gene,
#   vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
# )

Idents(obj.np) <- obj.np$mVC
avg.clust <- AverageExpression(obj.np, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(def.genes.100$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% def.genes.100$gene), ]
# sort expression by gene

# none of the essential genes are in the top 100 DEG
gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% ess$Gene)]
)
gene.labels <- c("OVOL2", "SOX10", "GFAP")

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
  row.names = c("mVC1", "mVC2", "mVC3"),
  mVC = c("mVC1", "mVC2", "mVC3")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

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
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/07_DEG.100.pHeatmap_2023.05.05.pdf", height = 8, width = 8)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()

################################################################################
# VS module scores

# top 100 deg module scoring and expr in each cell state

dir.create("Output/Figures/05_analysis_2024.05.03/deg.by.vs.100")

for (i in 1:length(unique(def.genes.100$cluster))) {

  nm <- unique(def.genes.100$cluster)[i]
  vs <- def.genes.100[which(def.genes.100$cluster == nm), ]
  obj.np <- AddModuleScore(obj.np, features = list(vs$gene), name = paste0(nm, "_100_MS_"))

  pdf(paste0("Output/Figures/05_analysis_2024.05.03/ms.by.vs.100/06_MS_Vln_", nm, "_2024.02.27.pdf"), height = 8, width = 8)
  print(
    VlnPlot(
      obj.np,
      features = paste0(nm, "_100_MS_1"),
      group.by = "mVC",
      cols = viridis_pal(option = "turbo")(length(unique(obj.np$mVC)))
    ) +
      geom_boxplot() +
      ylab("Module Score \n") +
      # stat_compare_means(
      #   method = "anova"
      # ) +
      mrd.theme +
      theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      )
  )
  dev.off()

}

pdf("Output/Figures/05_analysis_2024.05.03/08_DEGxVS_100_dotplot_2024.02.27.pdf", height = 3.25, width = 4.75)
DotPlot(
  obj.np,
  features = c("mVC1_100_MS_1", "mVC2_100_MS_1", "mVC3_100_MS_1"),
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

# VS signatures in each cell state

dir.create("Output/Figures/05_analysis_2024.05.03/sigs.by.vs")
for (i in 1:length(sigs)) {

  nm <- names(sigs)[i]
  vs <- sigs[[i]]$gene

  pdf(paste0("Output/Figures/05_analysis_2024.05.03/sigs.by.vs/06_SIGS_Vln_", nm, "_2024.02.27.pdf"), height = 8, width = 8)
  print(
    VlnPlot(
      obj.np,
      features = paste0("VS_", i, "_1"),
      group.by = "mVC",
      cols = viridis_pal(option = "turbo")(length(unique(obj.np$mVC)))
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

pdf("Output/Figures/05_analysis_2024.05.03/08_SIGSxVS_dotplot_2024.02.27.pdf", height = 3.25, width = 4.25)
DotPlot(
  obj.np,
  features = c("VS_1_1", "VS_2_1", "VS_3_1"),
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
# VIOLIN PLOT SHOWING GDS GENES THAT ARE ALSO DEG EXPRESSED IN EACH VS

# genes that are p < 0.05, log2FC > 0
dep.genes <- fin.deg.filt[which(fin.deg.filt$gene %in% ess$Gene), ]

obj <- AddModuleScore(obj, features = list(dep.genes$gene), name = "top.GDS_MS_")

pdf("Output/Figures/05_analysis_2024.05.03/09_VS_Vln_TopGDS.MS_2024.03.25.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "top.GDS_MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# w/ stats ----
test.df <- data.frame(
  MS = obj$top.GDS_MS_1,
  VS = obj$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

test.res <- aov(MS ~ VS, data = test.df)
print(test.res)
# Terms:
#   VS Residuals
# Sum of Squares   43.8833  330.8705
# Deg. of Freedom        3      5685
#
# Residual standard error: 0.241248
# Estimated effects may be unbalanced

tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided

# $VS
#                            diff         lwr         upr     p adj
# mVC2-mVC1            0.16855022  0.14228888  0.19481155 0.0000000
# mVC3-mVC1           -0.00513946 -0.02657305  0.01629413 0.9269531
# Non-Neoplastic-mVC1 -0.13600209 -0.16391118 -0.10809299 0.0000000
# mVC3-mVC2           -0.17368968 -0.19658296 -0.15079640 0.0000000
# Non-Neoplastic-mVC2 -0.30455231 -0.33359746 -0.27550715 0.0000000
# Non-Neoplastic-mVC3 -0.13086263 -0.15562877 -0.10609648 0.0000000

pdf("Output/Figures/05_analysis_2024.05.03/09_VS_Vln_TopGDS.MS_stats_2024.03.25.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "top.GDS_MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

write.csv(
  dep.genes,
  "Output/Rdata/05_analysis_2024.05.03/05_Dep.Sig_TOP_2024.03.25.csv"
)

################################################################################
# HEATMAPS SHOWING GDS GENES THAT ARE ALSO DEG EXPRESSED IN EACH VS

# DoHeatmap ----

dep.genes <- dep.genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 200)

avg.clust <- AggregateExpression(
  obj,
  assays = "RNA",
  features = dep.genes$gene,
  group.by = "mVC",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/05_analysis_2024.05.03/10_TopGDS.Heatmap_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  avg.clust,
  features = dep.genes$gene,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank() # element_text(size = 10)
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/11_TopGDS.Heatmap_single.cell_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  obj,
  features = dep.genes$gene,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank() # element_text(size = 10)
  )
dev.off()

# pHeatmap ----

# DON'T NEED TO DO? (DID NOT DO 5/5/2025)
# obj <- ScaleData(
#   obj,
#   assay = "RNA",
#   features = dep.genes$gene,
#   vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
# )

Idents(obj) <- obj$mVC
avg.clust <- AverageExpression(obj, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(dep.genes$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% dep.genes$gene), ]
# sort expression by gene

# none of the essential genes are in the top 100 DEG
gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% ess$Gene)]
)
gene.labels <- c("CDC42", "MTOR", "CDK2", "RCC1", "ITGAV1", "LMNA")

an.col <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic"),
  mVC = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

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
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/12_TopGDS.pHeatmap_2024.05.05.pdf", height = 8, width = 8)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()

# clustering and all labels

pdf("Output/Figures/05_analysis_2024.05.03/13_TopGDS.pHeatmap_clustered_2024.05.05.pdf", height = 12, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/14_TopGDS.pHeatmap_not.clustered_2024.05.05.pdf", height = 12, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

# pHeatmap, non-neoplastic cells only ----

Idents(obj.np) <- obj.np$mVC
avg <- AverageExpression(obj.np, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(dep.genes$gene), rownames(avg$RNA))
sorted_avg.clust.def <- avg$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% dep.genes$gene), ]
# sort expression by gene

an.col <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3"),
  mVC = c("mVC1", "mVC2", "mVC3")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

pdf("Output/Figures/05_analysis_2024.05.03/14_TopGDS.pHeatmap_NP.only_2024.05.05.pdf", height = 16, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

# pHeatmap, non-neoplastic cells only, all essential genes ----

order_vec <- match(unique(ess$Gene), rownames(avg$RNA))
sorted_avg.clust.def <- avg$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% ess$Gene), ]
# sort expression by gene

an.col <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3"),
  mVC = c("mVC1", "mVC2", "mVC3")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

pdf("Output/Figures/05_analysis_2024.05.03/14_AllGDS.pHeatmap_NP.only_2024.05.05.pdf", height = 16, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()


################################################################################
# LOOKING AT ESSENTIAL GENES THAT ARE W/IN WHOLE VS SIGNATURES

signatures <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/04_VS.signatures.all_2024.04.05.RDS"
)

wholesig.genes <- c()
for (i in 1:length(signatures)) {
  gn <- signatures[[i]]$gene
  wholesig.genes <- c(wholesig.genes, gn)
}

sig.ess <- unique(wholesig.genes)[which(unique(wholesig.genes) %in% ess$Gene)]
write.csv(
  sig.ess,
  "Output/Rdata/05_analysis_2024.05.03/06_sig.GDS_2024.04.05.csv"
)

obj <- AddModuleScore(obj, features = list(sig.ess), name = "Sig.GDS_MS_")

pdf("Output/Figures/05_analysis_2024.05.03/15_VS_Vln_SigGDS.MS_2024.03.25.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  group.by = "mVC",
  features = "Sig.GDS_MS_1",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

# w/ stats ----
test.df <- data.frame(
  MS = obj$Sig.GDS_MS_1,
  VS = obj$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

test.res <- aov(MS ~ VS, data = test.df)
print(test.res)
# Terms:
#   VS Residuals
# Sum of Squares   74.7823  552.3235
# Deg. of Freedom        3      5685
#
# Residual standard error: 0.3116962
# Estimated effects may be unbalanced

tukey <- TukeyHSD(test.res)
print(tukey) # bonferroni-adjusted p-values provided
# $VS
# diff         lwr         upr     p adj
# mVC2-mVC1            0.1989020  0.16497198  0.23283208 0.0000000
# mVC3-mVC1           -0.0172986 -0.04499113  0.01039394 0.3756582
# Non-Neoplastic-mVC1 -0.2022349 -0.23829390 -0.16617593 0.0000000
# mVC3-mVC2           -0.2162006 -0.24577909 -0.18662216 0.0000000
# Non-Neoplastic-mVC2 -0.4011369 -0.43866374 -0.36361016 0.0000000
# Non-Neoplastic-mVC3 -0.1849363 -0.21693457 -0.15293807 0.0000000

pdf("Output/Figures/05_analysis_2024.05.03/16_VS_Vln_SigGDS.MS_stats_2024.03.25.pdf", height = 8, width = 8)
VlnPlot(
  obj,
  features = "Sig.GDS_MS_1",
  group.by = "mVC",
  cols = viridis_pal(option = "turbo")(length(unique(obj$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  mrd.theme +
  stat_compare_means(
    method = "anova",
    step.increase = 0.125) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
dev.off()

################################################################################
# HEATMAPS SHOWING GDS GENES THAT ARE ALSO VS SIGNATURE GENES

# DoHeatmap

avg.clust <- AggregateExpression(
  obj,
  assays = "RNA",
  features = sig.ess,
  group.by = "mVC",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/05_analysis_2024.05.03/17_SigGDS.Heatmap_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  avg.clust,
  features = sig.ess,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
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

pdf("Output/Figures/05_analysis_2024.05.03/18_SigGDS.Heatmap_single.cell_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  obj,
  features = sig.ess,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
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

# pHeatmap ----

# DON'T NEED TO DO? (DID NOT DO 5/5/2025)
# obj <- ScaleData(
#   obj,
#   assay = "RNA",
#   features = sig.ess,
#   vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
# )

Idents(obj) <- obj$mVC
avg.clust <- AverageExpression(obj, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(sig.ess), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% sig.ess), ]
# sort expression by gene

an.col <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic"),
  mVC = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

# clustering and all labels

pdf("Output/Figures/05_analysis_2024.05.03/19_SigGDS.pHeatmap_clustered_2024.05.05.pdf", height = 12, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/20_SigGDS.pHeatmap_not.clustered_2024.05.05.pdf", height = 12, width = 8)
pheatmap::pheatmap(
  filt_avg.clust.def,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # annotation_row = def.genes$gene,
  annotation_col = an.col,
  annotation_colors = an.colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

################################################################################\
# VS MARKER GENE EXPRESSION / MODULE

# evaluate number of clusters w/in each VS ----
c <- c()
for (i in 1:length(signatures)) {
  b <- unique(signatures[[i]]$k_clust)
  c <- c(c, b)
  c <- unique(c)
}
c # there are 25 unique clusters w/in signatures

# make exhaustive list of all top 100 marker genes per state, including ----
# duplicated genes per signature
mark.cor.list <- list()
for (i in 1:length(signatures)) {
  sub <- signatures[[i]][which(signatures[[i]]$gene %in% sigs[[i]]$gene), ]
  sub2 <- sub[, -c(4, 6:10)]
  mark.cor.list[[i]] <- sub2
}
glimpse(mark.cor.list)

c <- c()
for (i in 1:length(mark.cor.list)) {
  b <- unique(mark.cor.list[[i]]$k_clust)
  c <- c(c, b)
  c <- unique(c)
}
c # now only 24 clusters -- meaning for the top 100 genes / signature, 1 cluster is not represented

# confirm all genes used as 100-gene signatures are in this mark.cor.list ----

all.genes <- c()
for (i in 1:length(mark.cor.list)) {
  g <- unique(mark.cor.list[[i]]$gene)
  all.genes <- c(all.genes, g)
  all.genes <- unique(all.genes)
}

genes <- read.csv( # read in 300-gene list
  "Output/Rdata/04_states.2_2024.05.02/02_all.VS.markers_noANN_2024.04.05.csv",
  row.names = 1
)
genes <- genes$x
unsig.genes <- unique(genes)

# confirm match
length(unsig.genes[unsig.genes %in% all.genes]) # 292

# add average expression of genes in cells w/in k_clust ----

files <- list.files(
  path = "Output/Rdata/03_states.1_2024.05.02",
  pattern = "^03_cluster.data"
)
path = "Output/Rdata/03_states.1_2024.05.02/"

# start the loop for extraction of expression data
expr.df <- data.frame()
final <- data.frame()
for (i in 1:length(mark.cor.list)) {
  tums <- unique(mark.cor.list[[i]]$tumor)
  data.list <- list()
  for (t in 1:length(tums)) {
    file <- files[grep(tums[t], files)]
    data <- read.csv(
      paste0(path, file),
      row.names = 1
    )
    data.list[[t]] <- data
  }
  names(data.list) <- tums

  norm.list <- list()
  for (clust in 1:length(unique(mark.cor.list[[i]]$k_clust))) {
    k_clus <- unique(mark.cor.list[[i]]$k_clust)[clust]
    pre <- strsplit(unique(mark.cor.list[[i]]$k_clust)[clust], "_")
    tumor <- pre[[1]][4]
    k1 <- pre[[1]][2]
    k <- pre[[1]][3]
    dat.inx <- which(names(data.list) == tumor)
    index.k <- which(colnames(data.list[[dat.inx]]) == paste0("k_", k1))

    genes <- all.genes
    cells <- rownames(data.list[[dat.inx]])[which(data.list[[dat.inx]][, index.k] == k)]

    sub.obj <- subset(obj.np, cells = cells)

    norm <- sub.obj@assays$RNA$data[which(rownames(sub.obj@assays$RNA$data) %in% genes), ]

    if (length(rownames(norm)) > 1) {
      norm.avg <- rowMeans(norm)
    }
    else {
      norm.avg <- mean(norm)
      names(norm.avg) <- genes
    }

    norm.list[[clust]] <- norm.avg

  }
  names(norm.list) <- unique(mark.cor.list[[i]]$k_clust)

  df <- data.frame()
  for (sig in 1:length(norm.list)) {
    clus <- names(norm.list)[sig]
    merge <- paste0(names(norm.list[sig]), "_", names(norm.list[[sig]]))
    new <- data.frame(
      genes = names(norm.list[[sig]]),
      avg.norm.expr = as.double(norm.list[[sig]]),
      merge = merge,
      clust = clus
    )
    df <- rbind(df, new)
  }
  df2 <- df
  df$clust <- NULL
  df$genes <- NULL
  fin <- merge(mark.cor.list[[i]], df, by = "merge")
  expr.df <- rbind(expr.df, df2)
  final <- rbind(final, fin)
}

glimpse(final)
glimpse(expr.df)

# plot expression heatmap ----

# clean up data
expr.df$merge <- NULL
mtx <- expr.df %>%
  pivot_wider(names_from = clust, values_from = avg.norm.expr, values_fill = NA) %>%
  as.data.frame()
# confirm no NA
mtx[is.na(mtx)]
rownames(mtx) <- mtx$genes
mtx$genes <- NULL
for (col in 1:length(colnames(mtx))) {
  mtx[, col] <- as.double(mtx[, col])
}
glimpse(mtx)

# function for string extraction
fun <- function(x) {
  word <- strsplit(x, "_")[[1]][4]
  return(word)
}

# start plot
ann <- data.frame(
  clust.ID = colnames(mtx),
  tumor = sapply(colnames(mtx), fun)
)
smol <- final[, -c(1:2, 4, 6:7)]
smol <- smol %>%
  distinct()

colnames(smol) <- c("clust.ID", "clust")
ann2 <- merge(ann, smol, by = "clust.ID", all = FALSE)
glimpse(ann2)
colnames(ann2) <- c("clust.ID", "tumor.ID", "state")
rownames(ann2) <- ann2$clust.ID
ann2$clust.ID <- NULL
ann2$state <- factor(ann2$state, levels = c("1", "2", "3"))

ann.colors.col <- list(
  # clust.ID = generate.colors(ann2)[[1]],
  tumor.ID = generate.colors(ann2)[[1]],
  state = generate.colors_states(ann2)[[2]]
)

# need to assign each gene to a signature
smol2 <- final[, -c(1, 3:4, 6:7)]
smol2 <- smol2 %>%
  distinct()
# remove duplicates
dup <- which(duplicated(smol2$gene))
smol2 <- smol2[-dup, ]
colnames(smol2) <- c("gene", "state")
rownames(smol2) <- smol2$gene
smol2$gene <- NULL

scaled <- ScaleData(mtx)

my.breaks <- c(seq(min(scaled), 0, length.out=ceiling(100/2) + 1),
               seq(max(scaled)/100, max(scaled), length.out=floor(100/2)))

# generate heatmap object
pdf(paste0("Output/Figures/05_analysis_2024.05.03/21_VS.Markers.Heatmap_unsorted_2023.09.25.pdf"), height = 12, width = 8)
pheatmap::pheatmap(
  as.matrix(scaled),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  annotation_row = smol2,
  annotation_col = ann2,
  annotation_colors = ann.colors.col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

# need to order genes appropriately (order by signature)
markers <- sigs %>%
  purrr::reduce(full_join, by = "gene")
hmm <- scaled[markers$gene, ] # markers should be ordered by log2FC
# hmm <- hmm[rownames(smol2), ] # smol is ordered by tumor

pdf(paste0("Output/Figures/05_analysis_2024.05.03/22_VS.Markers.Heatmap_sorted_2023.09.25.pdf"), height = 12, width = 8)
pheatmap::pheatmap(
  as.matrix(hmm),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)), #paletteer_c("ggthemes::Red-Blue-White Diverging", 100, direction = -1)
  annotation_row = smol2,
  annotation_col = ann2,
  annotation_colors = ann.colors.col,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf(paste0("Output/Figures/05_analysis_2024.05.03/23_VS.Markers.Heatmap_sorted.noclust_2023.09.25.pdf"), height = 12, width = 8)
pheatmap::pheatmap(
  as.matrix(hmm),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)), #paletteer_c("ggthemes::Red-Blue-White Diverging", 100, direction = -1)
  annotation_row = smol2,
  # annotation_col = ann2,
  annotation_colors = ann.colors.col[2],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

################################################################################
# ALL DEG LOGFC > 0, P < 0.05 AS HEATMAP (TOP 1500 genes)
# 1500 better visually

# DoHeatmap ----

to.hm <- fin.deg.filt %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1500)

avg.clust <- AggregateExpression(
  obj,
  assays = "RNA",
  features = to.hm$gene,
  group.by = "mVC",
  normalization.method = "CLR",
  return.seurat = TRUE
)

pdf("Output/Figures/05_analysis_2024.05.03/24_1500.DEG.Heatmap_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  avg.clust,
  features = to.hm$gene,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank() # element_text(size = 10)
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/25_1500.DEG.Heatmap_single.cell_2024.05.05.pdf", height = 10, width = 5)
DoHeatmap(
  obj,
  features = to.hm$gene,
  group.by = "mVC",
  draw.lines = FALSE, group.colors = viridis_pal(option = "turbo")(length(unique(avg.clust$mVC)))
) +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank() # element_text(size = 10)
  )
dev.off()

# pHeatmap ----

# DON'T NEED TO DO? (DID NOT DO 5/5/2025)
# obj <- ScaleData(
#   obj,
#   assay = "RNA",
#   features = to.hm$gene,
#   vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
# )

Idents(obj) <- obj$mVC
avg.clust <- AverageExpression(obj, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(to.hm$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% to.hm$gene), ]
# sort expression by gene

# none of the essential genes are in the top 100 DEG
gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% ess$Gene)]
)

an.col <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic"),
  mVC = c("mVC1", "mVC2", "mVC3", "Non-Neoplastic")
)
an.colors <- list(mVC = generate.colors_states(an.col)[[1]])

my.breaks <- c(seq(min(filt_avg.clust.def), 0, length.out=ceiling(100/2) + 1),
               seq(max(filt_avg.clust.def)/100, max(filt_avg.clust.def), length.out=floor(100/2)))

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
  annotation_names_col = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/26_1500.DEG.pHeatmap_2024.05.05.pdf", height = 12, width = 8)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()

################################################################################
# IMPROVED ENRICHMENT ANALYSIS X2 FINAL FINAL FINAL FINAL - ON MARKERS

# H, C2, C3 ----

glimpse(sigs)

markers.list <- list(
  VS1 = sigs[[1]]$gene,
  VS2 = sigs[[2]]$gene,
  VS3 = sigs[[3]]$gene
)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb2 <- msigdbr(species = "Homo sapiens", category = "H")
msigdb6 <- msigdbr(species = "Homo sapiens", category = "C3")

msigdb7 <- rbind(msigdb, msigdb2, msigdb6)

msigdb_ref <- msigdb7 %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  markers.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

glimpse(comp@compareClusterResult)
table(comp@compareClusterResult$Cluster)

View(comp@compareClusterResult)

# remove irrelevant enrichment results
filt <- comp[-grep("MCLACHLAN", comp@compareClusterResult$ID), ]
filt <- filt[-grep("ROSTY", filt$ID), ]
filt <- filt[-grep("DUTERTRE", filt$ID), ]
filt <- filt[-grep("BLANCO", filt$ID), ]
filt <- filt[-grep("UNKNOWN", filt$ID), ]
filt <- filt[-grep("MIR33$", filt$ID), ]
filt <- filt[-grep("KOBAYASHI", filt$ID), ]
filt <- filt[-grep("SOTIRIOU", filt$ID), ]
filt <- filt[-grep("SHEDDEN", filt$ID), ]
filt <- filt[-grep("HORIUCHI", filt$ID), ]
filt <- filt[-grep("CHIANG", filt$ID), ]
filt <- filt[-grep("24HR$", filt$ID), ]
filt <- filt[-grep("HOFFMAN", filt$ID), ]
filt <- filt[-grep("KONG", filt$ID), ]

comp.filt <- comp
comp.filt@compareClusterResult <- filt
# glimpse(comp.filt@compareClusterResult)

dir.create("Output/Figures/05_analysis_2024.05.03/marker_enrich")

pdf("Output/Figures/05_analysis_2024.05.03/marker_enrich/01_msigdb.C2.C3.H_filt_ORA_2024.04.24.pdf", height = 8, width = 10)
dotplot(
  comp.filt,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
  label_format = 50
) +
  # facet_grid("DB") +
  theme_bw() +
  # geom_text(aes(label = GeneRatio), hjust = -0.2, colour = "black") +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, colour = "black")
  )
dev.off()

# C5 ----

msigdb <- msigdbr(species = "Homo sapiens", category = "C5")
# msigdb2 <- msigdbr(species = "Homo sapiens", category = "C2")
# msigdb3 <- msigdbr(species = "Homo sapiens", category = "C3")

# msigdb7 <- rbind(msigdb, msigdb3)

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  markers.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
# table(comp@compareClusterResult$DB)

glimpse(comp@compareClusterResult)
table(comp@compareClusterResult$Cluster)

# filt <- comp[-grep("MCLACHLAN", comp@compareClusterResult$ID), ]
#
# comp.filt <- comp
# comp.filt@compareClusterResult <- filt
#
# View(comp.filt@compareClusterResult)

pdf("Output/Figures/05_analysis_2024.05.03/marker_enrich/02_msigdb.C5_filt_ORA_2024.04.24.pdf", height = 8, width = 10)
dotplot(
  comp.filt,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
  label_format = 50
) +
  # facet_grid("DB") +
  theme_bw() +
  # geom_text(aes(label = GeneRatio), hjust = -0.2, colour = "black") +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, colour = "black")
  )
dev.off()

################################################################################
# ANNOTATE AND SAVE NEFTEL OBJ

obj@meta.data <- obj@meta.data %>%
  mutate(
    mVC.types = case_when(
      mVC == "Non-Neoplastic" ~ deconv,
      mVC != "Non-Neoplastic" ~ mVC
    )
  )

saveRDS(
  obj,
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed_VC.ANN_2024.05.06.RDS"
)
saveRDS(
  obj.np,
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed.NP_VC.ANN_2024.05.06.RDS"
)
