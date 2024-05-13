library(tidyverse)
library(ggpubr)
library(Seurat)
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
library(viridis)
library(RColorBrewer)
library(paletteer)
library(SPOTlight)

mrd.theme <- readRDS("Output/mrd.fig.theme.RDS")
umap.theme <- readRDS("Output/mrd.umap.theme.RDS")

dir.create("Output/Figures/XX_MAIN.FIGS.enhanced")

################################################################################

generate.colors_tumors <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- colorRampPalette(paletteer_d(palette = "tidyquant::tq_light", 12))(n_levels)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

generate.colors_states <- function(df) {
  myne <- list()
  for (vec in 1:length(colnames(df))) {
    n_levels <- length(unique(df[, vec]))
    colors <- colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(4)
    names(colors) <- unique(df[, vec])
    myne[[length(myne) + 1]] <- colors
  }
  return(myne)
}

################################################################################
# FIGURE 1 A, B, C

efx <- read.delim(
  file = "data/r.input_CRISPRGeneEffect_23Q2.csv",
  sep = ",",
  row.names = 1,
  check.names = FALSE
)
colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))

gbm.stats_m <- read.csv(
  gbm.stats_m,
  file = "Output/Rdata/01_essential.genes/02_depmap.gbm.stats_mean.effect.scores_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# effect scores for all genes in all cell lines
all <- t(as.data.frame(lapply(efx, mean)))
all <- as.data.frame(all)
colnames(all) <- c("mEffectScore")

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_mean.effect.score_all.genes_all.lines_2024.02.07.pdf", height = 4, width = 4)
all %>%
  ggplot(
    aes(x = mEffectScore, y = after_stat(density)) # compute kernel density estimates and replaces counts w/ density in each bucket
  ) +
  geom_density(color = "#053061", fill = "#A7CFE4") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  xlim(-3, 0.25) +
  scale_x_continuous(n.breaks = 8) +
  theme_classic() +
  mrd.theme
dev.off()

# add Essential / Not Essential annotation column
gbm.stats_m$Essential <- "Non-Essential"
gbm.stats_m$Essential[which(gbm.stats_m$Gene %in% killing$Gene)] <- "Essential"
gbm.stats_m$Essential <- as.factor(gbm.stats_m$Essential)

# mean effect scores for all genes in GBM cell lines
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_mean.effect.score_all.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
gbm.stats_m %>%
  ggplot(
    aes(x = mEffectScore, y = after_stat(density)) # compute kernel density estimates and replaces counts w/ density in each bucket
  ) +
  geom_density(color = "#053061", fill = "#A7CFE4") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  xlim(-3, 0.25) +
  theme_classic() + # theme_bw()
  mrd.theme
dev.off()

# mean effect scores for all genes significantly lower in GBM cell ines
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_mean.effect.score_sig.lower.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
gbm.stats_m[which(gbm.stats_m$q.value < 0.05 & gbm.stats_m$EffectSize < 0), ] %>%
  ggplot(
    aes(x = mEffectScore, y = after_stat(density))
  ) +
  geom_density(color = "#053061", fill = "#A7CFE4") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  xlim(-3, 0.25) +
  theme_classic() + # theme_bw()
  mrd.theme
dev.off()

# mean effect scores for all essential genes
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_mean.effect.score_killing.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
gbm.stats_m[which(gbm.stats_m$Essential == "Essential" ), ] %>%
  ggplot(
    aes(x = mEffectScore, y = after_stat(density))
  ) +
  geom_density(color = "#053061", fill = "#A7CFE4") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  xlim(-3, -0.5) +
  theme_classic() + # theme_bw()
  mrd.theme
dev.off()

################################################################################
# FIGURE 1 D

nef <- readRDS(
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed_VC.ANN_2024.05.06.RDS"
)

comp.2 <- list(
  c("Neoplastic", "Non-Neoplastic")
)
ast <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  symbols = c("****", "***", "**", "*", "ns")
)
test.df <- data.frame(
  MS = nef$ess.MS_1,
  CellTypes = nef$neoplastic.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_essential.genes_module.score_npVnon_2024.02.05.pdf", height = 8, width = 6)
VlnPlot(
  nef,
  features = "ess.MS_1",
  group.by = "neoplastic.state",
  pt.size = 0.01,
  cols = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(2)
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(nef$ess.MS_1)[1], 1.1) +
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

################################################################################
# FIGURE 1 E

killing <- read.table(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv",
  row.names = 1, header = TRUE, sep = ","
)

killing <- killing %>%
  arrange(mEffectScore)
glimpse(killing)

gs.list <- list(
  ess = killing$Gene
)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb2 <- msigdbr(species = "Homo sapiens", category = "H")

msigdb3 <- rbind(msigdb, msigdb2)

msigdb_ref <- msigdb3 %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  gs.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

funni <- function(X) {
  db <- str_split(X, pattern = "_")[[1]][1]
  return(db)
}

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

# filter only for databases of interest
filt <- comp[which(comp@compareClusterResult$DB == "REACTOME" | comp@compareClusterResult$DB == "KEGG" | comp@compareClusterResult$DB == "HALLMARK" | comp@compareClusterResult$DB == "BIOCARTA"), ]
comp.filt <- comp
comp.filt@compareClusterResult <- filt
glimpse(comp.filt@compareClusterResult)

filt <- comp.filt[-grep("HIV", comp.filt@compareClusterResult$ID), ]
comp.filt2 <- comp.filt
comp.filt2@compareClusterResult <- filt
glimpse(comp.filt2@compareClusterResult)

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/01_msigdb.C2.H_filt_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp.filt2,
  x = "Count",
  color = "p.adjust",
  showCategory = 25,
  label_format = 50
) +
  theme_bw() +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, colour = "black")
  )
dev.off()

################################################################################################################################################################
# FIGURE 2 B

to.plot <- readRDS(
  "Output/Rdata/03_states.1_2024.05.02/02.1_to.plot_MGH106_2024.05.02.RDS"
)
clust.rows <- readRDS(
  "Output/Rdata/03_states.1_2024.05.02/02.5_clust.rows_MGH106_2024.03.26.RDS"
)

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/02_MGH106.Heatmap.pdf", height = 8, width = 9)
pheatmap(
  to.plot,
  color = colorRampPalette(paletteer_d(palette = "RColorBrewer::PuOr", 11))(100),
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  use_raster = FALSE
)
dev.off()

################################################################################
# FIGURE 2 C

pivot <- readRDS("Output/Rdata/04_states.2_2024.05.02/01_pivot.mtx_2024.04.05.RDS")

to.plot <- readRDS("Output/Rdata/04_states.2_2024.05.02/02_to.plot_2024.04.05.RDS")


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
  tumor.ID = generate.colors_tumors(ann)[[1]],
  # state = generate.colors(ann.row)[[2]],
  # phase = generate.colors(ann.row)[[3]],
  Cluster = generate.colors_states(ann.col)[[1]]
)
ann.colors.col$Cluster <- ann.colors.col$Cluster[-4]

scaled <- scale(to.plot, center = TRUE, scale = FALSE)
my.breaks <- c(seq(min(to.plot), 0, length.out=ceiling(100/2) + 1),
               seq(max(to.plot)/100, max(to.plot), length.out=floor(100/2)))

diag(to.plot) <- NA
# generate heatmap object
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/02_Vulnerability.States.Heatmap.pdf", height = 8, width = 9)
pheatmap(
  to.plot,
  color = colorRampPalette(paletteer_c(palette = "viridis::magma", 100, direction = -1))(100),
  annotation_row = ann,
  annotation_col = ann.col,
  annotation_colors = ann.colors.col,
  cluster_rows = clust.rows,
  cluster_cols = clust.rows,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE,
  # breaks = my.breaks
)
dev.off()

################################################################################################################################################################
# FIGURE 3 A

sigs <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)

signatures <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/04_VS.signatures.all_2024.04.05.RDS"
)

nef.np <- readRDS(
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed.NP_VC.ANN_2024.05.06.RDS"
)

# make exhaustive list of all top 100 marker genes per state, including ----
# duplicated genes per signature
mark.cor.list <- list()
for (i in 1:length(signatures)) {
  sub <- signatures[[i]][which(signatures[[i]]$gene %in% sigs[[i]]$gene), ]
  sub2 <- sub[, -c(4, 6:10)]
  mark.cor.list[[i]] <- sub2
}
glimpse(mark.cor.list)

genes <- read.csv( # read in 300-gene list
  "Output/Rdata/04_states.2_2024.05.02/02_all.VS.markers_noANN_2024.04.05.csv",
  row.names = 1
)
genes <- genes$x
genes <- unique(genes)

# add average expression of genes in cells w/in k_clust ----

files <- list.files(
  path = "Output/Rdata/03_states.1_2024.05.02",
  pattern = "^03_cluster.data"
)
path = "Output/Rdata/03_states.1_2024.05.02/"

# start the loop for extraction of expression data of genes of interest
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

    sub.obj <- subset(nef.np, cells = cells)

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
  tumor.ID = generate.colors_tumors(ann2)[[1]],
  state = generate.colors_states(ann2)[[2]]
)
# VS out of order for some reason and can't figure out how to rearrange
ann.colors.col$state <- ann.colors.col$state[-4]
num2 <- ann.colors.col$state[3]
ann.colors.col$state[3] <- ann.colors.col$state[2]
ann.colors.col$state[2] <- num2

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
smol2$state <- as.factor(smol2$state)

scaled <- ScaleData(mtx)

my.breaks <- c(seq(min(scaled), 0, length.out=ceiling(100/2) + 1),
               seq(max(scaled)/100, max(scaled), length.out=floor(100/2)))

# need to order genes appropriately (order by signature)
markers <- sigs %>%
  purrr::reduce(full_join, by = "gene")
hmm <- scaled[markers$gene, ] # markers should be ordered by log2FC

labels <- c(
  "MYT1L", "EPHB6", "SCN3", "ANK3", "PTPRN", "CLDN18",
  "PLP1", "MYRF", "CLDN", "RPS5", "RPS13", "RPS29", "RPL27", "RPL34",
  "AURKA", "AURKB", "CDC20", "CDK1", "EZH2", "BARD1"
)

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

ph <- pheatmap::pheatmap(
  as.matrix(hmm),
  color = colorRampPalette(paletteer_d(palette = "RColorBrewer::RdBu", 11, direction = -1))(100),
  annotation_row = smol2,
  annotation_colors = ann.colors.col[2],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  use_raster = FALSE,
  breaks = my.breaks
)
dev.off()

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_VS.Markers.Heatmap_GENES_2023.09.25.pdf", height = 12, width = 8)
add.flag(ph,
         kept.labels = labels,
         repel.degree = 0)
dev.off()

################################################################################
# FIGURE 3 B

markers.list <- list(
  VS1 = sigs[[1]]$gene,
  VS2 = sigs[[2]]$gene,
  VS3 = sigs[[3]]$gene
)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb2 <- msigdbr(species = "Homo sapiens", category = "H")
msigdb6 <- msigdbr(species = "Homo sapiens", category = "C3")
# plot Enrichment score in HM and arrange w/ clustering?

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
table(comp@compareClusterResult$Cluster)

# remove cancer-irrelevant enrichment results from top 5 / cluster
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

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_msigdb.C2.C3.H_filt_ORA_2024.04.24.pdf", height = 8, width = 10)
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
# FIGURE 3 C (GDS x DEG)

to.hm <- read.csv (
  "Output/Rdata/05_analysis_2024.05.03/07_VC.top.1500_2024.05.06.csv",
  row.names = 1, header = TRUE
)

Idents(nef) <- nef$mVC
avg.clust <- AverageExpression(nef, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(to.hm$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% to.hm$gene), ]
# sort expression by gene

gene.labels <- c(
  rownames(filt_avg.clust.def)[which(rownames(filt_avg.clust.def) %in% killing$Gene)]
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
  color = colorRampPalette(paletteer_d(palette = "RColorBrewer::RdBu", 11, direction = -1))(100),
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

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_1500.DEG.pHeatmap_2024.05.05.pdf", height = 6, width = 8)
add.flag(ph,
         kept.labels = gene.labels,
         repel.degree = 0)
dev.off()

################################################################################
# FIGURE 3 D, E, F (VLN)

# neftel ----

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_VS.MS_Neftel_Vln_2023.10.18.pdf", height = 4, width = 8)
VlnPlot(
  nef,
  features = "ess.MS_1",
  group.by = "mVC", pt.size = 0.01,
  cols = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(nef$mVC)))
) +
  geom_boxplot() +
  ylab("Module Score \n") +
  ylim(range(nef$ess.MS_1)[1], range(nef$ess.MS_1)[2]+0.2) +
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
  MS = nef$ess.MS_1,
  VS = nef$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

test.res <- aov(MS ~ VS, data = test.df)
print(test.res)

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

# johnson ----

john <- readRDS(
  "data/johnson_2024.10.04.RDS"
)

test.df <- data.frame(
  MS = john$ess.MS_1,
  State = john$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_VS.MS_Johnson_Vln_2023.10.18.pdf", height = 4, width = 8)
VlnPlot(
  john,
  features = "ess.MS_1",
  group.by = "mVC", pt.size = 0.01,
  cols = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(nef$mVC)))
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

# suter ----

suter <- readRDS(
  "data/Suter_scGBM_2023.10.04.RDS"
)

test.df <- data.frame(
  MS = suter$ess.MS_1,
  State = suter$mVC
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/03_VS.MS_Suter_Vln_2023.10.18.pdf", height = 4, width = 8)
VlnPlot(
  suter,
  features = "ess.MS_1",
  group.by = "mVC", pt.size = 0.01,
  cols = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(nef$mVC)))
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


################################################################################################################################################################
# FIGURE 4 BARPLOTS

# neftel ----

nef.np <- readRDS(
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed.NP_VC.ANN_2024.05.06.RDS"
)

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/04_bar_VS.prop_2024.05.09.pdf", height = 4, width = 8)
dittoBarPlot(
  nef.np,
  var = "mVC",
  group.by = "tumor.id",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Cells \n",
  main = NULL,
  color.panel = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(nef.np$mVC))+1)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

# johnson ----

john.np <- readRDS(
  "data/johnson.NP_2024.10.04.RDS"
)

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/04_JOHN_bar_VS.prop_2024.05.09.pdf", height = 4, width = 8)
dittoBarPlot(
  john.np,
  var = "mVC",
  group.by = "patient.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(john.np$mVC))+1)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

# suter ----

sut.np <- readRDS(
  "data/Suter_scGBM.NP_2023.10.04.RDS"
)

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/04_SUTER_bar_VS.prop_2024.05.09.pdf", height = 4, width = 8)
dittoBarPlot(
  sut.np,
  var = "mVC",
  group.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(unique(sut.np$mVC))+1)
) +
  mrd.theme +
  theme(legend.position = "right")
dev.off()


################################################################################################################################################################
# FIGURE 5, SPATIAL FIGURES

# dir.create("Output/Figures/15_spatial_2024.05.09")

spatial <- readRDS("Output/Rdata/15_spatial_2024.05.09/08_deconv_v5_2024.05.09.RDS")

props.state <- data.frame(
  row.names = rownames(spatial@meta.data),
  VS1 = spatial$propVS1,
  VS2 = spatial$propVS2,
  VS3 = spatial$propVS3,
  Non = spatial$propNON
)

img <- Images(spatial)
run <- img[grep("_C_", img, invert = TRUE)] # remove entry cortex slices
Idents(spatial) <- spatial$orig.ident

dir3 <- "Output/Figures/XX_MAIN.FIGS.enhanced/05_spatial_he"
dir.create(dir3)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(spatial, idents = slice)
  pdf(paste0(dir3, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = spatial,
      y = props.state[spots, ],
      cell_types = colnames(props.state),
      img = TRUE, slice = run[i],
      pie_scale = 0.4, scatterpie_alpha = 0.9
    ) +
      scale_fill_manual(
        values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(colnames(props.state))),
        breaks = colnames(props.state)
      )
  )
  dev.off()
}

dir3.1 <- "Output/Figures/XX_MAIN.FIGS.enhanced/05_spatial_emty"
dir.create(dir3.1)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(spatial, idents = slice)
  pdf(paste0(dir3.1, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = spatial,
      y = props.state[spots, ],
      cell_types = colnames(props.state),
      img = FALSE, slice = run[i], ####
      pie_scale = 0.40, scatterpie_alpha = 1.0
    ) +
      scale_fill_manual(
        values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(colnames(props.state))),
        breaks = colnames(props.state)
      )
  )
  dev.off()
}

# 243, 260, 269 : less spaces between spots (increase pie_scale)

slice <- "UKF243_T_ST"
spots <- WhichCells(spatial, idents = slice)
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/05_UKF243_spatial.emty.pie_2024.05.09.pdf")
plotSpatialScatterpie(
  x = spatial,
  y = props.state[spots, ],
  cell_types = colnames(props.state),
  img = FALSE, slice = slice, ####
  pie_scale = 0.45, scatterpie_alpha = 1.0
) +
  scale_fill_manual(
    values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(colnames(props.state))),
    breaks = colnames(props.state)
  )
dev.off()

slice <- "UKF260_T_ST"
spots <- WhichCells(spatial, idents = slice)
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/05_UKF260_spatial.emty.pie_2024.05.09.pdf")
plotSpatialScatterpie(
  x = spatial,
  y = props.state[spots, ],
  cell_types = colnames(props.state),
  img = FALSE, slice = slice, ####
  pie_scale = 0.47, scatterpie_alpha = 1.0
) +
  scale_fill_manual(
    values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(colnames(props.state))),
    breaks = colnames(props.state)
  )
dev.off()

slice <- "UKF269_T_ST"
spots <- WhichCells(spatial, idents = slice)
pdf("Output/Figures/XX_MAIN.FIGS.enhanced/05_UKF269_spatial.emty.pie_2024.05.09.pdf")
plotSpatialScatterpie(
  x = spatial,
  y = props.state[spots, ],
  cell_types = colnames(props.state),
  img = FALSE, slice = slice, ####
  pie_scale = 0.45, scatterpie_alpha = 1.0
) +
  scale_fill_manual(
    values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(colnames(props.state))),
    breaks = colnames(props.state)
  )
dev.off()

################################################################################
# FIGURE 5, BARPLOT

quant <- data.frame(
  row.names = rownames(spatial@meta.data),
  sliceID = spatial@meta.data$orig.ident,
  VS1 = spatial@meta.data$propVS1,
  VS2 = spatial@meta.data$propVS2,
  VS3 = spatial@meta.data$propVS3,
  Astrocytes = spatial@meta.data$propAstro,
  Endothelial = spatial@meta.data$propEndo,
  Fibroblast = spatial@meta.data$propFibro,
  Myeloid = spatial@meta.data$propMyelo,
  Oligodendrocyte = spatial@meta.data$propOligo,
  T.Cell = spatial@meta.data$propT
)

# remove "C" slices ("NA" values in proportions)
quant <- quant[complete.cases(quant), ]
glimpse(quant)

quant.long <- quant %>%
  pivot_longer(cols = -sliceID, names_to = "CellType", values_to = "Proportion")
glimpse(quant.long)

quant.ag <- quant.long %>%
  aggregate(Proportion ~ sliceID + CellType, FUN = mean)
glimpse(quant.ag)
quant.ag$CellType <- factor(quant.ag$CellType, levels = c("VS1", "VS2", "VS3", "Astrocytes", "Endothelial", "Fibroblast", "Myeloid", "Oligodendrocyte", "T.Cell"))

#### can plot all cell types here

non <- c("Astrocytes", "Endothelial", "Fibroblast", "Myeloid", "Oligodendrocyte", "T.Cell")
quant.ag <- quant.ag[-which(quant.ag$CellType %in% non), ]

quant.np <- quant.ag %>%
  pivot_wider(names_from = CellType, values_from = Proportion)
quant.np <- quant.np %>%
  mutate(
    sum = VS1+VS2+VS3
  )
quant.np <- quant.np %>%
  mutate(
    fct = 1/sum
  )
quant.np <- quant.np %>%
  mutate(
    newVS1 = VS1*fct,
    newVS2 = VS2*fct,
    newVS3 = VS3*fct
  )
glimpse(quant.np)
quant.np <- quant.np[, -c(2:6)]
colnames(quant.np) <- c("sliceID", "VS1", "VS2", "VS3")

quant.long <- quant.np %>%
  pivot_longer(cols = -sliceID, names_to = "CellType", values_to = "Proportion")
glimpse(quant.long)

quant.ag <- quant.long %>%
  aggregate(Proportion ~ sliceID + CellType, FUN = mean)
glimpse(quant.ag)
quant.ag$CellType <- factor(quant.ag$CellType, levels = c("VS1", "VS2", "VS3"))

pdf("Output/Figures/XX_MAIN.FIGS.enhanced/05_spatial.bar.VS_2024.05.09.pdf", height = 4, width = 8)
quant.ag %>%
  ggplot(aes(x = sliceID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(
    values = colorRampPalette(paletteer_d(palette = "RColorBrewer::Dark2", 8))(length(levels(quant.ag$CellType))+1),
    breaks = levels(quant.ag$CellType)
  ) +
  theme_minimal() +
  mrd.theme +
  theme(legend.position = "right") +
  RotatedAxis()
dev.off()
