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

mrd.theme <- readRDS("Output/mrd.fig.theme.RDS")
umap.theme <- readRDS("Output/mrd.umap.theme.RDS")

################################################################################

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
# FIGURE 1 A, B, C

efx <- read.delim(
  file = "data/r.input_CRISPRGeneEffect_23Q2.csv",
  sep = ",",
  row.names = 1,
  check.names = FALSE
)
colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))

write.table(
  gbm.stats_m,
  file = "Output/Rdata/01_essential.genes/02_depmap.gbm.stats_mean.effect.scores_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# effect scores for all genes in all cell lines
all <- t(as.data.frame(lapply(efx, median)))
all <- as.data.frame(all)
colnames(all) <- c("mEffectScore")

pdf("Output/Figures/01_essential.genes/density/00_median.effect.score_all.genes_all.lines_2024.02.07.pdf", height = 4, width = 4)
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

all <- t(as.data.frame(lapply(efx, mean)))
all <- as.data.frame(all)
colnames(all) <- c("mEffectScore")

pdf("Output/Figures/01_essential.genes/density/00_mean.effect.score_all.genes_all.lines_2024.02.07.pdf", height = 4, width = 4)
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
pdf("Output/Figures/01_essential.genes/density/01_mean.effect.score_all.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
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
pdf("Output/Figures/01_essential.genes/density/02_mean.effect.score_sig.lower.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
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
pdf("Output/Figures/01_essential.genes/density/03_mean.effect.score_killing.genes_GBM.lines_2024.02.07.pdf", height = 4, width = 4)
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

saveRDS(
  obj,
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
  MS = obj$ess.MS_1,
  CellTypes = obj$neoplastic.state
)
hist(test.df$MS) # MS roughly normally distributed
dev.off()

pdf("Output/Figures/02_scRNA.neftel/03_essential.genes_module.score_npVnon_2024.02.05.pdf", height = 8, width = 4)
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

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

# filter only for databases of interest
filt <- comp[which(comp@compareClusterResult$DB == "REACTOME" | comp@compareClusterResult$DB == "KEGG" | comp@compareClusterResult$DB == "HALLMARK" | comp@compareClusterResult$DB == "BIOCARTA"), ]
comp.filt <- comp
comp.filt@compareClusterResult <- filt

glimpse(comp.filt@compareClusterResult)
View(comp.filt@compareClusterResult)

filt <- comp.filt[-grep("HIV", comp.filt@compareClusterResult$ID), ]
comp.filt2 <- comp.filt
comp.filt2@compareClusterResult <- filt
glimpse(comp.filt2@compareClusterResult)

pdf("Output/Figures/01_essential.genes/pan_enrich/05_msigdb.C2.H_filt_ORA_2024.03.13.pdf", height = 8, width = 10)
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

saveRDS(
  to.plot,
  paste0(dir, "02.1_to.plot_", names(input.list)[i], "_2024.05.02.RDS")
)

pdf(paste0(dir2, "01_", names(input.list)[i], "_clustered.pheatmap_2024.03.27.pdf"), height = 8, width = 9)
print(
  pheatmap(
    to.plot,
    color = viridis_pal(option = "inferno")(100),
    cluster_rows = clust.rows,
    cluster_cols = clust.rows,
    show_rownames = FALSE,
    show_colnames = FALSE,
    use_raster = FALSE
  )
)
dev.off()

################################################################################
# FIGURE 2 C

saveRDS(pivot, "Output/Rdata/04_states.2_2024.05.02/01_pivot.mtx_2024.04.05.RDS")

saveRDS(to.plot, "Output/Rdata/04_states.2_2024.05.02/02_to.plot_2024.04.05.RDS")


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

################################################################################################################################################################
# FIGURE 3 A

sigs <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)

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
# names(ann.colors.col$state) <- c("1", "2", "3")
# ann.colors.col$state[[1]] <- "#30123BFF"
# ann.colors.col$state[[2]] <- "#A2FC3CFF"
# ann.colors.col$state[[3]] <- "#7A0403FF"

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
# hmm <- hmm[rownames(smol2), ] # smol is ordered by tumor


pdf(paste0("Output/Figures/05_analysis_2024.05.03/23_VS.Markers.Heatmap_sorted.noclust_2023.09.25.pdf"), height = 12, width = 8)
pheatmap::pheatmap(
  as.matrix(hmm),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)), #paletteer_c("ggthemes::Red-Blue-White Diverging", 100, direction = -1)
  annotation_row = smol2,
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
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)), #paletteer_c("ggthemes::Red-Blue-White Diverging", 100, direction = -1)
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

pdf("Output/Figures/05_analysis_2024.05.03/24_VS.Markers.Heatmap_GENES_2023.09.25.pdf", height = 12, width = 8)
add.flag(ph,
         kept.labels = labels,
         repel.degree = 0)
dev.off()

################################################################################
# FIGURE 3 B

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

################################################################################
# FIGURE 3 C (GDS x DEG)

to.hm <- read.csv (
  "Output/Rdata/05_analysis_2024.05.03/07_VC.top.1500_2024.05.06.csv",
  row.names = 1, header = TRUE
)

Idents(obj) <- obj$mVC
avg.clust <- AverageExpression(obj, group.by = "ident", assay = "RNA", layer = "scale.data")

order_vec <- match(unique(to.hm$gene), rownames(avg.clust$RNA))
sorted_avg.clust.def <- avg.clust$RNA[order_vec, ]

filt_avg.clust.def <- sorted_avg.clust.def[which(rownames(sorted_avg.clust.def) %in% to.hm$gene), ]
# sort expression by gene

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
  color = paletteer_c("ggthemes::Red-Blue-White Diverging", 100, direction = -1),
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
# FIGURE 3 D, E, F (VLN)

# neftel ----

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

# johnson ----

# suter ----

saveRDS(
  suter,
  "data/Suter_scGBM_2023.10.04.RDS"
)

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


################################################################################################################################################################
# FIGURE 4 BARPLOTS

# neftel ----

saveRDS(
  obj.np,
  "Output/Rdata/05_analysis_2024.05.03/00_neftel.noPed.NP_VC.ANN_2024.05.06.RDS"
)

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

# johnson ----


# suter ----

saveRDS(
  obj.np,
  "data/Suter_scGBM.NP_2023.10.04.RDS"
)

pdf("Output/Figures/07_suter_2024.05.06/04_Essential.Clusters_Barplot_Tumor.NonNP_2023.12.08.pdf", height = 4, width = 8)
dittoBarPlot(
  obj.np,
  var = "mVC",
  group.by = "orig.ident",
  scale = "percent",
  y.breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  xlab = NULL,
  ylab = "Percent of Neoplastic Cells \n",
  main = NULL,
  color.panel = viridis_pal(option = "turbo")(length(unique(obj.np$mVC)))
) +
  mrd.theme
dev.off()


################################################################################
