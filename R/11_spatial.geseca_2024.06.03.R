library(tidyverse)
library(Seurat)
library(cowplot)
library(msigdbr)
library(fgsea)
library(corrplot)
library(paletteer)
library(ggpubr)
library(data.table)

dir.create("Output/Figures/11_spatial/06_geseca")
dir.create("Output/Rdata/11_spatial/06_geseca")
dir.create("Output/Rdata/11_spatial/06_geseca/objects")

################################################################################
# READ IN DATA

obj <- readRDS("Output/Rdata/11_spatial/08_deconv_v5_2024.05.09.RDS")
glimpse(obj@meta.data)
glimpse(obj@assays)

# select images to run
img <- Images(obj)
run <- img[grep("_C_", img, invert = TRUE)] # remove entry cortex slices

# read in VS signatures
sigs <- readRDS(
  "Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS"
)
vs <- list()
for (i in 1:length(sigs)) {
  genes <- sigs[[i]]$gene
  vs[[i]] <- genes
}
names(vs) <- names(sigs)
glimpse(vs)

################################################################################
# SET UP ENRICHMENT

# geseca function
score.geseca <- function(object, pathways, assay = DefaultAssay(object), prefix = "", scale = FALSE) {
  dat <- GetAssay(object, assay)
  expr <- dat$scale.data
  res <- object
  for (i in seq_along(pathways)) {
    genes <- pathways[[i]]
    overlap <- intersect(unique(genes), rownames(expr))
    score <- colSums(expr[overlap, , drop = FALSE])/sqrt(length(overlap))
    scaled <- scale(score, center = TRUE, scale = scale)
    res <- AddMetaData(res, metadata = score, col.name = paste0("raw_", names(pathways)[i]))
    res <- AddMetaData(res, metadata = scaled, col.name = paste0(prefix, names(pathways)[i]))
    # res@meta.data[[paste0(prefix, names(pathways)[i])]] <- score # BAD
  }
  return(res)
}

# score VS function
score.vs <- function(object, set, assay = DefaultAssay(object), prefix = "sing", scale = FALSE) {
  dat <- GetAssay(object, assay)
  expr <- as.data.frame(dat$scale.data)
  rank.data <- singscore::rankGenes(expr)
  res <- object
  for (i in seq_along(set)) {
    genes <- set[[i]]
    state <- paste0(names(set)[i])
    new <- GeneSet()
    new@geneIds <- as.character(genes)
    scored <- singscore::simpleScore(rankData = rank.data, upSet = new)
    res <- AddMetaData(res, metadata = scored[, "TotalScore"], col.name = paste0(prefix, "_", names(set)[i]))
  }
  return(res)
}

# select pathways for enrichment
# pathwaysDF <- msigdbr("human", category = "H")
pathwaysDF.2 <- msigdbr("human", category = "C5")
# pathwaysdf <- rbind(pathwaysDF, pathwaysDF.2)
pathwaysdf <- pathwaysDF.2[grep("^GOBP_", pathwaysDF.2$gs_name), ] # GOBP enrichment
pathways <- split(pathwaysdf$gene_symbol, pathwaysdf$gs_name)
# pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

################################################################################
# START LOOP

# scaling already done w/ regressing out sct across all images in dataset,
# no need to rescale w/in loop

for (i in 1:length(run)) {

  inx <- which(Images(obj) == run[i]) # get image indices to keep

  image.name <- run[i]
  image.metadata <- obj@images[[inx]]
  spot.ids <- rownames(image.metadata@coordinates)

  subset.obj <- subset(obj, cells = spot.ids)
  keys.to.remove <- setdiff(Images(subset.obj), run[i])
  for (key in keys.to.remove) {
    subset.obj@images[[key]] <- NULL
  }

  DefaultAssay(subset.obj) <- "Spatial"
  subset.obj <- RunPCA(
    subset.obj,
    assay = "Spatial",
    rev.pca = TRUE,
    reduction.name = "pca.rev",
    reduction.key = "PCr_", npcs = 50,
    verbose = FALSE
  )

  E <- subset.obj@reductions$pca.rev@feature.loadings

  # enrichment analysis
  set.seed(707)
  results <- geseca(
    pathways, E,
    minSize = 2, # ferroptosis gene set only has 5 genes in it
    maxSize = 500,
    center = FALSE
  )
  saveRDS(results, paste0("Output/Rdata/08_geseca_no.scale/", i, "_", image.name, "_geseca.results_2024.06.03.RDS")) # supplemental file

  results <- results[which(results$padj < 0.05), ] # filter for significant terms
  topPathways <- results[, pathway] |> head(50) # change number for top n

  # geseca score
  scored.obj <- score.geseca(subset.obj, pathways[topPathways], assay = "Spatial", prefix = "", scale = FALSE)
  scored.obj <- score.vs(scored.obj, vs, assay = "Spatial") # singscore vs
  scored.obj <- score.vs(scored.obj, pathways[topPathways], assay = "Spatial") # singscore pathways
  saveRDS(scored.obj, paste0("Output/Rdata/08_geseca_no.scale/objects/", image.name, "_geseca.scored.obj_2024.06.03.RDS"))

  # enrichment plots
  for (j in seq_along(pathways[topPathways])) {
    pth <- names(pathways[topPathways])[j]
    plot <- SpatialFeaturePlot(
      scored.obj,
      features = names(pathways[topPathways])[j],
      combine = FALSE,
      image.alpha = 0, pt.size.factor = 2.5, stroke = 0.8,
    )[[1]]
    plot$scales$scales[plot$scales$find("fill")] <- NULL
    suppressMessages({
      p2 <- plot +
        scale_fill_gradientn(
          limits = c(-4, 4), breaks = c(-4, 0, 4),
          oob = scales::squish,
          colors = c("darkblue", "lightgrey", "darkred"),
          guide = "colourbar", name = "z-score"
        ) +
        theme(legend.position = theme_get()$legend.position)
    })

    dir = paste0("Output/Figures/08_geseca_no.scale/", image.name)
    dir.create(dir)
    pdf(paste0(dir, "/", i, "_", image.name, "_", pth, "_H.C5_2024.05.30.pdf"))
    print(p2)
    dev.off()

  }

  # correlation analysis

  cor.df <- data.frame(
    row.names = rownames(scored.obj@meta.data),
    VS1 = scored.obj$sing_VS1,
    VS2 = scored.obj$sing_VS2,
    VS3 = scored.obj$sing_VS3,
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[1]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[2]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[3]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[4]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[5]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[6]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[7]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[8]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[9]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[10]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[11]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[12]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[13]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[14]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[15]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[16]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[17]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[18]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[19]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[20]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[21]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[22]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[23]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[24]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[25]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[26]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[27]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[28]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[29]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[30]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[31]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[32]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[33]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[34]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[35]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[36]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[37]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[38]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[39]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[40]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[41]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[42]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[43]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[44]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[45]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[46]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[47]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[48]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[49]))],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == paste0("sing_", names(pathways[topPathways])[50]))]
  )
  colnames(cor.df)[4:length(colnames(cor.df))] <- c(names(pathways[topPathways])[1:length(names(pathways[topPathways]))])

  mtx <- as.matrix(cor.df)
  # saveRDS(mtx, paste0("Output/Rdata/08_geseca_no.scale/01_", image.name, "_corr.df_2024.05.09.RDS"))

  # make correlation
  cor.mtx <- cor(
    mtx,
    use = "everything",
    method = "spearman"
  )
  saveRDS(
    cor.mtx,
    paste0("Output/Rdata/08_geseca_no.scale/", i, "_", image.name, "_corr.mtx_2024.05.09.RDS")
  )

  # corrplot

  cor.mtx[is.na(cor.mtx)] <- 0 # set for any NA values (if no VS in one section)
  pdf(paste0("Output/Figures/08_geseca_no.scale/", image.name, "/01_sigs.corrplot_2024.05.09.pdf"), height = 8, width = 8)
  print(
    corrplot::corrplot(
      cor.mtx,
      tl.pos = "l", tl.col = "black", #tl.cex = 0.2,
      order = "hclust",
      hclust.method = "ward.D2",
      method = "color",
      type = "full", # addrect = 3, # rect.col = "red",
      col = colorRampPalette(paletteer_d(palette = "RColorBrewer::BrBG", 11))(100),
      is.corr = TRUE,
      diag = TRUE, cl.pos = "b",
      na.label = "square", na.label.col = "grey"
    )
  )
  dev.off()

}

################################################################################
# EXAMINE TOP ENRICHED PATHWAYS PER VULNERABILITY STATE

files <- list.files(
  "Output/Rdata/11_spatial/06_geseca/objects"
)

meta.list <- list()
dir = "Output/Rdata/11_spatial/06_geseca/objects/"
for (f in seq_along(files)) {
  tmp <- readRDS(paste0(dir, files[f]))
  meta <- tmp@meta.data

  meta.list[[f]] <- meta
}

new <- rbindlist(meta.list, fill = TRUE)
glimpse(new)

# only use spots with < 10% non-neoplastic content
filt <- new[which(new$propNON < 0.2), ]

# filter only for columns with GOBP singscores
terms <- colnames(filt)[grep("^sing_GOBP_", colnames(filt))]

# VS1 correlations
cor.res1 <- data.frame()
for (e in seq_along(terms)) {
  # n <- filt[which(filt$propVS1 > 0.5), ] # more permissive with VS1 spots bc so few
  x <- as.data.frame(filt)[, which(colnames(filt) == terms[e])]
  c <- cor(x = x, y = filt$sing_VS1, use = "complete.obs")
  add <- c(terms[e], c)
  cor.res1 <- rbind(cor.res1, add)
}
colnames(cor.res1) <- c("Term", "Corr")
cor.res1$Corr <- as.double(cor.res1$Corr)
range(cor.res1$Corr, na.rm = TRUE)
vs1.keep <- cor.res1[which(cor.res1$Corr > 0.3), ]
vs1.paths <- vs1.keep$Term
vs1.paths

# VS2 correlations
cor.res2 <- data.frame()
for (e in seq_along(terms)) {
  n <- filt[which(filt$propVS2 > 0.5), ] # only use spots with predicted VS prop > 50%
  x <- as.data.frame(n)[, which(colnames(n) == terms[e])]
  if (length(x[complete.cases(x)]) > 0) {
    c <- cor(x = x, y = n$sing_VS2, use = "complete.obs")
    add <- c(terms[e], c)
    cor.res2 <- rbind(cor.res2, add)
  }
  else { print(paste(terms[e], "does not have any singscores for VS2")) }
}
colnames(cor.res2) <- c("Term", "Corr")
cor.res2$Corr <- as.double(cor.res2$Corr)
range(cor.res2$Corr, na.rm = TRUE)
vs2.keep <- cor.res2[which(cor.res2$Corr > 0.3), ]
vs2.paths <- vs2.keep$Term
vs2.paths

# VS3 correlations
cor.res3 <- data.frame()
for (e in seq_along(terms)) {
  n <- filt[which(filt$propVS3 > 0.5), ]
  x <- as.data.frame(n)[, which(colnames(n) == terms[e])]
  if (length(x[complete.cases(x)]) > 0) {
    c <- cor(x = x, y = n$sing_VS3, use = "complete.obs")
    add <- c(terms[e], c)
    cor.res3 <- rbind(cor.res3, add)
  }
  else { print(paste(terms[e], "does not have any singscores for VS3")) }
}
colnames(cor.res3) <- c("Term", "Corr")
cor.res3$Corr <- as.double(cor.res3$Corr)
range(cor.res3$Corr, na.rm = TRUE)
vs3.keep <- cor.res3[which(cor.res3$Corr > 0.3), ]
vs3.paths <- vs3.keep$Term
vs3.paths




