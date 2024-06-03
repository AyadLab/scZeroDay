library(tidyverse)
library(Seurat)
library(cowplot)
library(msigdbr)
library(fgsea)
library(corrplot)
library(paletteer)

dir.create("Output/Figures/06_geseca")
dir.create("Output/Rdata/06_geseca")
dir.create("Output/Rdata/06_geseca/objects")

################################################################################
# READ IN DATA

obj <- readRDS("Output/Rdata/11_spatial/08_deconv_v5_2024.05.09.RDS")
glimpse(obj@meta.data)
glimpse(obj@assays)

# non-regressed scale.data
non.regr <- readRDS("Output/Rdata/06_scale.data.not.regressed_2024.04.05.RDS")
obj@assays$Spatial$scale.data <- non.regr

# select images to run
img <- Images(obj)
run <- img[grep("_C_", img, invert = TRUE)] # remove entry cortex slices

################################################################################
# SET UP ENRICHMENT

pathwaysDF <- msigdbr("human", category = "H")
pathwaysDF.2 <- msigdbr("human", category = "C5")
pathwaysdf <- rbind(pathwaysDF, pathwaysDF.2)
pathwaysdf <- pathwaysdf[-grep("^HP_", pathwaysdf$gs_name), ]
pathways <- split(pathwaysdf$gene_symbol, pathwaysdf$gs_name)
# pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

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

################################################################################
# START LOOP

for (i in 1:length(run)) {

  inx <- which(Images(obj) == run[i]) # get image indices to keep

  image.name <- run[i]
  image.metadata <- obj@images[[inx]]
  spot.ids <- rownames(image.metadata@coordinates)

  # subset seurat object / image
  subset.obj <- subset(obj, cells = spot.ids)
  keys.to.remove <- setdiff(Images(subset.obj), run[i])
  for (key in keys.to.remove) {
    subset.obj@images[[key]] <- NULL
  }

  # run reverse pca (on genes vs. on cells)
  subset.obj <- RunPCA(
    subset.obj,
    assay = "Spatial", # SCT
    rev.pca = TRUE,
    reduction.name = "pca.rev",
    reduction.key = "PCr_", npcs = 50,
    verbose = FALSE
  )

  # isolate pc results
  E <- subset.obj@reductions$pca.rev@feature.loadings

  # enrichment analysis
  set.seed(707)
  results <- geseca(
    pathways, E,
    minSize = 15,
    maxSize = 500,
    center = FALSE
  )
  saveRDS(results, paste0("Output/Rdata/06_geseca/", i, "_", image.name, "_geseca.results_2024.06.03.RDS"))

  # top pathways
  topPathways <- results[, pathway] |> head(10) # change number for top n

  # geseca score
  scored.obj <- score.geseca(subset.obj, pathways[topPathways], assay = "Spatial", prefix = "", scale = FALSE)
  saveRDS(scored.obj, paste0("Output/Rdata/06_geseca/objects/", image.name, "_geseca.scored.obj_2024.06.03.RDS"))

  # enrichment plots
  for (j in seq_along(pathways[topPathways])) {
    pth <- names(pathways[topPathways])[j]
    plot <- SpatialFeaturePlot(
      scored.obj,
      features = names(pathways[topPathways])[j],
      combine = FALSE,
      image.alpha = 0, pt.size.factor = 2.5, stroke = 0.8
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

    dir = paste0("Output/Figures/06_geseca/", image.name)
    dir.create(dir)
    pdf(paste0(dir, "/", i, "_", image.name, "_", pth, "_H.C5_2024.05.30.pdf"))
    print(p2)
    dev.off()

  }

  # correlation analysis: vs proportions x geseca enrichment z-score (spearman)

  cor.df <- data.frame(
    row.names = rownames(scored.obj@meta.data),
    VS1 = scored.obj$propVS1,
    VS2 = scored.obj$propVS2,
    VS3 = scored.obj$propVS3,
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[1])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[2])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[3])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[4])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[5])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[6])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[7])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[8])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[9])],
    scored.obj@meta.data[, which(colnames(scored.obj@meta.data) == names(pathways[topPathways])[10])]
  )
  colnames(cor.df)[4:13] <- c(names(pathways[topPathways])[1:10])

  mtx <- as.matrix(cor.df)
  # saveRDS(mtx, paste0("Output/Rdata/06_geseca/01_", image.name, "_corr.df_2024.05.09.RDS"))

  # make correlation
  cor.mtx <- cor(
    mtx,
    use = "everything",
    method = "spearman"
  )
  saveRDS(
    cor.mtx,
    paste0("Output/Rdata/06_geseca/", i, "_", image.name, "_corr.mtx_2024.05.09.RDS")
  )

  # corrplot

  cor.mtx[is.na(cor.mtx)] <- 0 # set for any NA values (if no VS in one section)
  pdf(paste0("Output/Figures/06_geseca/", image.name, "/01_sigs.corrplot_2024.05.09.pdf"), height = 8, width = 8)
  print(
    corrplot::corrplot(
      cor.mtx,
      tl.pos = "l", tl.col = "black", tl.cex = 0.2,
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





