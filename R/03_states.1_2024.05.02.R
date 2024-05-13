library(tidyverse)
library(Seurat)
library(viridis)
library(ggplot2)
library(viridis)
library(pheatmap)
library(corrplot)
library(proxy)
library(factoextra)
library(FactoMineR)
library(NbClust)
library(cowplot)
library(cluster)
library(ComplexHeatmap)
library(clustree)

setwd("/data/mrd/cognitive.seeds_zero")
dir.create("Output/Rdata/03_states.1_2024.05.02")
dir.create("Output/Figures/03_states.1_2024.05.02")

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
# READ DATA

# dependency gene list
ess <- read.csv(
  "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv"
)

obj <- readRDS("Output/Rdata/02_scRNA.neftel/02_neftel.noPed.NP_2024.05.02.RDS")
DefaultAssay(obj) <- "RNA"

# extract annotation data for heatmaps downstream
meta.data <- as.data.frame(
  cbind(
    as.character(obj$tumor.id),
    as.character(obj$cell.state),
    as.character(obj$cell.cycle.phase)
  )
)
colnames(meta.data) <- c("tumor", "state", "phase")
rownames(meta.data) <- rownames(obj@meta.data)

write.csv(
  meta.data,
  "Output/Rdata/03_states.1_2024.05.02/00_meta.data_2023.09.15.csv"
)

################################################################################
# CREATE EXPRESSION MATRIX FOR CLUSTERING

# create function to binarize data
binarize <- function(val) {
  ifelse(val > 0, 1, 0)
}

# set identity to tumor ID
Idents(obj) <- obj$tumor.id
tumorID <- unique(obj$tumor.id)

# loop through tumors and calculate essential gene expression matrix, binarized
input.list <- list()
scale.list <- list()

for (i in 1:length(tumorID)) {
  # subset seurat object, scale data
  sub <- subset(obj, idents = tumorID[i])
  sub <- ScaleData(
    sub,
    features = ess$Gene,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    verbose = FALSE
  )
  # extract scaled expression data
  input <- as.data.frame(
    sub@assays$RNA@scale.data[which(rownames(sub@assays$RNA@scale.data) %in% ess$Gene), ]
  )
  scale.list[[i]] <- input
  # binarize and transpose
  input[] <- lapply(input, binarize)
  input <- t(input)
  # store
  input.list[[i]] <- input
}
names(input.list) <- tumorID
names(scale.list) <- tumorID

# save
saveRDS(
  input.list,
  "Output/Rdata/03_states.1_2024.05.02/01_expr.mtx_scaled.gds.trans.bin_xPT_2024.03.26.RDS"
)
saveRDS(
  scale.list,
  "Output/Rdata/03_states.1_2024.05.02/02_expr.mtx_scaled.gds_xPT_2024.03.26.RDS"
)

################################################################################
# CLUSTERING

# define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

dir <- paste0("Output/Rdata/03_states.1_2024.05.02/")
dir2 <- paste0("Output/Figures/03_states.1_2024.05.02/")

sub.list <- list()
for (i in 1:length(input.list)) {

  #### SET UP ####
  # function to compute agglomerative coefficient
  ac <- function(x) {
    agnes(input.list[[i]], method = x)$ac
  }
  # determine appropriate clustering linkage method (highest agglomerative coef)
  meth <- sapply(m, ac)
  meth.list <- c(meth[[1]], meth[[2]], meth[[3]], meth[[4]])
  x <- which(meth.list == max(meth.list))
  clm <- names(meth[x])

  # correlation/jaccard matrix for plotting
  set.seed(707)
  to.plot <- input.list[[i]] %>%
    simil(method = "Jaccard") %>%
    as.matrix()
  saveRDS(
    to.plot,
    paste0(dir, "02.1_to.plot_", names(input.list)[i], "_2024.05.02.RDS")
  )

  #### CLUSTER AND FILTERING ####
  set.seed(707)
  # calculate distance metric and perform clustering
  if (clm == "ward") {
    # cluster data
    clust.rows <- input.list[[i]] %>%
      dist(method = "binary") %>%
      hclust(method = "ward.D2")
    saveRDS(
      clust.rows,
      paste0(dir, "02.5_clust.rows_", names(input.list)[i], "_2024.03.26.RDS")
    )
    # iterate over multiple values of k ----
    set.seed(707)
    clusters <- as.data.frame(cutree(clust.rows, k = seq(2, 30, 1)))
    colnames(clusters) <- paste0("k_", colnames(clusters))
    # clustree(clusters, prefix = "k_") # plot clustertree
    # exclude clusters w/ < 5 cells or > 80% of cells in tumor
    num.cells <- length(rownames(input.list[[i]]))
    tbls <- sapply(clusters, table)
    for (t in 1:length(tbls)) {
      tst <- tbls[[t]] < 5
      if ("TRUE" %in% tst) {
        remove <- names(tbls[t])
        clusters[, which(colnames(clusters) == remove)] <- NULL
      }
      tst2 <- tbls[[t]]/num.cells > 0.8 # SHOULD I REMOVE THIS RESTRICTION #
      if ("TRUE" %in% tst2) {
        remove <- names(tbls[t])
        clusters[, which(colnames(clusters) == remove)] <- NULL
      }
    }
    # aggregate cluster info and metadata
    clust.data <- cbind(
      meta.data[which(rownames(meta.data) %in% rownames(input.list[[i]])), ],
      clusters
    )
    write.csv(
      clust.data,
      paste0(dir, "03_cluster.data_", names(input.list)[i], "_2024.03.26.csv")
    )

    # add various cluster data to seurat objects
    sub <- subset(obj, idents = names(input.list)[i])
    sub <- AddMetaData(sub, clusters)
    sub.list[[i]] <- sub

    # find all markers
    ks <- colnames(clusters)
    sigs.df <- data.frame()
    for (val in 1:length(ks)) {
      Idents(sub) <- sub@meta.data[, which(colnames(sub@meta.data) == ks[val])]
      deg <- FindAllMarkers(
        sub,
        assay = "RNA",
        logfc.threshold = 0, # ADJUST #
        test.use = "MAST",
        only.pos = TRUE,
        verbose = FALSE
      )
      write.csv(
        deg,
        paste0(dir, "04_deg_", names(input.list)[i], "_", ks[val], "_2024.03.26.csv")
      )
      # filter for p < 0.05 and p < 0.001
      # can't use same parameters as Neftel, bc they used t.test for deg
      deg.filt.p5 <- deg[which(deg$p_val_adj < 0.05), ] # ADJUST #
      print(length(rownames(deg.filt.p5)))
      deg.filt.p1 <- deg[which(deg$p_val_adj < 0.001), ] # ADJUST #
      # filter for np5 > 50 and np1 > 10
      t5 <- table(deg.filt.p5$cluster)
      t5.keep <- names(which(t5 >= 50)) # ADJUST #
      t1 <- table(deg.filt.p1$cluster)
      t1.keep <- names(which(t1 >= 1)) # ADJUST #
      # identify clusters to retain
      keep <- intersect(t5.keep, t1.keep)
      if (length(keep) != 0) {
        # isolate signatures for retained clusters (n5)
        sigs <- deg.filt.p5[which(deg.filt.p5$cluster %in% keep), ]
        sigs$k_value <- ks[val]
        sigs$k_clust <- paste0(sigs$k_value, "_", sigs$cluster)
        # make df of signatures for all k w/in tumor
        sigs.df <- rbind(sigs.df, sigs)
        write.csv(
          sigs,
          paste0(dir, "05_signatures.N5_", names(input.list)[i], "_", ks[val], "_2024.03.26.csv")
        )
      }
      else {
        print(paste(ks[val], "in tumor", names(input.list)[i], "has no significant signature."))
      }
    }
    if (length(rownames(sigs.df) != 0)) {

      write.csv(
        sigs.df,
        paste0(dir, "06_signatures.N5_", names(input.list)[i], "_all.k_2024.03.26.csv")
      )

      # sigs.df <- sigs.df[, -c(1:6, 8)]

      sigs.df <- sigs.df[, -c(1, 3:6, 8)]
      new <- sigs.df
      new$k_clust2 <- paste0(new$k_clust, "_", names(input.list)[i])
      new$k_clust <- new$k_clust2
      new$k_clust2 <- NULL
      write.csv(
        new,
        paste0(dir, "06.5_signatures.N5_w.FC_", names(input.list)[i], "_all.k_2024.03.26.csv")
      )

      sigs.df <- sigs.df[, -1]

      if (length(unique(sigs.df$k_clust)) > 1) {

        # pairs of clusters w/ jaccard > 0.75 -- removed cluster w/ lower N5
        pivot <- sigs.df %>%
          mutate(present = 1) %>%
          pivot_wider(names_from = k_clust, values_from = present, values_fill = 0) %>%
          as.data.frame()
        rownames(pivot) <- pivot$gene
        pivot$gene <- NULL
        for (col in 1:length(colnames(pivot))) {
          pivot[, col] <- as.integer(pivot[, col])
        }
        write.csv(
          pivot,
          paste0(dir, "07_jac.input_signatures.N5_", names(input.list)[i], "_all.k_2024.03.26.csv")
        )
        # similarities between retained clusters w/in tumor
        clust.sim <- t(pivot) %>%
          simil(method = "Jaccard") %>%
          as.matrix()
        find <- which(clust.sim > 0.75, arr.ind = TRUE) %>%
          as.data.frame()
        if (length(rownames(find)) > 2) {
          find <- find[-c(grep("\\.", rownames(find))), ]
          cols <- as.integer(names(which(table(find$col) > 1)))
          cols.to.keep <- colnames(clust.sim)[cols]
          sigs.df <- sigs.df[which(sigs.df$k_clust %in% cols.to.keep), ]
        }
        else { # i.e. if only two clusters (if only one cluster, then will be taken care of below)
          kp <- unique(sigs.df$k_clust)[1]
          sigs.df <- sigs.df[which(sigs.df$k_clust == kp), ]
        }
      }

      else {
        print(paste("Tumor", names(input.list)[i], "has only 1 significant signature."))
      }

      sigs.df$tumor <- names(input.list)[i]
      # FINAL, FILTERED N5 SIGANTURES / CLUSTER IN TUMOR 1
      write.csv(
        sigs.df,
        paste0(dir, "08_filt.fin.signatures.N5_", names(input.list)[i], "_2024.03.26.csv")
      )
    }
    else {
      print(
        paste(names(input.list)[i], "does not have any significant signatures.")
      )
    }

    # PLOT HEATMAP OF JACCARD SCORES, PHEATMAP ----
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
    # PLOT HEATMAP OF JACCARD SCORES, CORRPLOT
    pdf(paste0(dir2, "02_", names(input.list)[i], "_clustered.corrplot_2024.03.27.pdf"), height = 8, width = 9)
    corrplot(
      to.plot,
      tl.pos = "n",
      order = "hclust",
      hclust.method = "ward.D2",
      method = "shade",
      type = "full",
      col = viridis_pal(option = "inferno")(100),
      is.corr = FALSE,
      diag = TRUE,
      na.label = "square", na.label.col = "grey"
    )
    dev.off()

  }
  else {
    print(
      paste0(
        "The best linkage method for ", names(input.list)[i], " is ", clm, "not Ward"
      )
    )
  }
}
