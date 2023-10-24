library(tidyverse)
library(viridis)
library(pheatmap)
library(proxy)
library(corrplot)
library(factoextra)
library(FactoMineR)
library(NbClust)
library(cowplot)
library(cluster)

setwd("/data/mrd/cognitive.seeds")

tumor.mtx <- readRDS(
  "Output/Rdata/jaccard/03_jaccard.mtx_xPT_2023.010.10.RDS"
)

################################################################################
# IDENTIFY OPTIMAL K FOR EACH TUMOR

set.seed(707)

# determine optimal number of clusters for kmeans clustering for every tumor

jac <- tumor.mtx[[4]]
jac[is.na(jac)] <- 1

elbow <- fviz_nbclust(jac, kmeans, method = "wss", k.max = 25) +
  theme_minimal() +
  ggtitle("Elbow Method")

# gap statistic calculation takes too long when too many samples
gap_stat <- clusGap(jac, FUN = kmeans, K.max = 20, B = 50)
gap <- fviz_gap_stat(gap_stat) +
  theme_minimal() +
  ggtitle("Gap Statistic Method")


# silhouette <- fviz_nbclust(jac, kmeans, method = "silhouette", k.max = 12) +
#   theme_minimal() +
#   ggtitle("Silhouette Method")

pdf(paste0("Output/Figures/QC/jaccard.optimal.clusters_", names(tumor.mtx)[4], "_2023.09.20.pdf"), height = 7, width = 7)
print(plot_grid(elbow, gap))
dev.off()



