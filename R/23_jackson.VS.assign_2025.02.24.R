# WHAT IS THIS SCRIPT FOR?? NEEDED FOR JACKSON DATA??

library(tidyverse)
library(Seurat)
library(future)
library(singscore)
library(GSEABase)
library(presto)
library(cowplot)

options(future.globals.maxSize = +Inf)

################################################################################
obj <- readRDS(
  "/data/mrd/misc/data/jackson_GSE278456_TumorSeurat.RDS"
)

# ess <- read.csv(
#   "/data/mrd/cognitive.seeds_zero/Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv",
#   row.names = 1
# )
# 
# # vs markers
# markers <- readRDS("/data/mrd/cognitive.seeds_zero/Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS")
# names(markers)
# 
# # gds genes in venn
# gdsXvenn <- read.csv(
#   "/data/mrd/cognitive.seeds_zero/Output/Rdata/19_gds.dive/gdsXvenn_genes_2024.11.21.csv",
#   row.names = 1
# )
# 
# top 100 venn genes
topVenn <- read.csv(
  "/data/mrd/misc/data/04_VS.top100.overlapping_2024.11.12.csv",
  row.names = 1
)
glimpse(topVenn)

# # all venn genes
# venn <- read.csv(
#   "/data/mrd/cognitive.seeds_zero/Output/Rdata/10_deg.corr/04_VS.overlapping.genes_2024.11.12.csv",
#   row.names = 1
# )


################################################################################
# CAN'T RE-PROCESS DATA...DON'T HAVE TUMOR-LEVEL INFORMATION

plan("multisession", workers = 10)
obj <- ScaleData(
  obj,
  features = rownames(obj[["RNA"]]$counts),
  vars.to.regress = c("pct.mito", "S.Score", "G2M.Score")
)
obj.data <- obj[["RNA"]]$scale.data
saveRDS(obj.data, "/data/mrd/misc/Output/jackson_scale.data_2025.02.03.RDS")
plan("sequential")

################################################################################
# set up ranked data based on scRNA-seq signatures
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
# rm(obj.data)
rm(obj.data.df)
gc()

# score each cell for each cell state
for (i in 1:length(unique(topVenn$cluster))) {
  # capture unique gene signature
  sig <- topVenn$merge[which(topVenn$cluster == unique(topVenn$cluster)[i])]
  # capture name of cell state ,
  state <- unique(topVenn$cluster)[i]
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
  "scored.mVC1",
  "scored.mVC2",
  "scored.mVC3"
)

# append singscores to seurat object
for (n in singscores.list) {
  state <- unlist(strsplit(as.character(n), split = ".", fixed = TRUE))[2]
  scored <- eval(parse(text = n))
  obj <- AddMetaData(
    obj,
    metadata = scored[, "TotalScore"],
    col.name = paste0(state, "_singscore")
  )
}

# add state identities to each cell
singscore.df <- data.frame(
  row.names = rownames(obj@meta.data),
  mVC1 = obj@meta.data$mVC1_singscore,
  mVC2 = obj@meta.data$mVC2_singscore,
  mVC3 = obj@meta.data$mVC3_singscore
)

assignments <- data.frame(
  row.names = rownames(obj@meta.data)
)
assignments$cell.state <- "Ambiguous"

for (i in 1:length(rownames(singscore.df))) {
  state <- which.max(singscore.df[i, ])
  assignments$cell.state[i] <- names(state)
}

# add cell states to seurat object
obj <- AddMetaData(obj, metadata = assignments, col.name = "VS")
table(obj@meta.data$VS)
# majority VS1, aligns with RNAseq assignments!



# dir.create("UM002/Output")
saveRDS(
  obj@meta.data,
  "/data/mrd/misc/Output/jackson_meta.data_2025.02.03.RDS"
)

