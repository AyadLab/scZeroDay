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

dir.create("Output/Figures/09_signature.corr")
dir.create("Output/Rdata/09_signature.corr")

################################################################################
# READ IN SIGNATURES

# vs signatures ----
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

# neftel signatures ----
suva.sigs <- read.csv("data/Suva_States_2023.11.09.csv")
glimpse(suva.sigs)
suva <- list()
nm <- c()
for (i in 1:length(unique(suva.sigs$Cell.State))) {
  genes <- suva.sigs$Gene[which(suva.sigs$Cell.State == unique(suva.sigs$Cell.State)[i])]
  suva[[i]] <- genes
  nm <- c(nm, unique(suva.sigs$Cell.State)[i])
}
names(suva) <- nm
glimpse(suva)

# verhaak signatures ----
ver.sigs <- read.csv("data/Verhaak.States.csv")
glimpse(ver.sigs)
verhaak <- list(
  Mesenchymal = ver.sigs$Mesenchymal[nzchar(ver.sigs$Mesenchymal)],
  Proneural = ver.sigs$Proneural[nzchar(ver.sigs$Proneural)],
  Neural = ver.sigs$Neural[nzchar(ver.sigs$Neural)],
  Classical = ver.sigs$Classical[nzchar(ver.sigs$Classical)]
)
glimpse(verhaak)

# ravi signatures ----
ravi.sigs <- read.csv("data/ravi_sigs.csv")
glimpse(ravi.sigs)
ravi <- list()
nm <- c()
for (i in 1:length(unique(ravi.sigs$Regional.Programm))) {
  genes <- ravi.sigs$Gene[which(ravi.sigs$Regional.Programm == unique(ravi.sigs$Regional.Programm)[i])]
  ravi[[i]] <- genes
  nm <- c(nm, unique(ravi.sigs$Regional.Programm)[i])
}
names(ravi) <- nm
glimpse(ravi)

# make master list ----
master <- c(vs, suva, verhaak, ravi)

# kandoth signatures ----
kan.sigs <- read.csv("data/kandoth_sigs.csv")
glimpse(kan.sigs)
kandoth <- list(
  TF.Regulator = kan.sigs$TF.Regulator[nzchar(kan.sigs$TF.Regulator)],
  Histone.Mod = kan.sigs$Histone.Mod[nzchar(kan.sigs$Histone.Mod)],
  Genome.Integ = kan.sigs$Genome.Integ[nzchar(kan.sigs$Genome.Integ)],
  RTK.Sig = kan.sigs$RTK.Sig[nzchar(kan.sigs$RTK.Sig)],
  Cell.Cycle = kan.sigs$Cell.Cycle[nzchar(kan.sigs$Cell.Cycle)],
  MAPK.Sig = kan.sigs$MAPK.Sig[nzchar(kan.sigs$MAPK.Sig)],
  PI3K.Sig = kan.sigs$PI3K.Sig[nzchar(kan.sigs$PI3K.Sig)],
  TGFB.Sig = kan.sigs$TGFB.Sig[nzchar(kan.sigs$TGFB.Sig)],
  WNT.Bcat.Sig = kan.sigs$WNT.Bcat.Sig[nzchar(kan.sigs$WNT.Bcat.Sig)],
  Histone = kan.sigs$Histone[nzchar(kan.sigs$Histone)],
  Proteolysis = kan.sigs$Proteolysis[nzchar(kan.sigs$Proteolysis)],
  Splicing = kan.sigs$Splicing[nzchar(kan.sigs$Splicing)],
  HIPPO.Sig = kan.sigs$HIPPO.Sig[nzchar(kan.sigs$HIPPO.Sig)],
  DNA.me = kan.sigs$DNA.me[nzchar(kan.sigs$DNA.me)],
  Metabolism = kan.sigs$Metabolism[nzchar(kan.sigs$Metabolism)],
  NFE2L = kan.sigs$NFE2L[nzchar(kan.sigs$NFE2L)],
  Protein.Phosphatase = kan.sigs$Protein.Phosphatase[nzchar(kan.sigs$Protein.Phosphatase)],
  Ribosome = kan.sigs$Ribosome[nzchar(kan.sigs$Ribosome)],
  TOR.Sig = kan.sigs$TOR.Sig[nzchar(kan.sigs$TOR.Sig)],
  Other = kan.sigs$Other[nzchar(kan.sigs$Other)]
)
glimpse(kandoth)

################################################################################
# SINGLE CELL DATASET(S)

sut.np <- readRDS(
  "data/Suter_scGBM.NP_2023.10.04.RDS"
)

################################################################################
# MODULE SCORE

# DO IN ISO OBJECT?
for (i in 1:length(master)) {
  sut.np <- AddModuleScore(
    sut.np,
    features = list(master[[i]]),
    name = paste0("corr_", names(master)[i], "_"),
    slot = "data"
  )
}

################################################################################
# SINGSCORE

obj.data <- GetAssayData(sut.np, assay = "RNA", layer = "scale.data")
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
rm(obj.data)
rm(obj.data.df)
gc()

sing.names.list <- c()
# score each cell for each cell state
for (i in 1:length(master)) {
  # capture unique gene signature
  sig <- master[[i]]
  # capture name of cell state ,
  state <- paste0(names(master)[i])
  # singscore scoring
  set <- GeneSet()
  set@geneIds <- as.character(sig)
  scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
  # rename $Sig to appropriate cell state
  scored$Sig <- as.character(state)
  assign(paste0("scored.", state), scored)

  sing.names.list[i] <- paste0("scored.", state)

  # append singscores to seurat object
  sut.np <- AddMetaData(
    sut.np,
    metadata = scored[, "TotalScore"],
    col.name = paste0(state, "_singscore")
  )

}

saveRDS(
  sut.np,
  "data/Suter_scGBM.NP_2023.10.04.RDS"
)

################################################################################
# CORRELATION ANALYSIS

cor.df <- data.frame(
  row.names = rownames(sut.np@meta.data),
  VS1 = sut.np$VS1_singscore,
  VS2 = sut.np$VS2_singscore,
  VS3 = sut.np$VS3_singscore,
  Neftel.MES2 = sut.np$MES2_singscore,
  Neftel.MES1 = sut.np$MES1_singscore,
  Neftel.AC = sut.np$AC_singscore,
  Neftel.OPC = sut.np$OPC_singscore,
  Neftel.NPC1 = sut.np$NPC1_singscore,
  Neftel.NPC2 = sut.np$NPC2_singscore,
  Verhaak.Mes = sut.np$Mesenchymal_singscore,
  Verhaak.Pro = sut.np$Proneural_singscore,
  Verhaak.Neu = sut.np$Neural_singscore,
  Verhaak.Cla = sut.np$Classical_singscore,
  Ravi.Glia = sut.np$Radial.Glia_singscore,
  Ravi.rImm = sut.np$Reactive.Immune_singscore,
  Ravi.rNPC = sut.np$Regional.NPC_singscore,
  Ravi.rOPC = sut.np$Regional.OPC_singscore,
  Ravi.rHyp = sut.np$Reactive.Hypoxia_singscore
)
mtx <- as.matrix(cor.df)
saveRDS(mtx, "Output/Rdata/09_signature.corr/01_corr.df_2024.05.09.RDS")

# make correlation
cor.mtx <- cor(
  mtx,
  use = "everything",
  method = "spearman"
)
saveRDS(cor.mtx, "Output/Rdata/09_signature.corr/02_corr.mtx_2024.05.09.RDS")

# corrplot, lower

pdf("Output/Figures/09_signature.corr/01_sigs.corrplot_lower_2024.05.09.pdf", height = 8, width = 8)
corrplot(
  cor.mtx,
  tl.pos = "l", tl.col = "black",
  order = "hclust",
  hclust.method = "ward.D2",
  method = "circle",
  type = "lower", # addrect = 7, rect.col = "red",
  col = colorRampPalette(paletteer_d(palette = "RColorBrewer::PiYG", 11))(100),
  is.corr = TRUE,
  diag = TRUE,
  na.label = "square", na.label.col = "grey"
)
dev.off()

# corrplot, complete w/ boxes
pdf("Output/Figures/09_signature.corr/02_sigs.corrplot_full.boxed_2024.05.09.pdf", height = 8, width = 8)
corrplot(
  cor.mtx,
  tl.pos = "l", tl.col = "black",
  order = "hclust",
  hclust.method = "ward.D2",
  method = "color",
  type = "full", addrect = 5, rect.col = "red",
  col = colorRampPalette(paletteer_d(palette = "RColorBrewer::PiYG", 11))(100),
  is.corr = TRUE,
  diag = TRUE,
  na.label = "square", na.label.col = "grey"
)
dev.off()

# pheatmap

# my.breaks <- c(seq(min(cor.mtx), 0, length.out=ceiling(100/2) + 1),
#                seq(max(cor.mtx)/100, max(cor.mtx), length.out=floor(100/2)))
# pheatmap::pheatmap(
#   cor.mtx,
#   color = colorRampPalette(paletteer_d(palette = "RColorBrewer::PRGn", 11))(100),
#   # annotation_row = rownames(cor.mtx),
#   # annotation_col = an.col,
#   # annotation_colors = an.colors,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   show_rownames = TRUE,
#   show_colnames = FALSE,
#   annotation_names_col = FALSE,
#   use_raster = FALSE,
#   breaks = my.breaks
# )

################################################################################
# KANDOTH SINGSCORE

obj.data <- GetAssayData(sut.np, assay = "RNA", layer = "scale.data")
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
rm(obj.data)
rm(obj.data.df)
gc()

sing.names.list <- c()
# score each cell for each cell state
for (i in 1:length(kandoth)) {
  # capture unique gene signature
  sig <- kandoth[[i]]
  # capture name of cell state ,
  state <- paste0(names(kandoth)[i])
  # singscore scoring
  set <- GeneSet()
  set@geneIds <- as.character(sig)
  scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
  # rename $Sig to appropriate cell state
  scored$Sig <- as.character(state)
  assign(paste0("scored.", state), scored)

  sing.names.list[i] <- paste0("scored.", state)

  # append singscores to seurat object
  sut.np <- AddMetaData(
    sut.np,
    metadata = scored[, "TotalScore"],
    col.name = paste0(state, "_singscore")
  )

}

saveRDS(
  sut.np,
  "data/Suter_scGBM.NP_kandoth_2024.05.10.RDS"
)

################################################################################
# KANDOTH CORRELATION ANALYSIS

names(kandoth)
kan.df <- data.frame(
  row.names = rownames(sut.np@meta.data),
  VS1 = sut.np$VS1_singscore,
  VS2 = sut.np$VS2_singscore,
  VS3 = sut.np$VS3_singscore,
  TF.Regulator = sut.np$TF.Regulator_singscore,
  Histone.Mod = sut.np$Histone.Mod_singscore,
  Genome.Integ = sut.np$Genome.Integ_singscore,
  RTK.Sig = sut.np$RTK.Sig_singscore,
  Cell.Cycle = sut.np$Cell.Cycle_singscore,
  MAPK.Sig = sut.np$MAPK.Sig_singscore,
  PI3K.Sig = sut.np$PI3K.Sig_singscore,
  TGFB.Sig = sut.np$TGFB.Sig_singscore,
  WNT.Bcat.Sig = sut.np$WNT.Bcat.Sig_singscore,
  Histone = sut.np$Histone_singscore,
  Proteolysis = sut.np$Proteolysis_singscore,
  Splicing = sut.np$Splicing_singscore,
  HIPPO.Sig = sut.np$HIPPO.Sig_singscore,
  DNA.me = sut.np$DNA.me_singscore,
  Metabolism = sut.np$Metabolism_singscore,
  NFE2L = sut.np$NFE2L_singscore,
  Protein.Phosphatase = sut.np$Protein.Phosphatase_singscore,
  Ribosome = sut.np$Ribosome_singscore,
  TOR.Sig = sut.np$TOR.Sig_singscore,
  Other = sut.np$Other_singscore
)
mtx <- as.matrix(kan.df)
saveRDS(mtx, "Output/Rdata/09_signature.corr/01_kan.df_2024.05.09.RDS")

# make correlation
kan.mtx <- cor(
  mtx,
  use = "everything",
  method = "spearman"
)
saveRDS(kan.mtx, "Output/Rdata/09_signature.corr/02_kan.mtx_2024.05.09.RDS")

# corrplot, lower

pdf("Output/Figures/09_signature.corr/03_kan.corrplot_lower_2024.05.09.pdf", height = 8, width = 8)
corrplot(
  kan.mtx,
  tl.pos = "l", tl.col = "black",
  order = "hclust",
  hclust.method = "ward.D2",
  method = "circle",
  type = "lower", # addrect = 7, rect.col = "red",
  col = colorRampPalette(paletteer_d(palette = "RColorBrewer::PiYG", 11))(100),
  is.corr = TRUE,
  diag = TRUE,
  na.label = "square", na.label.col = "grey"
)
dev.off()

# corrplot, complete w/ boxes
pdf("Output/Figures/09_signature.corr/04_kan.corrplot_full.boxed_2024.05.09.pdf", height = 8, width = 8)
corrplot(
  kan.mtx,
  tl.pos = "l", tl.col = "black",
  order = "hclust",
  hclust.method = "ward.D2",
  method = "color",
  type = "full", addrect = 4, rect.col = "red",
  col = colorRampPalette(paletteer_d(palette = "RColorBrewer::PiYG", 11))(100),
  is.corr = TRUE,
  diag = TRUE,
  na.label = "square", na.label.col = "grey"
)
dev.off()













################################################################################


# # add state identities to each cell for each gene set (4 sets)
# singscore.df <- data.frame(
#   row.names = rownames(obj.np@meta.data),
#   mVC1 = obj.np@meta.data$Cluster_1singscore,
#   mVC2 = obj.np@meta.data$Cluster_2singscore,
#   mVC3 = obj.np@meta.data$Cluster_3singscore
# )
#
# assignments <- data.frame(
#   row.names = rownames(obj.np@meta.data)
# )
# assignments$cell.state <- "Ambiguous"
#
# for (i in 1:length(rownames(singscore.df))) {
#   state <- which.max(singscore.df[i, ])
#   assignments$cell.state[i] <- names(state)
# }
#
# # add cell states to seurat object (np only)
# obj.np <- AddMetaData(obj.np, metadata = assignments, col.name = "mVC")
# table(obj.np@meta.data$mVC)
#
# # add cell states to seurat object (all)
# obj <- AddMetaData(obj, metadata = assignments, col.name = "mVC")
# obj@meta.data <- obj@meta.data %>%
#   mutate(
#     mVC = case_when(
#       neoplastic.state == "Neoplastic" ~ mVC,
#       neoplastic.state == "Non-Neoplastic" ~ "Non-Neoplastic"
#     )
#   )
# table(obj@meta.data$mVC)
#
#
