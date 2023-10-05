library(tidyverse)
library(Seurat)
# library(ggplot2)
# library(ggsci)
# library(ggpubr)
# library(ggrepel)
# library(viridis)
# library(EnhancedVolcano)
# library(MAST)
# library(SummarizedExperiment)

dir.create("DEPMAP/Output/Rdata/scRNA")

################################################################################
# READ DATA

gbm.seurat <- readRDS("DepMap.Scored_scGBM_2023.02.21.RDS")

# lower <- read.table(
#  "Output/Rdata/02_gbm.all.lower_2023.09.07.csv",
#  row.names = 1, header = TRUE, sep = ","
# )
# 
# sig.lower <- read.table(
#  "Output/Rdata/03_gbm.all.sig.lower_2023.09.07.csv",
#  row.names = 1, header = TRUE, sep = ","
# )
# 
# killing <- read.table(
#  "Output/Rdata/04_gbm.killing_2023.09.07.csv",
#  row.names = 1, header = TRUE, sep = ","
# )
# 
# true <- read.table(
#  "Output/Rdata/05_gbm.true.ess_2023.09.07.csv",
#  row.names = 1, header = TRUE, sep = ","
# )
# 
# ################################################################################
# # ANNOTATE
# 
# Idents(gbm.seurat)  <- gbm.seurat$ClassificationNPvNon
# gbm.seurat <- RenameIdents(gbm.seurat, 'MES1' = 'MES')
# gbm.seurat <- RenameIdents(gbm.seurat, 'MES2' = 'MES')
# gbm.seurat <- RenameIdents(gbm.seurat, 'NPC1' = 'NPC')
# gbm.seurat <- RenameIdents(gbm.seurat, 'NPC2' = 'NPC')
# gbm.seurat[["CellStateClassification4"]] <- Idents(gbm.seurat)
# 
# gbm.seurat@meta.data <- gbm.seurat@meta.data %>%
#  mutate(
#    pt.split = case_when(
#    NPvNon == "Neoplastic" ~ paste0(orig.ident, "_", CellStateClassification4),
#    NPvNon == "Non-Neoplastic" ~ "Non-Neoplastic"
#  )
# )
# 
# ################################################################################
# # DIFFERENTIAL EXPRESSION
# 
# Idents(gbm.seurat) <- gbm.seurat@meta.data$NPvNon
# 
# de <- FindMarkers(
#    gbm.seurat,
#    ident.1 = "Neoplastic", # finding markers for neoplastic cells (positive FC)
#    ident.2 = "Non-Neoplastic",
#    test.use = "MAST",
#    logfc.threshold = 0.25,
#    only.pos = TRUE,
#    assay = "RNA"
# )
# 
# write.table(
#    de,
#    file = "Output/Rdata/scRNA/DEG_gbm_suter_2023.09.07.csv",
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
# 
# de <- de[which(de$p_val_adj < 0.01), ]
# 
# # filtering essential gene lists to contain only DEG
# lower.ref <- lower[which(lower$Gene %in% de$gene), ]
# write.table(
#    lower.ref,
#    file = "Output/Rdata/scRNA/02_REF_gbm.all.lower_2023.09.07.csv",
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
# 
# sig.lower.ref <- sig.lower[which(sig.lower$Gene %in% de$gene), ]
# write.table(
#    sig.lower.ref,
#    file = "Output/Rdata/scRNA/03_REF_gbm.all.sig.lower_2023.09.07.csv",
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
# killing.ref <- killing[which(killing$Gene %in% de$gene), ]
# write.table(
#    killing.ref,
#    file = "Output/Rdata/scRNA/04_REF_gbm.killing_2023.09.07.csv",
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
# 
# true.ref <- true[which(true$Gene %in% de$gene), ]
# write.table(
#    true.ref,
#    file = "Output/Rdata/scRNA/05_REF_gbm.true.ess_2023.09.07.csv",
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
# 
# ################################################################################
# # SCORE GENE SETS
# 
# # score unrefined sets
# 
# set.list <- list(
#    lower$Gene,
#    sig.lower$Gene,
#    killing$Gene,
#    true$Gene
# )
# name.list <- list(
#    "DM_lower_score_",
#    "DM_sig.lower_score_",
#    "DM_killing_score_",
#    "DM_true_score_"
# )
# 
# for (i in 1:length(set.list)) {
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat, features = list(set.list[[i]]),
#        name = paste0(name.list[i])
#    )
# }
# 
# # score refined sets ("ALL")
# 
# set.list.all <- list(
#    lower.ref$Gene,
#    sig.lower.ref$Gene,
#    killing.ref$Gene,
#    true.ref$Gene
# )
# 
# name.list.all <- list(
#    "DM_ref_lower_score_",
#    "DM_ref_sig.lower_score_",
#    "DM_ref_killing_score_",
#    "DM_ref_true_score_"
# )
# 
# for (i in 1:length(set.set.list.all)) {
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat, features = list(set.list.all[[i]]),
#        name = paste0(name.list.all[i])
#    )
# }

################################################################################
# OUTPUT SCALED GENE EXPRESSION, NP CELLS ONLY

Idents(gbm.seurat) <- gbm.seurat$NPvNon
gbm.np <- subset(gbm.seurat, idents = "Neoplastic")

# extract scaled expression data and transpose (cells = rows, genes = columns)
scaled <- t(gbm.np@assays$RNA@scale.data)

write.table(
    scaled, 
    file = "Output/Rdata/scRNA/scaled.rna.expr_gbm.np_suter_2023.09.28.csv",
    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# ################################################################################
# # BY PATIENT
# 
# dir.create("Output/Rdata/scRNA/xPT")
# 
# ## CREATE LIST OF TUMOR OBJECTS
# Idents(gbm.seurat) <- gbm.seurat@meta.data$pt.split
# 
# suter.tumors <- c("GBM21", "GBM41", "GBM47", "GBM49", "GBM51", "GBM53")
# suter.obj.list <- list()
# for (i in 1:length(suter.tumors)) {
#  tum <- subset(
#    gbm.seurat,
#    idents = c(
#      paste0(suter.tumors[i], "_MES"),
#      paste0(suter.tumors[i], "_NPC"),
#      paste0(suter.tumors[i], "_OPC"),
#      paste0(suter.tumors[i], "_AC"),
#      "Non-Neoplastic"
#    )
#  )
#  suter.obj.list[[length(suter.obj.list) + 1]] <- assign(suter.tumors[i], tum)
#  rm(tum)
# }
# 
# saveRDS(suter.obj.list, "Output/Rdata/scRNA/xPT/xPT_Suter_obj.list_2023.06.01.RDS")
# 
# # DIFFERENTIAL EXPRESSION
# 
# markers.list <- list()
# for (i in 1:length(suter.obj.list)) {
#  obj <- suter.obj.list[[i]]
#  Idents(obj) <- obj@meta.data$NPvNon
#  markers <- FindMarkers(
#    obj,
#    ident.1 = "Neoplastic",
#    ident.2 = "Non-Neoplastic",
#    test.use = "MAST",
#    logfc.threshold = 0.25,
#    only.pos = TRUE,
#    assay = "RNA"
#  )
#  write.table(
#    markers,
#    file = paste0("Output/Rdata/scRNA/xPT/DEG_", suter.tumors[i], "_2023.05.26.csv"),
#    quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
#  )
#  markers <- markers[which(markers$p_val_adj < 0.01), ] # remove any genes with p > 0.01
#  markers.list[[length(markers.list) + 1]] <- assign(
#    paste0(suter.tumors[i], "_markers"),
#    markers
#  )
#  rm(obj)
#  rm(markers)
# }
# saveRDS(
#    markers.list,
#    file = "Output/Rdata/scRNA/xPT/DEG_xPT_list_2023.09.11.rds"
# )
# 
# ## ESSENTIAL GENES FOR EACH TUMOR
# 
# ess.xPT.list <- list()
# for (i in 1:length(markers.list)) {
# 
#    ess.list <- list()
# 
#    lower.genes <- lower[which(lower$Gene %in% markers.list[[i]]$gene), ]
#    write.table(
#        lower.genes,
#        file = paste0("Output/Rdata/scRNA/xPT/02_REF_gbm.all.lower_", suter.tumors[i], "_2023.09.11.csv"),
#        quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
#    )
#    ess.list[[1]] <- assign("lower_ref", lower.genes)
# 
#    sig.lower.genes <- sig.lower[which(sig.lower$Gene %in% markers.list[[i]]$gene), ]
#    write.table(
#        sig.lower.genes,
#        file = paste0("Output/Rdata/scRNA/xPT/03_REF_gbm.all.sig.lower_", suter.tumors[i], "_2023.09.11.csv"),
#        quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
#    )
#    ess.list[[2]] <- assign("sig.lower_ref", sig.lower.genes)
# 
#    killing.genes <- killing[which(killing$Gene %in% markers.list[[i]]$gene), ]
#    write.table(
#        killing.genes,
#        file = paste0("Output/Rdata/scRNA/xPT/04_REF_gbm.killing_", suter.tumors[i], "_2023.09.11.csv"),
#        quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
#    )
#    ess.list[[3]] <- assign("killing_ref", killing.genes)
# 
#    true.genes <- true[which(true$Gene %in% markers.list[[i]]$gene), ]
#    write.table(
#        true.genes,
#        file = paste0("Output/Rdata/scRNA/xPT/05_REF_gbm.true.ess_", suter.tumors[i], "_2023.09.11.csv"),
#        quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
#    )
#    ess.list[[4]] <- assign("true_ref", true.genes)
# 
# 
#    ess.xPT.list[[length(markers.list) + 1]] <- assign(
#        paste0(suter.tumors[i], "_ref"),
#        ess.list
#    )
# 
# }
# saveRDS(
#    ess.xPT.list,
#    file = "Output/Rdata/scRNA/xPT/00_REF_xPT_list_2023.09.11.rds"
# )
# 
# ## ADD MODULE SCORE
# 
# name.list.xPT <- names(ess.xPT.list)
# name.list.xPT.2 <- names(ess.xPT.list[[1]])
# 
# for (i in 1:length(ess.xPT.list)) {
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat,
#        features = list(ess.xPT.list[[i]][[1]]),
#        name = paste0(name.list.xPT[i], "_", name.list.xPT.2[1])
#    )
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat,
#        features = list(ess.xPT.list[[i]][[2]]),
#        name = paste0(name.list.xPT[i], "_", name.list.xPT.2[2])
#    )
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat,
#        features = list(ess.xPT.list[[i]][[3]]),
#        name = paste0(name.list.xPT[i], "_", name.list.xPT.2[3])
#    )
#    gbm.seurat <- AddModuleScore(
#        gbm.seurat,
#        features = list(ess.xPT.list[[i]][[4]]),
#        name = paste0(name.list.xPT[i], "_", name.list.xPT.2[4])
#    )
# }

# ###############################################################################
# SAVE RDS
# 
# saveRDS(
#  gbm.seurat,
#  file = "DepMap.Scored_scGBM_2023.02.21.RDS"
# )
