library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(viridis)
library(limma)
library(EnhancedVolcano)

source("data/r.input_run_lm_stats_limma_FUNCTION_2023.03.23.R")
source("data/r.input_run_lm_fit_only_2023.05.31.R")

# dir.create("Output")
# dir.create("Output/Rdata")
# dir.create("Output/Figures")
# dir.create("Output/Figures/QC")

################################################################################
# READ DATA

efx <- read.delim(
  file = "data/r.input_CRISPRGeneEffect.csv",
  sep = ",",
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  file = "data/r.input_metadata_Model.csv",
  sep = ",",
  row.names = 1
)

# com.ess <- read.delim(
#   file = "r.input_CRISPRInferredCommonEssentials.csv",
#   sep = ","
# )

################################################################################
# CLEAN DATA

# rewrite column names of efx data to be Gene ID only (remove numbers trailing)
colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))

# filter metadata for only those with effect score data
meta.filt <- meta[match(rownames(efx), rownames(meta)), ]

# # rewrite common essential genes to be Gene ID only (remove numbers trailing)
# com.ess <- gsub(" \\([0-9]*\\)", "", com.ess$Essentials)

################################################################################
# GLIOBLASTOMA

# make vector of 0/1 indicating which cell lines are Glioblastoma
gbm.ind.var <- as.numeric(meta.filt$OncotreeSubtype == "Glioblastoma")

# run limma, just lmfit to acquire residuals
gbm.fit <- run_lm_fit(
  mat = efx %>% as.matrix(),
  vec = gbm.ind.var
)
# run limma, complete analysis to acquire essential genes
gbm.stats <- run_lm_stats_limma(
  mat = efx %>% as.matrix(),
  vec = gbm.ind.var
)

saveRDS(
  gbm.fit,
  file = "Output/Rdata/00_depmap.gbm.fit_2023.09.07.rds"
)

write.table(
  gbm.stats,
  file = "Output/Rdata/00_depmap.gbm.stats_2023.03.23.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# qq plot for quality control
pdf("Output/Figures/QC/GBM_QQ_2023.05.31.pdf", height = 7, width = 7)
qqt(gbm.fit$t, df = gbm.fit$df.prior+gbm.fit$df.residual, pch = 16, cex = 0.2)
abline(0, 1)
dev.off()

# volcano plot of limma results
pdf("Output/Figures/01_GBM_Volcano.Essential.Genes_2023.05.31.pdf", height = 7, width = 7)
EnhancedVolcano(
  gbm.stats,
  lab = gbm.stats$Gene,
  x = "EffectSize",
  y = "q.value",
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  boxedLabels = TRUE,
  pCutoff = 1e-2,
  FCcutoff = 0,
  pointSize = 2,
  labSize = 6,
  labCol = "black",
  labFace = "bold"
)
dev.off()

#### acquire mean effect score for every gene in gbm cell lines

# subset metadata
gbm <- subset(
  meta.filt,
  meta.filt$OncotreeSubtype == "Glioblastoma"
)
# subset effects data
gbm.efx <- subset(
  efx,
  rownames(efx) %in% rownames(gbm)
)

# calculate mean effect score
gbm.mtx <- data.frame()
for (i in 1:ncol(gbm.efx)) {
  mn <- mean(as.numeric(gbm.efx[, i]))
  gene <- colnames(gbm.efx)[i]
  gbm.mtx[i, 1] <- gene
  gbm.mtx[i, 2] <- mn
}
colnames(gbm.mtx) <- c("Gene", "mEffectScore")

write.table(
  gbm.mtx,
  file = "Output/Rdata/01_depmap.gbm.mean.effect.scores_2023.08.08.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

################################################################################
# SAVE GENE LISTS

# add mean effect score to gbm.stats output
gbm.stats_m <- merge(gbm.stats, gbm.mtx, by = "Gene")

# true essential genes (mean effect score <-0.5, lower in gbm lines, no common)

lower <- gbm.stats_m[which(gbm.stats_m$EffectSize < 0), ] # all lower 8930
sig.lower <- lower[which(lower$q.value < 5e-2), ] # all sig lower 377
killing <- sig.lower[which(sig.lower$mEffectScore < -0.5), ] # all sig lower w/ killing effect 168
# true <- killing[which(!killing$Gene %in% com.ess), ] # no common essential 23

# save all gene lists

write.table(
  lower,
  file = "Output/Rdata/02_gbm.all.lower_2023.09.07.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

write.table(
  sig.lower,
  file = "Output/Rdata/03_gbm.all.sig.lower_2023.09.07.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

write.table(
  killing,
  file = "Output/Rdata/04_gbm.killing_2023.09.07.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# write.table(
#   true,
#   file = "Output/Rdata/05_gbm.true.ess_2023.09.07.csv",
#   quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
# )
