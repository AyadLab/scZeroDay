library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(viridis)
library(limma)
library(EnhancedVolcano)
library(paletteer)
library(plotly)
library(msigdbr)
library(clusterProfiler)
library(RColorBrewer)
library(enrichR)

source("data/r.input_run_lm_stats_limma_FUNCTION_2023.03.23.R")
source("data/r.input_run_lm_fit_only_2023.05.31.R")

dir.create("Output")
dir.create("Output/Rdata")
dir.create("Output/Rdata/01_essential.genes")
dir.create("Output/Figures")
dir.create("Output/Figures/01_essential.genes")
dir.create("Output/Figures/01_essential.genes/enrichR")

# figure theme
mrd.theme <- theme(
  title = element_blank(),
  plot.title = element_blank(),
  plot.subtitle = element_blank(), plot.caption = element_blank(),
  legend.position = "none",
  axis.text = element_text(size = 10, colour = "black"),
  axis.text.x = element_text(size = 10, colour = "black"),
  axis.text.y = element_text(size = 10, colour = "black"),
  axis.title = element_text(size = 12, colour = "black"),
  axis.line = element_line(linewidth = 1, colour = "black"),
  axis.ticks.length = unit(0.1, "cm"),
  axis.ticks = element_line(linewidth = 1, colour = "black")
)
saveRDS(mrd.theme, "Output/mrd.fig.theme.RDS")

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

efx <- read.delim(
  file = "data/r.input_CRISPRGeneEffect_23Q2.csv",
  sep = ",",
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  file = "data/r.input_metadata_Model.csv",
  sep = ",",
  row.names = 1
)

################################################################################
# CLEAN DATA

# rewrite column names of efx data to be Gene ID only (remove numbers trailing)
colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))

# filter metadata for only those with effect score data
meta.filt <- meta[match(rownames(efx), rownames(meta)), ]

################################################################################
# FIND LINEAGE OF INTEREST

names(table(meta.filt$OncotreePrimaryDisease))
names(table(meta.filt$OncotreeLineage))
names(table(meta.filt$OncotreeSubtype))

"Glioblastoma" %in% unique(meta.filt$OncotreeSubtype)

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
  file = "Output/Rdata/01_essential.genes/00_depmap.gbm.fit_2024.01.31.RDS"
)

write.table(
  gbm.stats,
  file = "Output/Rdata/01_essential.genes/00_depmap.gbm.stats_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

# qq plot for quality control
pdf("Output/Figures/01_essential.genes/00_GBM_QQ_2024.07.08.pdf", height = 8, width = 8)
qqt(gbm.fit$t, df = gbm.fit$df.prior+gbm.fit$df.residual, pch = 16, cex = 0.2)
abline(0, 1)
dev.off()

# save metadata of n = 49 GBM cell lines
gbm.meta <- meta.filt[which(meta.filt$OncotreeSubtype == "Glioblastoma"), ]
write.table(
  gbm.meta,
  file = "Output/Rdata/01_essential.genes/00_gbm.meta.data_2024.05.01.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

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
  file = "Output/Rdata/01_essential.genes/01_depmap.gbm.mean.effect.scores_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

################################################################################
# SAVE GENE LISTS

# add mean effect score to gbm.stats output
gbm.stats_m <- merge(gbm.stats, gbm.mtx, by = "Gene")

# true essential genes (mean effect score <-0.5, lower in gbm lines)

lower <- gbm.stats_m[which(gbm.stats_m$EffectSize < 0), ] # all lower 8930
sig.lower <- lower[which(lower$q.value < 5e-2), ] # all sig lower 377
killing <- sig.lower[which(sig.lower$mEffectScore <= -0.5), ] # all sig lower w/ killing effect 168

# save all gene lists

write.table(
  lower,
  file = "Output/Rdata/01_essential.genes/03_essential.genes_all.lower_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

write.table(
  sig.lower,
  file = "Output/Rdata/01_essential.genes/04_essential.genes_sig.lower_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

write.table(
  killing,
  file = "Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)


################################################################################
# GENERATE RELEVANT FIGURES

# SET UP ----

# make color palette for enhanced volcano
colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(4) # has to be 4 by default for EnhancedVolcano (i.e. multiple axis cutoffs)
# "#67001F" "#F7B698" "#A7CFE4" "#053061"
volc.colors <- c("#F7B698", "#67001F", "#A7CFE4", "#053061")

# PLOTS ----

# volcano plot of limma results
## 05.01.2024 -- boxed labels, draw connectors not working for some reason...
pdf("Output/Figures/01_essential.genes/01_GBM_Volcano.Essential.Genes_small_2023.05.31.pdf", height = 4, width = 4)
EnhancedVolcano(
  gbm.stats,
  lab = gbm.stats$Gene, x = "EffectSize", y = "q.value",
  selectLab = killing$Gene,
  ylim = c(0, -log10(range(gbm.stats$q.value)[1])),
  xlim = c(range(gbm.stats$EffectSize)[1], range(gbm.stats$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  pCutoff = 5e-2, FCcutoff = 0, cutoffLineWidth = 0.5,
  pointSize = 2,
  labSize = 4, labCol = "black", labFace = "bold", # boxedLabels = TRUE,
  col = volc.colors,
  # drawConnectors = TRUE, widthConnectors = 0.5, typeConnectors = "closed",
  max.overlaps = 10,
  # maxoverlapsConnectors = 10
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/01_GBM_Volcano.Essential.Genes_big_2023.05.31.pdf", height = 8, width = 8)
EnhancedVolcano(
  gbm.stats,
  lab = gbm.stats$Gene, x = "EffectSize", y = "q.value",
  selectLab = killing$Gene,
  ylim = c(0, -log10(range(gbm.stats$q.value)[1])),
  xlim = c(range(gbm.stats$EffectSize)[1], range(gbm.stats$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  pCutoff = 5e-2, FCcutoff = 0, cutoffLineWidth = 0.5,
  pointSize = 4,
  labSize = 4, labCol = "black", labFace = "bold", # boxedLabels = TRUE,
  col = volc.colors,
  # drawConnectors = TRUE, widthConnectors = 0.5, typeConnectors = "closed",
  max.overlaps = 20,
  # maxoverlapsConnectors = 10
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()

# without labels
pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_no.labels_2024.02.01.pdf", height = 4, width = 4)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_no.labels_big_2024.02.01.pdf", height = 8, width = 8)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 4,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()


# make keys to change coloring, shapes, etc.
keyvals <- ifelse(
  gbm.stats_m$EffectSize < 0 & gbm.stats_m$q.value < 0.05, 19, 4
)
names(keyvals)[keyvals == 19] <- "GBM Essential"
names(keyvals)[keyvals == 4] <- "Not GBM Essential"

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_change.shape_small_2024.02.01.pdf", height = 4, width = 6)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 3,
  shapeCustom = keyvals,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme +
  theme(
    legend.position = "right"
  )
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_change.shape_big_2024.02.01.pdf", height = 8, width = 10)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 6,
  shapeCustom = keyvals,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme +
  theme(
    legend.position = "right"
  )
dev.off()

# make keys to change coloring, etc.
keyvals.2 <- ifelse(
  gbm.stats_m$EffectSize < 0 & gbm.stats_m$q.value < 0.05, "#A7CFE4",
  ifelse(gbm.stats_m$EffectSize < 0 & gbm.stats_m$q.value > 0.05, "#F7B698",
         ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$q.value > 0.05, "#67001F",
                ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$q.value < 0.05, "#053061", "black")
         )
  )
)
names(keyvals.2)[keyvals.2 == "#A7CFE4"] <- "GBM Essential"
names(keyvals.2)[keyvals.2 == "#F7B698"] <- "ES < 0, NS"
names(keyvals.2)[keyvals.2 == "#67001F"] <- "ES > 0, NS"
names(keyvals.2)[keyvals.2 == "#053061"] <- "ES > 0, S"
names(keyvals.2)[keyvals.2 == "black"] <- "NA"

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring.LEGEND_small_2024.02.01.pdf", height = 4, width = 6)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 3,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme +
  theme(
    legend.position = "right"
  )
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring.LEGEND_big_2024.02.01.pdf", height = 8, width = 10)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 6,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme +
  theme(
    legend.position = "right"
  )
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring_small_2024.02.01.pdf", height = 4, width = 4)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 3,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring_big_2024.02.01.pdf", height = 8, width = 8)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 6,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  coord_flip() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring_normal.axes_small_2024.02.01.pdf", height = 4, width = 4)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 3,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/02_GBM_Volcano.Essential.Genes_better.coloring_normal.axes_big_2024.02.01.pdf", height = 8, width = 8)
EnhancedVolcano(
  gbm.stats_m,
  lab = gbm.stats_m$Gene,
  labSize = 0,
  x = "EffectSize",
  y = "q.value",
  cutoffLineWidth = 0.5,
  pCutoff = 5e-2,
  FCcutoff = 0,
  pointSize = 6,
  colCustom = keyvals.2,
  ylim = c(0, -log10(range(gbm.stats_m$q.value)[1])),
  xlim = c(range(gbm.stats_m$EffectSize)[1], range(gbm.stats_m$EffectSize)[2]),
  xlab = "Effect Size",
  ylab = bquote(~-log[10] ~italic((p-value))),
  col = volc.colors
) +
  theme_minimal() +
  mrd.theme
dev.off()

# PLOTS ----
# updated volcano, easier customizaiton, vertical

gbm.stats_m <- gbm.stats_m %>%
  mutate(
    Coloring = ifelse(
      gbm.stats_m$EffectSize < 0 & gbm.stats_m$q.value < 0.05, "GBM Essential",
      ifelse(gbm.stats_m$EffectSize < 0 & gbm.stats_m$q.value > 0.05, "ES < 0, ns",
             ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$q.value > 0.05, "ES > 0, ns",
                    ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$q.value < 0.05, "ES > 0, s", "NA")
             )
      )
    )
  )
table(gbm.stats_m$Coloring)
plot.cols <- c("#F7B698", "#67001F", "#053061", "#A7CFE4")

pdf("Output/Figures/01_essential.genes/03_GBM_Volcano.Essential.Genes_vert_small_2024.02.01.pdf", height = 4, width = 4)
gbm.stats_m %>%
  ggplot(
    aes(x = EffectSize, y = -log10(q.value), fill = Coloring)
  ) +
  geom_point(size = 3, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols) +
  theme_bw() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/03_GBM_Volcano.Essential.Genes_vert_big_2024.02.01.pdf", height = 8, width = 8)
gbm.stats_m %>%
  ggplot(
    aes(x = EffectSize, y = -log10(q.value), fill = Coloring)
  ) +
  geom_point(size = 6, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols) +
  theme_bw() +
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/03_GBM_Volcano.Essential.Genes_vert_LEGEND_small_2024.02.01.pdf", height = 4, width = 6)
gbm.stats_m %>%
  ggplot(
    aes(x = EffectSize, y = -log10(q.value), fill = Coloring)
  ) +
  geom_point(size = 3, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols) +
  theme_bw() +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/01_essential.genes/03_GBM_Volcano.Essential.Genes_vert_LEGEND_big_2024.02.01.pdf", height = 8, width = 10)
gbm.stats_m %>%
  ggplot(
    aes(x = EffectSize, y = -log10(q.value), fill = Coloring)
  ) +
  geom_point(size = 6, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols) +
  theme_bw() +
  mrd.theme +
  theme(legend.position = "right")
dev.off()

# plot of mean effect scores vs. Effect Size

gbm.stats_m <- gbm.stats_m %>%
  mutate(
    Coloring.2 = ifelse(
      gbm.stats_m$EffectSize < 0 & gbm.stats_m$mEffectScore < -0.5, "GBM Essential",
      ifelse(gbm.stats_m$EffectSize < 0 & gbm.stats_m$mEffectScore > -0.5, "ES < 0, nk",
             ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$mEffectScore > -0.5, "ES > 0, nk",
                    ifelse(gbm.stats_m$EffectSize > 0 & gbm.stats_m$mEffectScore < -0.5, "ES > 0, k", "NA")
             )
      )
    )
  )
table(gbm.stats_m$Coloring.2)
# colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(5)
# "#67001F" "#E58267" "#F7F7F7" "#6AACD0" "#053061"
plot.cols2 <- c("#A7CFE4", "#053061", "#053061", "#F7F7F7")

pdf("Output/Figures/01_essential.genes/04_GBM_Volcano.EfxScore_invert_small_2024.02.01.pdf", height = 4, width = 4)
gbm.stats_m[which(gbm.stats_m$q.value < 0.05), ] %>%
  ggplot(
    aes(x = EffectSize, y = mEffectScore, fill = Coloring.2)
  ) +
  geom_point(size = 3, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols2) +
  theme_bw() + # theme_classic()
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/04_GBM_Volcano.EfxScore_invert_big_2024.02.01.pdf", height = 8, width = 8)
gbm.stats_m[which(gbm.stats_m$q.value < 0.05), ] %>%
  ggplot(
    aes(x = EffectSize, y = mEffectScore, fill = Coloring.2)
  ) +
  geom_point(size = 6, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols2) +
  theme_bw() + # theme_classic()
  mrd.theme
dev.off()

pdf("Output/Figures/01_essential.genes/04_GBM_Volcano.EfxScore_invert_LEGEND_small_024.02.01.pdf", height = 4, width = 6)
gbm.stats_m[which(gbm.stats_m$q.value < 0.05), ] %>%
  ggplot(
    aes(x = EffectSize, y = mEffectScore, fill = Coloring.2)
  ) +
  geom_point(size = 3, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols2) +
  theme_bw() + # theme_classic()
  mrd.theme +
  theme(legend.position = "right")
dev.off()

pdf("Output/Figures/01_essential.genes/04_GBM_Volcano.EfxScore_invert_LEGEND_big_024.02.01.pdf", height = 8, width = 10)
gbm.stats_m[which(gbm.stats_m$q.value < 0.05), ] %>%
  ggplot(
    aes(x = EffectSize, y = mEffectScore, fill = Coloring.2)
  ) +
  geom_point(size = 6, alpha = 0.8, shape = 21) +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = plot.cols2) +
  theme_bw() + # theme_classic()
  mrd.theme +
  theme(legend.position = "right")
dev.off()


write.table(
  gbm.stats_m,
  file = "Output/Rdata/01_essential.genes/02_depmap.gbm.stats_mean.effect.scores_2024.01.31.csv",
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)

################################################################################
# ESSENTIAL GENE ENRICHMENT ANALYSIS ####

dbs <- listEnrichrDbs()
dbs <- c(
  "GO_Biological_Process_2021",
  "GO_Molecular_Function_2021",
  "KEGG_2021_Human",
  "MSigDB_Hallmark_2020",
  "WikiPathway_2021_Human"
)

enriched <- enrichr(killing$Gene, dbs)

for (j in 1:length(enriched)) {

  pdf(paste0("Output/Figures/01_essential.genes/enrichR/p0", j, "_essential.genes_enrichr.Ratio_", dbs[j], "_2024.02.01.pdf"), height = 8, width = 10)
  print(
    plotEnrich(
      enriched[[j]],
      showTerms = 10,
      numChar = 40,
      y = "Ratio",
      orderBy = "P.value",
      title = ""
    ) +
      theme_classic() +
      mrd.theme +
      theme(
        legend.position = "right",
        title = element_text(size = 12, colour = "black")
      )
  )
  dev.off()

  pdf(paste0("Output/Figures/01_essential.genes/enrichR/p0", j, "_essential.genes_enrichr.Count_", dbs[j], "_2024.02.01.pdf"), height = 8, width = 10)
  print(
    plotEnrich(
      enriched[[j]],
      showTerms = 10,
      numChar = 40,
      y = "Count",
      orderBy = "P.value",
      title = ""
    ) +
      theme_classic() +
      mrd.theme +
      theme(
        legend.position = "right",
        title = element_text(size = 12, colour = "black")
      )
  )
  dev.off()

  pdf(paste0("Output/Figures/01_essential.genes/enrichR/p0", j, "_essential.genes_enrichr.Count_50.terms_", dbs[j], "_2024.02.01.pdf"), height = 8, width = 10)
  print(
    plotEnrich(
      enriched[[j]][which(enriched[[j]]$Adjusted.P.value < 0.05), ],
      showTerms = 50,
      numChar = 40,
      y = "Count",
      orderBy = "P.value",
      title = ""
    ) +
      theme_classic() +
      mrd.theme +
      theme(
        legend.position = "right",
        title = element_text(size = 12, colour = "black")
      )
  )
  dev.off()

  pdf(paste0("Output/Figures/01_essential.genes/enrichR/p0", j, "_essential.genes_enrichr.Ratio_50.terms_", dbs[j], "_2024.02.01.pdf"), height = 8, width = 10)
  print(
    plotEnrich(
      enriched[[j]][which(enriched[[j]]$Adjusted.P.value < 0.05), ],
      showTerms = 50, # max(ish) to get all with sig p-value
      numChar = 40,
      y = "Ratio",
      orderBy = "P.value",
      title = ""
    ) +
      geom_text(aes(label = Overlap), hjust = -0.2, colour = "black") +
      theme_classic() +
      mrd.theme +
      theme(
        legend.position = "right",
        title = element_text(size = 12, colour = "black")
      )
  )
  dev.off()

}

saveRDS(
  enriched,
  file = "Output/Rdata/01_essential.genes/06_essential.gene_enrichment.results_2024.02.05.RDS"
)

################################################################################
# ESSENTIAL GENE ENRICHMENT ANALYSIS, UPGRADED ####

dir.create("Output/Figures/01_essential.genes/pan_enrich")

killing <- killing %>%
  arrange(mEffectScore)
glimpse(killing)

gs.list <- list(
  ess = killing$Gene
)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")

msigdb_ref <- msigdb %>%
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

glimpse(comp@compareClusterResult)

funni <- function(X) {
  db <- str_split(X, pattern = "_")[[1]][1]
  return(db)
}

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

pdf("Output/Figures/01_essential.genes/pan_enrich/01_msigdb.C2_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Count",
  color = "p.adjust",
  # split = "DB",
  showCategory = 25,
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

# TRANSCRIPTION FACTOR ENRICHMENT ANALYSIS

msigdb <- msigdbr(species = "Homo sapiens", category = "C3")

msigdb_ref <- msigdb %>%
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

pdf("Output/Figures/01_essential.genes/pan_enrich/02_msigdb.C3.TF_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Count",
  color = "p.adjust",
  # split = "DB",
  showCategory = 25,
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

# COMPUTATIONAL CANCER

msigdb <- msigdbr(species = "Homo sapiens", category = "C4")

msigdb_ref <- msigdb %>%
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

pdf("Output/Figures/01_essential.genes/pan_enrich/03_msigdb.C4_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Count",
  color = "p.adjust",
  # split = "DB",
  showCategory = 25,
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

# # ONCOGENIC GENE SETS # NO ENRICHMENT FOUND #
#
# msigdb <- msigdbr(species = "Homo sapiens", category = "C6")
#
# msigdb_ref <- msigdb %>%
#   distinct(gs_name, gene_symbol) %>%
#   as.data.frame()
#
# comp <- compareCluster(
#   gs.list,
#   fun = "enricher",
#   TERM2GENE = msigdb_ref,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
# )
#
# pdf("Output/Figures/01_essential.genes/pan_enrich/04_msigdb.C6.ONC_ORA_2024.03.13.pdf", height = 8, width = 10)
# dotplot(
#   comp,
#   x = "Count",
#   color = "p.adjust",
#   # split = "DB",
#   showCategory = 25,
#   label_format = 50
# ) +
#   # facet_grid("DB") +
#   theme_bw() +
#   # geom_text(aes(label = GeneRatio), hjust = -0.2, colour = "black") +
#   mrd.theme +
#   theme(
#     legend.position = "right",
#     legend.title = element_text(size = 11, colour = "black")
#   )
# dev.off()

# FILTERED, FINAL, POLISHED

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

################################################################################
# ESSENTIAL GENE DENSITY PLOTS ####

dir.create("Output/Figures/01_essential.genes/density")

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
