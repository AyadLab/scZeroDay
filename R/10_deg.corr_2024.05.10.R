library(tidyverse)
library(ggpubr)
library(corrplot)
library(cowplot)
library(dittoSeq)
library(RColorBrewer)
library(paletteer)
library(VennDiagram)

dir.create("Output/Figures/10_deg.corr")
dir.create("Output/Rdata/10_deg.corr")

################################################################################
# READ IN DEG

nef.deg <- read.csv(
  "Output/Rdata/04_states.2_2024.05.02/05_mVC_DEGs_2023.10.17.csv",
  row.names = 1
)
nef.sig <- nef.deg[which(nef.deg$p_val_adj < 0.05), ]
sut.deg <- read.csv(
  "Output/Rdata/07_suter_2024.05.06/02_mVC_DEGs_2024.02.26.csv",
  row.names = 1
)
sut.sig <- sut.deg[which(sut.deg$p_val_adj < 0.05), ]
jon.deg <- read.csv(
  "Output/Rdata/08_johnson_2024.05.06/02_mVC_DEGs_2024.02.26.csv",
  row.names = 1
)
jon.sig <- jon.deg[which(jon.deg$p_val_adj < 0.05), ]

################################################################################
# VS1 VENN

nef.v1 <- nef.sig$gene[which(nef.sig$cluster == "mVC1")]
sut.v1 <- sut.sig$gene[which(sut.sig$cluster == "mVC1")]
jon.v1 <- jon.sig$gene[which(jon.sig$cluster == "mVC1")]
nVs <- intersect(nef.v1, sut.v1)
nVj <- intersect(nef.v1, jon.v1)
sVj <- intersect(sut.v1, jon.v1)
nsVnj <- intersect(nVs, nVj)
njVsj <- intersect(nVj, sVj)
n123 <- intersect(nsVnj, njVsj)

pdf("Output/Figures/10_deg.corr/01_VS1.venn_2024.05.10.pdf", height = 4, width = 4)
draw.triple.venn(
  area1 = length(nef.v1),
  area2 = length(sut.v1),
  area3 = length(jon.v1),
  n12 = length(intersect(nef.v1, sut.v1)),
  n13 = length(intersect(nef.v1, jon.v1)),
  n23 = length(intersect(sut.v1, jon.v1)),
  n123 = length(n123),
  category = c("Neftel VS1", "Suter VS1", "Johnson VS1"),
  print.mode = c("raw", "percent"),
  fill = colorRampPalette(paletteer_d(palette = "RColorBrewer::Pastel2", 8))(3),
  alpha = rep(0.8, 3),
  cat.pos = c(0, 0, 180)
)
dev.off()

################################################################################
# VS2 VENN

nef.v2 <- nef.sig$gene[which(nef.sig$cluster == "mVC2")]
sut.v2 <- sut.sig$gene[which(sut.sig$cluster == "mVC2")]
jon.v2 <- jon.sig$gene[which(jon.sig$cluster == "mVC2")]
nVs <- intersect(nef.v2, sut.v2)
nVj <- intersect(nef.v2, jon.v2)
sVj <- intersect(sut.v2, jon.v2)
nsVnj <- intersect(nVs, nVj)
njVsj <- intersect(nVj, sVj)
n123 <- intersect(nsVnj, njVsj)

pdf("Output/Figures/10_deg.corr/02_VS2.venn_2024.05.10.pdf", height = 4, width = 4)
draw.triple.venn(
  area1 = length(nef.v2),
  area2 = length(sut.v2),
  area3 = length(jon.v2),
  n12 = length(intersect(nef.v2, sut.v2)),
  n13 = length(intersect(nef.v2, jon.v2)),
  n23 = length(intersect(sut.v2, jon.v2)),
  n123 = length(n123),
  category = c("Neftel VS2", "Suter VS2", "Johnson VS2"),
  print.mode = c("raw", "percent"),
  fill = colorRampPalette(paletteer_d(palette = "RColorBrewer::Pastel2", 8))(3),
  alpha = rep(0.8, 3),
  cat.pos = c(0, 0, 180)
)
dev.off()

################################################################################
# VS3 VENN

nef.v3 <- nef.sig$gene[which(nef.sig$cluster == "mVC3")]
sut.v3 <- sut.sig$gene[which(sut.sig$cluster == "mVC3")]
jon.v3 <- jon.sig$gene[which(jon.sig$cluster == "mVC3")]
nVs <- intersect(nef.v3, sut.v3)
nVj <- intersect(nef.v3, jon.v3)
sVj <- intersect(sut.v3, jon.v3)
nsVnj <- intersect(nVs, nVj)
njVsj <- intersect(nVj, sVj)
n123 <- intersect(nsVnj, njVsj)

pdf("Output/Figures/10_deg.corr/03_VS3.venn_2024.05.10.pdf", height = 4, width = 4)
draw.triple.venn(
  area1 = length(nef.v3),
  area2 = length(sut.v3),
  area3 = length(jon.v3),
  n12 = length(intersect(nef.v3, sut.v3)),
  n13 = length(intersect(nef.v3, jon.v3)),
  n23 = length(intersect(sut.v3, jon.v3)),
  n123 = length(n123),
  category = c("Neftel VS3", "Suter VS3", "Johnson VS3"),
  print.mode = c("raw", "percent"),
  fill = colorRampPalette(paletteer_d(palette = "RColorBrewer::Pastel2", 8))(3),
  alpha = rep(0.8, 3),
  cat.pos = c(0, 0, 180)
)
dev.off()
