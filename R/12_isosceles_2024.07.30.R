# devtools::install_local("/data/mrd/isosceles/ISOSCELES_0.0.0.9000.tar.gz")
#
# library(shiny)
# library(ISOSCELES) # load the ISOSCELES library including necessary drug signature data
# ISOSCELES::runISOSCELES() # launch the ISOSCELES shiny GUI
#
# # have to run off-quorra bc app data pulls from local

################################################################################

library(tidyverse)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
# library(FSA)
library(limma)


# dir.create("Output/Rdata/12_isosceles")
# dir.create("Output/Figures/12_isosceles")

################################################################################
# READ IN DATA
obj <- readRDS("data/Suter_scGBM_2023.10.04.RDS")

corr.mtx <- read.csv(
  "Output/Rdata/12_isosceles/RDS_upload_L1000_consensus_corrmat.csv",
  row.names = 1
)
corr.mtx <- as.data.frame(t(corr.mtx))
colnames(corr.mtx)

# rep.hub <- read.table(
#   "Output/ISOSCELES/s3.amazonaws.com_data.clue.io_repurposing_downloads_repurposing_drugs_20200324.txt",
#   sep = "\t",
#   skip = 9,
#   header = TRUE,
#   quote = ""
# )
# rownames(rep.hub) <- rep.hub$pert_iname

rep.hub <- read.csv(
  "data/RepHub_Annotations_2023.10.27.csv",
  row.names = 1
)

################################################################################
# SUBSET CORRELATIION MATRIX

corr.mtx.ann <- corr.mtx %>%
  mutate(
    VS = case_when(
      rownames(corr.mtx) %in% rownames(obj@meta.data) ~ paste0(obj@meta.data$mVC)
    )
  )
table(corr.mtx.ann$VS)

VS.mtx <- subset(
  corr.mtx.ann,
  corr.mtx.ann$VS %in% c("mVC1", "mVC2", "mVC3")
)

################################################################################
# PROCESSING

avg.mtx <- VS.mtx %>%
  group_by(VS) %>%
  summarise_all(.funs = mean) %>%
  as.matrix()

# make dataframe (rows = drugs, cols = VS)
drugs <- as.data.frame(t(avg.mtx))
colnames(drugs) <- drugs[1, ]
drugs <- drugs[-1, ]
drugs2 <- as.data.frame(lapply(drugs, as.double))
rownames(drugs2) <- rownames(drugs)

# change formatting for identification of top drugs by VS
top.drugs.pre <- drugs2
top.drugs.pre$drug <- rownames(top.drugs.pre)

top.drugs <- reshape2::melt(top.drugs.pre, id = "drug")
colnames(top.drugs) <- c("drug", "VS", "correlation")
range(top.drugs$correlation)

################################################################################
# SUMMARIZE BY VULNERABILITY CLUSTER

# top 25
top.drugs.25 <- top.drugs %>%
  group_by(VS) %>%
  slice_min(n = 25, correlation)
glimpse(top.drugs.25)

top.25 <- subset(drugs2, rownames(drugs2) %in% top.drugs.25$drug)
top.25.inv <- top.25*-1

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.25.inv), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.25.inv)/paletteLength, max(top.25.inv), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/01_Drugs.mVC_top25.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.25.inv,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO
)
dev.off()

# top 50
top.drugs.50 <- top.drugs %>%
  group_by(VS) %>%
  slice_min(n = 50, correlation)
glimpse(top.drugs.50)

top.50 <- subset(drugs2, rownames(drugs2) %in% top.drugs.50$drug)
top.50.inv <- top.50*-1

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.50.inv), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.50.inv)/paletteLength, max(top.50.inv), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/02_Drugs.mVC_top50.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.50.inv,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdYlBu"))(100)),
  breaks = myBreaksISO, fontsize = 7
)
dev.off()

# scaled top 100
top.drugs.100 <- top.drugs %>%
  group_by(VS) %>%
  slice_min(n = 100, correlation)
glimpse(top.drugs.100)

top.100 <- subset(drugs2, rownames(drugs2) %in% top.drugs.100$drug)
top.100.inv <- top.100*-1
top.100.scale <- scale(top.100.inv)

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.100.scale), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.100.scale)/paletteLength, max(top.100.scale), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/03_Drugs.mVC_top100.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.100.scale,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO, fontsize = 3
)
dev.off()


################################################################################
# ANALYSIS BASED ON DRUG REPURPOSING HUB

# # which drugs missing from rep.hub (at least to make fig)
# rownames(top.100)[which(!rownames(top.100) %in% rep.hub$new_iname)]

# # add missing to repo hub df
# # XMD.16144
# rep.hub$new_iname[grep(".16144", rep.hub$new_iname)] # name in rep.hub AZD8330
# rownames(top.10)[grep("XMD.", rownames(top.10))] # look in drugs
# # add to rep.hub
# rep.hub <- rbind(rep.hub, c("XMD-16144", "Unknown", "kinase inhibitor", "", "", "", "XMD.16144"))
# # fix rowname
# rownames(rep.hub)[grep("XMD.16144", rep.hub$new_iname)] <- "XMD-16144"

# # add moa annotations
# row.ann <- as.data.frame(
#   cbind(
#     drugs = rownames(top.100),
#     MOA = rep.hub$moa[which(rep.hub$new_iname %in% rownames(top.100))]
#   )
# )
# rownames(row.ann) <- row.ann$drugs
# row.ann$drugs <- NULL

# # plot w/ annotations
# pdf("Output/ISOSCELES/05_Drugs.mVC_top10.INVERSE_2023.10.31.pdf", height = 10.5, width = 14)
# pheatmap(
#   top.10.inv,
#   annotation_row = row.ann
# )
# dev.off()

# may be missing some drugs if name in l1000 != exact name in rep.hub!!
# which drugs are "Launched"
repo.drugs <- top.drugs[which(top.drugs$drug %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched")]), ]

# top 100
top.repo.100 <- repo.drugs %>%
  group_by(VS) %>%
  slice_min(n = 100, correlation)
glimpse(top.repo.100)

top.100.rep <- subset(drugs2, rownames(drugs2) %in% top.repo.100$drug)
top.100.rep.inv <- top.100.rep*-1
top.100.rep.scale <- scale(top.100.rep.inv)

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.100.rep.scale), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.100.rep.scale)/paletteLength, max(top.100.rep.scale), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/04_Drugs.mVC_top100.RepHub_INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.100.rep.scale,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO, fontsize = 4
)
dev.off()

# which drugs are "Launched" & "Oncology"
repo.drugs <- top.drugs[which(top.drugs$drug %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]

# top 10
top.repo.10 <- repo.drugs %>%
  group_by(VS) %>%
  slice_min(n = 20, correlation)
glimpse(top.repo.10)

top.10.rep <- subset(drugs2, rownames(drugs2) %in% top.repo.10$drug)
top.10.rep.inv <- top.10.rep*-1
top.10.rep.scale <- scale(top.10.rep.inv)

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.10.rep.scale), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.10.rep.scale)/paletteLength, max(top.10.rep.scale), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/05_Drugs.mVC_top20.RepHub.Onc_INVERSE_2023.10.31.pdf", height = 10, width = 8)
pheatmap(
  top.10.rep.scale,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO, fontsize = 8
)
dev.off()

# ################################################################################
# # DIFFERENTIAL CORRELATION ANALYSIS
#
# limma

# # VS1/2/3 v. Non-Neoplastic
#
# # mtx <- corr.mtx.ann[,-c(length(colnames(corr.mtx.ann)))]
# # design <- data.frame(
# #   rownames = row.names(corr.mtx.ann),
# #   VS = corr.mtx.ann$VS
# # )
# # levels(design)
# # design$VS <- as.factor(design$VS)
# # design$VS <- relevel(design$VS, ref="Non-Neoplastic")
# #
# # dsgn <- model.matrix(~design$VS)
# #
# #
# # fit <- lmFit(t(mtx), design = dsgn)
# # fit <- eBayes(fit)
# # results <- topTable(fit, coef = "design$VSmVC1", number = Inf)
# # glimpse(results)
# # sig <- results[which(results$adj.P.Val < 0.05), ]
#
# # VS1 v. VS2/3
#
# # VS2 v. VS1/3
#
# # VS3 v. VS1/2
#
# all pairwise

# mtx <- VS.mtx[,-c(length(colnames(VS.mtx)))]
# mtx <- mtx*-1
# dgn <- data.frame(
#   rownames = row.names(VS.mtx),
#   VS = VS.mtx$VS
# )
# design <- model.matrix(~0+dgn$VS)
# colnames(design) <- c("VS1", "VS2", "VS3")
# design
#
# contrasts <- makeContrasts(
#   VS1-VS2, VS1-VS3, VS2-VS3,
#   levels=colnames(design)
# )
# contrasts
#
# fit <- lmFit(t(mtx), design = design)
# fit <- eBayes(fit)
#
# fit2 <- contrasts.fit(fit, contrasts = contrasts)
# fit2 <- eBayes(fit2)
#
# results_1.2 <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", p.value = 0.05)
# results_1.3 <- topTable(fit2, coef = 2, number = Inf, adjust.method = "BH", p.value = 0.05)
# results_2.3 <- topTable(fit2, coef = 3, number = Inf, adjust.method = "BH", p.value = 0.05)
#
# sig.1.2 <- results_1.2[which(results_1.2$logFC > 0), ]
# sig.1.3 <- results_1.3[which(results_1.3$logFC > 0), ]
# sig.1.2_25 <- sig.1.2 %>%
#   slice_max(n = 25, order_by = logFC)
# sig.1.3_25 <- sig.1.3 %>%
#   slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.1.2_25)), unique(rownames(sig.1.3_25))))
# # [1] "amoxapine"    "berberine"    "cevimeline"   "cholic.acid"  "D.serine"     "emedastine"
# # [7] "fenobam"      "GDC.0834"     "imiloxan"     "istaroxime"   "ketorolac"    "lodoxamide"
# # [13] "medetomidine" "TH.302"
#
#
# sig.2.1 <- results_1.2[which(results_1.2$logFC < 0), ]
# sig.2.3 <- results_2.3[which(results_2.3$logFC > 0), ]
# sig.2.1_25 <- sig.2.1 %>%
#   slice_min(n = 25, order_by = logFC)
# sig.2.3_25 <- sig.2.3 %>%
#   slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.2.1_25)), unique(rownames(sig.2.3_25))))
# # ?
#
#
# sig.3.1 <- results_1.3[which(results_1.3$logFC < 0), ]
# sig.3.2 <- results_2.3[which(results_2.3$logFC < 0), ]
# sig.3.1_25 <- sig.3.1 %>%
#   slice_min(n = 25, order_by = logFC)
# sig.3.2_25 <- sig.3.2 %>%
#   slice_min(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.3.1_25)), unique(rownames(sig.3.2_25))))
# # [1] "A.443654"      "alvocidib"     "BMS.387032"    "bortezomib"    "carfilzomib"
# # [6] "CGP.60474"     "defactinib"    "dinaciclib"    "lacidipine"    "MG.132"
# # [11] "mocetinostat"  "NVP.BEZ235"    "pentobarbital" "perampanel"    "torin.2"
# # [16] "triptolide"    "WZ.3105"       "YM.155"


# # # run ANOVA between VS across all compounds
# # anova <- list()
# # p <- c()
# # for (i in 1:(length(colnames(VS.mtx))-1)) {
# #   res <- aov(VS.mtx[,i] ~ VS, data = VS.mtx)
# #
# #   anova[[i]] <- res
# #   names(anova)[i] <- paste0(colnames(VS.mtx)[i])
# #
# #   p[i] <- summary(res)[[1]]$`Pr(>F)`[1]
# #   names(p)[i] <- paste0(colnames(VS.mtx)[i])
# # }
# #
# # tukey <- list()
# # for (i in 1:length(anova)) {
# #   tuk <- TukeyHSD(anova[[i]])
# #   tukey[[i]] <- tuk
# #   names(tukey)[i] <- names(anova)[i]
# # }
# #
# # # subset for only compounds with significant differences
# # sig.drugs <- p[which(p < 0.05)]
# # # 1659 / 1674 drugs significantly different
# # sig.tukey <- tukey[which(names(tukey) %in% names(sig.drugs))]
# # print(sig.tukey[[1]])
# # # diff         lwr         upr p adj
# # # mVC2-mVC1 0.009354567 0.008073580 0.010635555     0
# # # mVC3-mVC1 0.015676495 0.014717397 0.016635593     0
# # # mVC3-mVC2 0.006321928 0.005003519 0.007640337     0


#
# # kruskal-wallis test with dunn post-hoc
#
# glimpse(top.drugs)
#
# kw <- top.drugs %>%
#   group_by(drug) %>%
#   summarise(p.value = aov(correlation ~ VS)$p.value)
# glimpse(kw)
#
# # Filter for significant drugs
# sig <- kw %>%
#   filter(p.value < 0.05)  # Adjust threshold as needed
#
# # Perform Dunn's test for each significant drug
# dunn_results <- significant_drugs %>%
#   rowwise() %>%
#   mutate(dunn = list(dunnTest(SpearmanCorrelation ~ Group, data = data[data$Drug == Drug, ], method = "bonferroni"))) %>%
#   unnest(cols = c(dunn))



# applying linear fixed effects model -- just use limma
# glimpse(corr.mtx.ann)
#
# dat <- corr.mtx.ann
# dat$cell <- paste0(rownames(dat), "_", dat$VS)
# data <- reshape2::melt(dat, id = "cell")
#
# data$VS <- paste0(str_split_i(data$cell, "_", i = 3))
# colnames(data) <- c("cell", "drug", "correlation", "VS")
#
# library(lme4)
#
# # Fit the model for each drug
# results <- data %>%
#   group_by(drug) %>%
#   do(model = lm(correlation ~ VS + cell, data = .)) # cell included as extra fixed effect
# # results <- data %>%
# #   group_by(drug) %>%
# #   do(model = lm(correlation ~ VS, data = .)) # using fixed effect model
#
#
# library(broom.mixed)
#
# results <- results %>%
#   mutate(summary = map(model, tidy)) %>%
#   unnest(summary) %>%
#   filter(term == "Group") %>%
#   mutate(p.value = p.adjust(p.value, method = "bonferroni"))  # Adjust for multiple comparisons
#
#
# posthoc <- glht(model, linfct = mcp(Group = "Tukey"))
# summary(posthoc)

################################################################################

# limma implementation

# organize data matrix
mtx <- VS.mtx[,-c(length(colnames(VS.mtx)))]
mtx <- mtx*-1
dat.mtx <- as.matrix(t(mtx))

# create factor for groups
group <- factor(VS.mtx$VS)

# create design matrix
design <- model.matrix(~0 + group)  # ~0 removes the intercept, so each group is a separate column
colnames(design) <- levels(group)

fit <- lmFit(dat.mtx, design)

# Define contrasts to compare the groups, for example:
contrast_matrix <- makeContrasts(
  VS1_vs_VS2 = mVC1 - mVC2,
  VS1_vs_VS3 = mVC1 - mVC3,
  VS2_vs_VS3 = mVC2 - mVC3,
  levels = design
)

# Apply the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes moderation
fit2 <- eBayes(fit2)

saveRDS(fit2, "Output/Rdata/12_isosceles/limma_out_2024.11.08.RDS")

# Get the top table of results for each contrast
results_VS1.VS2 <- topTable(fit2, coef = "VS1_vs_VS2", number = Inf, adjust.method = "bonferroni", p.value = 0.05)
results_VS1.VS3 <- topTable(fit2, coef = "VS1_vs_VS3", number = Inf, adjust.method = "bonferroni", p.value = 0.05)
results_VS2.VS3 <- topTable(fit2, coef = "VS2_vs_VS3", number = Inf, adjust.method = "bonferroni", p.value = 0.05)

# VS1-specific drugs
sig.1.2 <- results_VS1.VS2[which(results_VS1.VS2$logFC > 0), ]
sig.1.3 <- results_VS1.VS3[which(results_VS1.VS3$logFC > 0), ]
sig.1.2_25 <- sig.1.2 %>%
  slice_max(n = 25, order_by = logFC)
sig.1.3_25 <- sig.1.3 %>%
  slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.1.2_25)), unique(rownames(sig.1.3_25))))

# VS2-specific drugs
sig.2.1 <- results_VS1.VS2[which(results_VS1.VS2$logFC < 0), ]
sig.2.3 <- results_VS2.VS3[which(results_VS2.VS3$logFC > 0), ]
sig.2.1_25 <- sig.2.1 %>%
  slice_min(n = 25, order_by = logFC)
sig.2.3_25 <- sig.2.3 %>%
  slice_max(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.2.1_25)), unique(rownames(sig.2.3_25))))

# VS3-specific drugs
sig.3.1 <- results_VS1.VS3[which(results_VS1.VS3$logFC < 0), ]
sig.3.2 <- results_VS2.VS3[which(results_VS2.VS3$logFC < 0), ]
sig.3.1_25 <- sig.3.1 %>%
  slice_min(n = 25, order_by = logFC)
sig.3.2_25 <- sig.3.2 %>%
  slice_min(n = 25, order_by = logFC)
# sort(intersect(unique(rownames(sig.3.1_25)), unique(rownames(sig.3.2_25))))

top.drugs <- c(rownames(sig.1.2_25), rownames(sig.1.3_25),
               rownames(sig.2.1_25), rownames(sig.2.3_25),
               rownames(sig.3.1_25), rownames(sig.3.2_25))
top.drugs <- unique(top.drugs) ####


res.1 <- results_VS1.VS2[which(rownames(results_VS1.VS2) %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]
res.2 <- results_VS1.VS3[which(rownames(results_VS1.VS3) %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]
res.3 <- results_VS2.VS3[which(rownames(results_VS2.VS3) %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]

launch.1.2 <- res.1[which(res.1$logFC > 0), ]
launch.1.3 <- res.2[which(res.2$logFC > 0), ]

launch.2.1 <- res.1[which(res.1$logFC < 0), ]
launch.2.3 <- res.3[which(res.3$logFC > 0), ]

launch.3.1 <- res.2[which(res.2$logFC < 0), ]
launch.3.2 <- res.3[which(res.3$logFC < 0), ]

launched <- c(rownames(launch.1.2), rownames(launch.1.3),
               rownames(launch.2.1), rownames(launch.2.3),
               rownames(launch.3.1), rownames(launch.3.2))
launched <- unique(launched) ####


avg.mtx <- VS.mtx %>%
  group_by(VS) %>%
  summarise_all(.funs = mean) %>%
  as.matrix()

# make dataframe (rows = drugs, cols = VS)
drugs <- as.data.frame(t(avg.mtx))
colnames(drugs) <- drugs[1, ]
drugs <- drugs[-1, ]
drugs2 <- as.data.frame(lapply(drugs, as.double))
rownames(drugs2) <- rownames(drugs)

top <- subset(drugs2, rownames(drugs2) %in% top.drugs)
top.inv <- top*-1
top.scale <- scale(top.inv)

lnch <- subset(drugs2, rownames(drugs2) %in% launched)
lnch.inv <- lnch*-1
lnch.scale <- scale(lnch.inv)

paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaksISO <- c(seq(min(top.scale), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.scale)/paletteLength, max(top.scale), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/07_VS_top.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.scale,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO
)
dev.off()

myBreaksISO <- c(seq(min(top.inv), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(top.inv)/paletteLength, max(top.inv), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/07_raw.corr_VS_top.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  top.inv,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO
)
dev.off()

pdf("Output/Figures/12_isosceles/07.1_VS.top.INVERSE_ann_scale.col_2023.10.31.pdf", height = 4, width = 12)
pheatmap(
  t(top.inv),
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  show_rownames = FALSE, angle_col = 45,
  scale = "column",
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  fontsize = 6
)
dev.off()

# for this one, I changed above to top 50 instead of top 25 manually, and reverted the code back to 25
pdf("Output/Figures/12_isosceles/07.2_VS.top50_INVERSE_ann_scale.col_2024.11.19.pdf", height = 4, width = 16)
pheatmap(
  t(top.inv),
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  show_rownames = FALSE, angle_col = 45,
  scale = "column",
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  fontsize = 4
)
dev.off()



myBreaksISO.lnch <- c(seq(min(lnch.scale), 0, length.out=ceiling(paletteLength/2) + 1),
                 seq(max(lnch.scale)/paletteLength, max(lnch.scale), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/08_VS.launched.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  lnch.scale,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO.lnch
)
dev.off()

ann <- data.frame(
  row.names = c("mVC1", "mVC2", "mVC3"),
  VS = c("VS1", "VS2", "VS3")
)

myBreaksISO.lnch <- c(seq(min(lnch.inv), 0, length.out=ceiling(paletteLength/2) + 1),
                      seq(max(lnch.inv)/paletteLength, max(lnch.inv), length.out=floor(paletteLength/2)))
pdf("Output/Figures/12_isosceles/08.1_VS.launched.INVERSE_ann_2023.10.31.pdf", height = 6, width = 12)
pheatmap(
  t(lnch.inv),
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  show_rownames = FALSE, angle_col = 45,
  # scale = "column",
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO.lnch
)
dev.off()
pdf("Output/Figures/12_isosceles/08.2_VS.launched.INVERSE_ann_2023.10.31.pdf", height = 6, width = 12)
pheatmap(
  t(lnch.inv),
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  show_rownames = FALSE, angle_col = 45,
  scale = "column",
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  # breaks = myBreaksISO.lnch
)
dev.off()



myBreaksISO.lnch <- c(seq(min(lnch.inv), 0, length.out=ceiling(paletteLength/2) + 1),
                      seq(max(lnch.inv)/paletteLength, max(lnch.inv), length.out=floor(paletteLength/2)))

pdf("Output/Figures/12_isosceles/08_raw.corr_VS.launched.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
pheatmap(
  lnch.inv,
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  breaks = myBreaksISO.lnch
)
dev.off()

# ################################################################################
# # ADD CORRELAIONS AS ASSAY OBJECT
#
# corr_assay <- CreateAssayObject(data = as.matrix(t(corr.mtx)*-1))
# obj[["ISOSCELES"]] <- corr_assay
# obj@assays$ISOSCELES$counts <- obj@assays$ISOSCELES$data
#
# Idents(obj) <- obj$mVC
# table(Idents(obj))
# neo <- WhichCells(obj, idents = "Non-Neoplastic", invert = TRUE)
# wil_limma <- FindAllMarkers(
#   subset(obj, cells = neo),
#   test.use = "wilcox",
#   only.pos = FALSE, logfc.threshold = -Inf,
#   assay = "ISOSCELES", slot = "data"
# )
# glimpse(wil_limma)
#
#
# top <- wil_limma[which(wil_limma$p_val_adj > 0.05), ] %>%
#   group_by(cluster) %>%
#   slice_max(n = 25, order_by = avg_log2FC)
# sort(top$gene[which(top$cluster == "mVC1")])
# sort(top$gene[which(top$cluster == "mVC2")])
# sort(top$gene[which(top$cluster == "mVC3")])
#
#
# launched <- wil_limma[which(wil_limma$gene %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]
# top.launched <- launched[which(launched$p_val_adj < 0.05), ] %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# sort(top.launched$gene[which(top.launched$cluster == "mVC1")])
# sort(top.launched$gene[which(top.launched$cluster == "mVC2")])
# sort(top.launched$gene[which(top.launched$cluster == "mVC3")])
#
#
# # obj <- ScaleData(obj, assay = "ISOSCELES")
# # DoHeatmap(
# #   obj,
# #   features = unique(top$gene[which(top$cluster != "Non-Neoplastic")]),
# #   cells = neo,
# #   group.by = "mVC",
# #   assay = "ISOSCELES",
# #   slot = "scale.data"
# # )
# # pseudo <- AggregateExpression(
# #   subset(obj, cells = neo),
# #   assays = "ISOSCELES", normalization.method =
# #   slot = "data",
# #   group.by = "ident",
# #   return.seurat = TRUE
# # )
# # pseudo.mtx <- pseudo@assays$ISOSCELES$data
#
# # plot <- as.matrix(t(corr.mtx)*-1)
# # plot <- plot[which(rownames(plot) %in% top$gene), ]
# # plot <- as.matrix(pseudo.mtx)[which(rownames(pseudo.mtx) %in% top$gene), ]
#
# # ann <- data.frame(
# #   row.names = rep.hub$new_iname[which(rep.hub$new_iname %in% top$gene)],
# #   moa = rep.hub$moa[which(rep.hub$new_iname %in% top$gene)]
# # )
# # lapply(ann, FUN = {})
# # col <- data.frame(
# #   row.names = rownames(obj@meta.data),
# #   vs = obj$mVC
# # )
#
# avg.mtx <- VS.mtx %>%
#   group_by(VS) %>%
#   summarise_all(.funs = mean) %>%
#   as.matrix()
#
# # make dataframe (rows = drugs, cols = VS)
# drugs <- as.data.frame(t(avg.mtx))
# colnames(drugs) <- drugs[1, ]
# drugs <- drugs[-1, ]
# drugs2 <- as.data.frame(lapply(drugs, as.double))
# rownames(drugs2) <- rownames(drugs)
#
# top.25 <- subset(drugs2, rownames(drugs2) %in% top$gene[which(top$cluster != "Non-Neoplastic")])
# top.25.inv <- top.25*-1
# top.25.scale <- scale(top.25.inv)
#
# paletteLength <- 100
# # use floor and ceiling to deal with even/odd length pallettelengths
# myBreaksISO <- c(seq(min(top.25.scale), 0, length.out=ceiling(paletteLength/2) + 1),
#                  seq(max(top.25.scale)/paletteLength, max(top.25.scale), length.out=floor(paletteLength/2)))
#
# pdf("Output/Figures/12_isosceles/06_DED.mVC_top25.INVERSE_2023.10.31.pdf", height = 10.5, width = 7)
# pheatmap(
#   top.25.scale,
#   color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
#   breaks = myBreaksISO
# )
# dev.off()

# ################################################################################

write.csv(
  rep.hub,
  "Output/Rdata/RepHub_Annotations_2023.10.27.csv"
)


################################################################################
# if you scale the matrix first, then look for the top drugs,
# you get more distinct results

scale.drugs <- as.data.frame(t(scale(t(drugs2), center = TRUE, scale = TRUE)))

ann.top <- scale.drugs
ann.top$drugs <- rownames(ann.top)

top.drugs <- reshape2::melt(ann.top)
glimpse(top.drugs)
colnames(top.drugs) <- c("drug", "mVC", "relative.corr")

# TOP 25 REGULONS PER CELL STATE
top.25 <- top.drugs %>%
  group_by(mVC) %>%
  slice_min(n = 25, order_by = relative.corr)
glimpse(top.25)

scale.10 <- subset(scale.drugs, rownames(scale.drugs) %in% unique(top.25$drug))

my.breaks.1 <- c(seq(min(scale.10), 0, length.out=ceiling(100/2) + 1),
                 seq(max(scale.10)/100, max(scale.10), length.out=floor(100/2)))

pdf("Output/Figures/12_isosceles/09_VS.drugs_scaled_2024.11.11.pdf", height = 6, width = 12)
pheatmap(
  t(scale.10)*-1,
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  show_rownames = FALSE,
  angle_col = "45", fontsize_row = 8, fontsize_col = 8,
  breaks = my.breaks.1
)
dev.off()

# TOP BY LAUNCHED
top.lnchd <- top.drugs[which(top.drugs$drug %in% rep.hub$new_iname[which(rep.hub$clinical_phase == "Launched" & rep.hub$disease_area == "oncology")]), ]
tippy <- top.lnchd %>%
  group_by(mVC) %>%
  slice_min(n = 25, order_by = relative.corr)
toppy <- subset(scale.drugs, rownames(scale.drugs) %in% unique(tippy$drug))
moa <- data.frame(
  row.names = rownames(toppy),
  moa = c(
    "antiandrogen", "ALK/RET.i", "PI3K.i", "aromatase.i", "antiandrogen",
    "MEK.i", "MEK.i", "antiandrogen", "BRAF.i", "nucleoside",
    "antiandrogen", "antiandrogen", "topII.i", "antiandrogen", "aromatase.i",
    "EGFR.i", "nucleoside", "VEGFR.i", "aromatase.i", "HER2/EGFR.i",
    "PARP.i", "CDK4/6.i", "kinase.i", "PARP.i", "SERM",
    "alkylator", "mTOR.i", "SERM", "MEK.i", "SMO.i"
  )
)
moa.cols <- c(colorRampPalette(paletteer_d("MoMAColors::Klein"))(18))
names(moa.cols) <- sort(unique(moa$moa))
library(RColorBrewer)
library(paletteer)
pdf("Output/Figures/12_isosceles/10_VS.drugs.launched_scaled_2024.11.11.pdf", height = 6, width = 6)
pheatmap(
  toppy*-1,
  annotation_col = ann,
  annotation_row = moa,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B"),
                           moa = moa.cols),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  show_colnames = FALSE,
  angle_col = "45", fontsize_row = 8, fontsize_col = 15
)
dev.off()

pdf("Output/Figures/12_isosceles/11_VS.drugs.launched_scaled_2024.11.11.pdf", height = 6, width = 10)
pheatmap(
  t(toppy*-1),
  # annotation_col = moa,
  annotation_row = ann,
  annotation_colors = list(VS = c("VS1" = "#1B9E77", "VS2" = "#9B58A5", "VS3" = "#BBA90B")),
  color = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(100)),
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = "45", fontsize_row = 8, fontsize_col = 8
)
dev.off()
# THIS REPLICATED 08.2 (GOOD VALIDATION)


