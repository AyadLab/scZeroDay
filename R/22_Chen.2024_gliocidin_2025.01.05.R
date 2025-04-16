library(tidyverse)
library(Seurat)
library(future)
library(singscore)
library(GSEABase)
library(presto)
library(cowplot)
library(ggpubr)


options(future.globals.maxSize = +Inf)

################################################################################
obj <- readRDS(
  "data/gliocidin/GSE277301_integr.Veh.Tmz.Gliocidin.Combo_Seurat.obj.rds"
)

ess <- read.csv(
  "/data/mrd/cognitive.seeds_zero/Output/Rdata/01_essential.genes/05_essential.genes_killing_2024.01.31.csv",
  row.names = 1
)

# vs markers
markers <- readRDS("/data/mrd/cognitive.seeds_zero/Output/Rdata/04_states.2_2024.05.02/03_VS.markers_2024.04.05.RDS")
names(markers)

# gds genes in venn
gdsXvenn <- read.csv(
  "/data/mrd/cognitive.seeds_zero/Output/Rdata/19_gds.dive/gdsXvenn_genes_2024.11.21.csv",
  row.names = 1
)

# top 100 venn genes
topVenn <- read.csv(
  "/data/mrd/cognitive.seeds_zero/Output/Rdata/10_deg.corr/04_VS.top100.overlapping_2024.11.12.csv",
  row.names = 1
)
glimpse(topVenn)

# all venn genes
venn <- read.csv(
  "/data/mrd/cognitive.seeds_zero/Output/Rdata/10_deg.corr/04_VS.overlapping.genes_2024.11.12.csv",
  row.names = 1
)

################################################################################
# CAN'T RE-PROCESS DATA...DON'T HAVE TUMOR-LEVEL INFORMATION

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

plan("multisession", workers = 2)
obj <- ScaleData(
  obj,
  features = rownames(obj[["RNA"]]$counts),
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
)
obj.data <- obj[["RNA"]]$scale.data
saveRDS(obj.data, "/data/mrd/seed.wrangler/Output/Rdata/gliociding_scale.data_2025.01.20.RDS")

################################################################################
# set up ranked data based on scRNA-seq signatures
# obj.data <- GetAssayData(um2, assay = "RNA", layer = "scale.data")
# saveRDS(obj.data, "Output/Rdata/20_um002/00_um002vDMSO_scale.data_cycle.regr_2024.11.24.RDS")
obj.data.df <- as.data.frame(obj.data)
rank.data <- singscore::rankGenes(obj.data.df)
# rm(obj.data)
rm(obj.data.df)
gc()

# # score each cell for each cell state
# for (i in 1:length(markers)) { # change to do top VENN genes
#   # capture unique gene signature
#   sig <- markers[[i]]$gene
#   # capture name of cell state ,
#   state <- names(markers)[i]
#   # singscore scoring
#   set <- GeneSet()
#   set@geneIds <- as.character(sig)
#   scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
#   # rename $Sig to appropriate cell state
#   scored$Sig <- as.character(state)
#   assign(paste0("scored.", state), scored)
# }
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

# can't really do alluvial -- should examine gene expression changes
table(obj$VS, obj$orig.ident)


# dir.create("UM002/Output")
saveRDS(
  obj@meta.data,
  "/data/mrd/seed.wrangler/Output/Rdata/gliocidin_meta.data_2025.01.21.RDS"
)

################################################################################
# differential expression between treatment groups
Idents(obj) <- obj$orig.ident
deg <- FindAllMarkers(
  obj,
  assay = "RNA",
  logfc.threshold = 0,
  test.use = "wilcox", only.pos = TRUE
)
sig <- deg[which(deg$p_val_adj < 0.05), ]

write.csv(
  deg,
  "/data/mrd/seed.wrangler/Output/Rdata/01_tx.deg_2025.01.21.csv",
  row.names = TRUE
)
write.csv(
  sig,
  "/data/mrd/seed.wrangler/Output/Rdata/02_tx.deg.sig_2025.01.21.csv",
  row.names = TRUE
)










################################################################################

meta <- readRDS(
  "/data/mrd/seed.wrangler/Output/Rdata/gliocidin_meta.data_2025.01.21.RDS"
)

plot <- data.frame(
  row.names = rownames(meta),
  Tx = meta$orig.ident,
  VS = meta$VS,
  VS1.score = meta$mVC1_singscore,
  VS2.score = meta$mVC2_singscore,
  VS3.score = meta$mVC3_singscore
)
glimpse(plot)

# vehicle mean and sd calculations
vs1.dmso.u <- mean(plot$VS1.score[which(plot$Tx == "Veh")])
vs1.dmso.sd <- sd(plot$VS1.score[which(plot$Tx == "Veh")])
vs2.dmso.u <- mean(plot$VS2.score[which(plot$Tx == "Veh")])
vs2.dmso.sd <- sd(plot$VS2.score[which(plot$Tx == "Veh")])
vs3.dmso.u <- mean(plot$VS3.score[which(plot$Tx == "Veh")])
vs3.dmso.sd <- sd(plot$VS3.score[which(plot$Tx == "Veh")])

new <- data.frame()
props <- data.frame()
tx <- c("Tmz", "Cpd", "Combo")
for (i in 1:length(tx)) {
  trt <- tx[i]
  tx.df <- plot[which(plot$Tx == trt), ]

  # adjusted VS enrichment scores by treatment (should be by tumor as well)
  tx.df$VS1 <- (tx.df$VS1.score-vs1.dmso.u)/vs1.dmso.sd
  tx.df$VS2 <- (tx.df$VS2.score-vs2.dmso.u)/vs2.dmso.sd
  tx.df$VS3 <- (tx.df$VS3.score-vs3.dmso.u)/vs3.dmso.sd
  new <- rbind(new, tx.df)

  # VS proportions by treatment (should be by tumor as well)
  prop.VS1 <- length(rownames(tx.df)[which(tx.df$VS == "mVC1")])/length(rownames(tx.df))
  prop.VS2 <- length(rownames(tx.df)[which(tx.df$VS == "mVC2")])/length(rownames(tx.df))
  prop.VS3 <- length(rownames(tx.df)[which(tx.df$VS == "mVC3")])/length(rownames(tx.df))
  add <- data.frame(
    Tx = trt,
    VS1 = prop.VS1,
    VS2 = prop.VS2,
    VS3 = prop.VS3
  )
  props <- rbind(props, add)
}
glimpse(new)
glimpse(props)

# adjust for tumor proportions
# (shouldn't be doing this for different treatments -- I don't have tumor data)
# may acutally have to do this to make comparable across tx groups?
new.adj <- new %>%
  left_join(props, by = "Tx", suffix = c("", "_prop")) %>%
  mutate(
    VS1.adj = VS1*(1/VS1_prop),
    VS2.adj = VS2*(1/VS2_prop),
    VS3.adj = VS3*(1/VS3_prop)
  )
glimpse(new.adj)

# mns <- new %>%
#   group_by(Tx) %>%
#   summarise(VS1 = mean(VS1), VS2 = mean(VS2), VS3 = mean(VS3))
# mns
#
# # create data frames for plotting
# dat.summ <- new %>%
#   group_by(Tx) %>%
#   summarise(
#     VS1.mean = mean(VS1), VS1.se = sd(VS1)/sqrt(n()),
#     VS2.mean = mean(VS2), VS2.se = sd(VS2)/sqrt(n()),
#     VS3.mean = mean(VS3), VS3.se = sd(VS3)/sqrt(n()),
#   )
# dat.summ
mns <- new.adj %>%
  group_by(Tx) %>%
  summarise(VS1 = mean(VS1.adj), VS2 = mean(VS2.adj), VS3 = mean(VS3.adj))
mns

# create data frames for plotting
dat.summ <- new.adj %>%
  group_by(Tx) %>%
  summarise(
    VS1.mean = mean(VS1.adj), VS1.se = sd(VS1.adj)/sqrt(n()),
    VS2.mean = mean(VS2.adj), VS2.se = sd(VS2.adj)/sqrt(n()),
    VS3.mean = mean(VS3.adj), VS3.se = sd(VS3.adj)/sqrt(n()),
  )
dat.summ

long_data <- dat.summ %>%
  pivot_longer(
    cols = -Tx,  # Exclude 'Tx' from reshaping
    names_to = c("VS", ".value"),  # Split column names
    names_sep = "\\."  # Separator between VS and metric
  )

for (i in 1:length(unique(long_data$Tx))) {
  trt <- unique(long_data$Tx)[i]
  fig <- long_data[which(long_data$Tx == trt), ] %>%
    ggplot(aes(x = VS, y = mean, fill = VS)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +  # Barplot
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black") +  # Error bars
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +  # Clean theme
    labs(
      title = paste0("Treatment: ", trt),
      x = "Vulnerability State",
      y = "Normalized Enrichment Score"
    )
  assign(paste0("fig.", trt), fig)
}


# plot
pdf("/data/mrd/seed.wrangler/Output/Figures/01_Tx.effect.by.VS_venn.markers_2025.01.21.pdf", height = 4, width = 12)
cowplot::plot_grid(
  fig.Tmz,
  fig.Cpd,
  fig.Combo,
  align = "h", nrow = 1, ncol = 3
)
dev.off()


################################################################################
# WITH STATISTICS

# anova.res <- new.adj %>%
#   group_by(Tx) %>%
#   summarise(
#     p_VS1 = list(aov(VS1.adj ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy()),
#     p_VS2 = list(aov(VS2.adj ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy()),
#     p_VS3 = list(aov(VS3.adj ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy())
#   )
# anova.res
# View(anova.res)




# anova.tx <- list(
#   VS1 = summary(aov(VS1.adj ~ Tx, data = new.adj))[[1]]$`Pr(>F)`[1],
#   VS2 = summary(aov(VS2.adj ~ Tx, data = new.adj))[[1]]$`Pr(>F)`[1],
#   VS3 = summary(aov(VS3.adj ~ Tx, data = new.adj))[[1]]$`Pr(>F)`[1]
# )
# anova.tx


new.adj <- new.adj %>%
  rename(VS = "VS_Category")
# Convert `new.adj` to long format for raw data points
long_raw_data <- new.adj %>%
  pivot_longer(cols = c(VS1.adj, VS2.adj, VS3.adj),
               names_to = "VS", values_to = "score") %>%
  mutate(VS = gsub(".adj", "", VS))  # Rename for consistency

for (i in 1:length(unique(long_raw_data$Tx))) {
  trt <- unique(long_raw_data$Tx)[i]

  fig <- long_raw_data %>%
    filter(Tx == trt) %>%
    ggplot(aes(x = VS, y = score, fill = VS)) +
    # Add jittered individual points
    geom_jitter(aes(color = VS), width = 0.2, size = 1.5, alpha = 0.6) +
    geom_boxplot(aes(fill = VS)) +
    # # Add bar plot for means
    # stat_summary(fun = mean, geom = "bar", width = 0.7, color = "black", alpha = 0.5) +
    # # Add error bars for standard error
    # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +  # Ensure consistent colors
    theme_minimal() +
    labs(
      title = paste0("Treatment: ", trt),
      x = "Vulnerability State",
      y = "Normalized Enrichment Score"
    ) +
    # Statistical comparisons for within-treatment differences
    stat_compare_means(
      method = "anova",
      label.y = max(long_raw_data$score, na.rm = TRUE) + 0.5
    ) +
    stat_compare_means(
      comparisons = list(c("VS1", "VS2"), c("VS1", "VS3"), c("VS2", "VS3")),
      method = "t.test", label = "p.signif", na.rm = TRUE
    )

  assign(paste0("fig.", trt), fig)
}

# fig.Tmz

# Save the updated plots
pdf("/data/mrd/seed.wrangler/Output/Figures/02_Tx.effect.by.VS_venn.markers_stats_2025.02.06.pdf", height = 8, width = 16)
cowplot::plot_grid(
  fig.Tmz,
  fig.Cpd,
  fig.Combo,
  align = "h", nrow = 1, ncol = 3
)
dev.off()


anova.res <- long_raw_data %>%
  group_by(Tx) %>%
  filter(Tx == "Tmz") %>%
  summarise(
    AOV = list(aov(score ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy())
  )
anova.res[[2]][1]
# term  contrast null.value estimate conf.low conf.high   adj.p.value
# <chr> <chr>         <dbl>    <dbl>    <dbl>     <dbl>         <dbl>
#   1 VS    VS2-VS1           0   -2.52    -2.70     -2.35  0
# 2 VS    VS3-VS1           0    0.473    0.295     0.651 0.00000000129
# 3 VS    VS3-VS2           0    3.00     2.82      3.18  0

anova.res <- long_raw_data %>%
  group_by(Tx) %>%
  filter(Tx == "Cpd") %>%
  summarise(
    AOV = list(aov(score ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy())
  )
anova.res[[2]][1]
# term  contrast null.value estimate conf.low conf.high adj.p.value
# <chr> <chr>         <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#   1 VS    VS2-VS1           0   -0.460   -0.755    -0.166    7.21e- 4
# 2 VS    VS3-VS1           0   -1.14    -1.43     -0.842    4.03e-14
# 3 VS    VS3-VS2           0   -0.676   -0.970    -0.382    2.17e- 7

anova.res <- long_raw_data %>%
  group_by(Tx) %>%
  filter(Tx == "Combo") %>%
  summarise(
    AOV = list(aov(score ~ VS, data = .) %>% TukeyHSD() %>% broom::tidy())
  )
anova.res[[2]][1]
# term  contrast null.value estimate conf.low conf.high adj.p.value
# <chr> <chr>         <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
#   1 VS    VS2-VS1           0    -2.19    -2.29     -2.08           0
# 2 VS    VS3-VS1           0    -5.36    -5.47     -5.26           0
# 3 VS    VS3-VS2           0    -3.18    -3.28     -3.07           0
