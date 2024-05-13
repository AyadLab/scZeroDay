library(tidyverse)
library(viridis)
library(cowplot)
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(reactome.db)
library(ComplexHeatmap)
library(RColorBrewer)

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
# READ IN DATA

def.genes.100 <- read.csv(
  "Output/Rdata/05_analysis_2024.05.03/03_VC.top.100_2023.09.26.csv",
  row.names = 1, header = TRUE
)

################################################################################
# ENRICHMENT ANALYSIS -- ALL SIG DEG

f1 <- def.genes.100[which(def.genes.100$cluster == "mVC1"), ] %>%
  arrange(desc(avg_log2FC))
f2 <- def.genes.100[which(def.genes.100$cluster == "mVC2"), ] %>%
  arrange(desc(avg_log2FC))
f3 <- def.genes.100[which(def.genes.100$cluster == "mVC3"), ] %>%
  arrange(desc(avg_log2FC))

library(org.Hs.eg.db)
f1.entrez <- mapIds(keys = f1$gene, org.Hs.eg.db, keytype = "SYMBOL", column = "ENTREZID")
f2.entrez <- mapIds(keys = f2$gene, org.Hs.eg.db, keytype = "SYMBOL", column = "ENTREZID")
f3.entrez <- mapIds(keys = f3$gene, org.Hs.eg.db, keytype = "SYMBOL", column = "ENTREZID")

clust.list <- list(
  mVC1 = f1$gene,
  mVC2 = f2$gene,
  mVC3 = f3$gene
)

clust.list.entrez <- list(
  mVC1 = f1.entrez,
  mVC2 = f2.entrez,
  mVC3 = f3.entrez
)

msigdb <- msigdbr(species = "Homo sapiens", category = "H")

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

comp.KEGG <- compareCluster(
  clust.list.entrez,
  fun = "enrichKEGG",
  organism = "hsa",
  qvalueCutoff = 0.05
)

comp.GO <- compareCluster(
  clust.list.entrez,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  qvalueCutoff = 0.05
)

dir.create("Output/Figures/05_analysis_2024.05.03/enrichment")

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/01_Essential.Clusters_MSigDB_2023.09.22.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  showCategory = 5,
  label_format = 50
) +
  theme_bw() +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, colour = "black")
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/02_Essential.Clusters_KEGG_2023.09.22.pdf", height = 8, width = 10)
dotplot(
  comp.KEGG,
  x = "Cluster",
  color = "p.adjust",
  showCategory = 5,
  label_format = 50
) +
  theme_bw() +
  mrd.theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, colour = "black")
  )
dev.off()

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/03_Essential.Clusters_GO_2023.09.22.pdf", height = 8, width = 10)
dotplot(
  comp.GO,
  x = "Cluster",
  color = "p.adjust",
  showCategory = 5,
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
# ENRICHMENT ANALYSIS, IMPROVED GENE SETS

dir.create("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich")

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# create function to isolate database name for filtering
funni <- function(X) {
  db <- str_split(X, pattern = "_")[[1]][1]
  return(db)
}

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/01_msigdb.C2_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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

# ONTOLOGY GENE SETS

msigdb <- msigdbr(species = "Homo sapiens", category = "C5")

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/05_msigdb.C5_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/02_msigdb.C3.TF_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/03_msigdb.C4_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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

# ONCOGENIC GENE SETS

msigdb <- msigdbr(species = "Homo sapiens", category = "C6")

msigdb_ref <- msigdb %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/04_msigdb.C6.ONC_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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


################################################################################
# IMPROVED ENRICHMENT ANALYSIS X2 FINAL FINAL FINAL

glimpse(clust.list)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb2 <- msigdbr(species = "Homo sapiens", category = "H")

msigdb3 <- rbind(msigdb, msigdb2)

msigdb_ref <- msigdb3 %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/06_msigdb.C2.H_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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

glimpse(comp@compareClusterResult)
table(comp@compareClusterResult$Cluster)

# filter only for databases of interest
filt <- comp[which(comp@compareClusterResult$DB == "REACTOME" | comp@compareClusterResult$DB == "KEGG" | comp@compareClusterResult$DB == "HALLMARK"), ]
comp.filt <- comp
comp.filt@compareClusterResult <- filt
# View(comp.filt@compareClusterResult)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/07_msigdb.C2.H_filt_ORA_2024.03.13.pdf", height = 8, width = 10)
dotplot(
  comp.filt,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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

glimpse(clust.list)

################################################################################
# ENRICHMENT ANALYSIS H AND C5

msigdb <- msigdbr(species = "Homo sapiens", category = "C5")
msigdb2 <- msigdbr(species = "Homo sapiens", category = "H")

msigdb3 <- rbind(msigdb, msigdb2)

msigdb_ref <- msigdb3 %>%
  distinct(gs_name, gene_symbol) %>%
  as.data.frame()

comp <- compareCluster(
  clust.list,
  fun = "enricher",
  TERM2GENE = msigdb_ref,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

comp@compareClusterResult$DB <- as.character(lapply(comp@compareClusterResult$ID, FUN = funni))
table(comp@compareClusterResult$DB)

pdf("Output/Figures/05_analysis_2024.05.03/enrichment/pan_enrich/08_msigdb.C5.H_ORA_2024.05.06.pdf", height = 8, width = 10)
dotplot(
  comp,
  x = "Cluster",
  color = "p.adjust",
  # split = "DB",
  showCategory = 5,
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

################################################################################
# ENRICHMENT ANALYSIS EXPLORATION





