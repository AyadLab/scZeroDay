library(Seurat)
library(tidyverse)
library(cowplot)
library(viridis)
library(singscore)
library(scran)
library(ggplot2)
library(SPOTlight)
library(RColorBrewer)
library(SPOTlight)
library(Biobase)
suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

################################################################################

i <- Sys.getenv("i")
i <- as.numeric(i)
print(paste("This is iteration", i, "of 18")) # 18 tumor sections (non-ctx)

################################################################################

setwd("/data/mrd/spatial")

obj <- readRDS("Output/Rdata/04_merged.joined.scored_2024.04.04.RDS")

# non-regressed scale.data
non.regr <- readRDS("Output/Rdata/06_scale.data.not.regressed_2024.04.05.RDS")
obj@assays$Spatial$scale.data <- non.regr

################################################################################
img <- Images(obj)
run <- img[grep("_C_", img, invert = TRUE)] # remove entry cortex slices

inx <- which(Images(obj) == run[i])

image.name <- run[i]
image.metadata <- obj@images[[inx]]
spot.ids <- rownames(image.metadata@coordinates)

subset.obj <- subset(obj, cells = spot.ids)
keys.to.remove <- setdiff(Images(subset.obj), run[i])
for (key in keys.to.remove) {
  subset.obj@images[[key]] <- NULL
}

################################################################################

# ref <- readRDS("/data/mrd/spatial/Suter_scGBM_2023.10.04.RDS") # updated new suter obj 5/7/2024
# 
# ref@meta.data <- ref@meta.data %>%
#   mutate(
#     deconv = case_when(
#       NPvNon == "Neoplastic" ~ mVC, # "Neoplastic" <- replace "Neoplastic" w/ this to do state-based deconv
#       NPvNon == "Non-Neoplastic" ~ CellType
#     )
#   )
# table(ref$deconv)
# 
# eset.ref <- ExpressionSet(
#   assayData = as.matrix(ref@assays$RNA$counts),
#   phenoData = AnnotatedDataFrame(ref@meta.data)
# )
# saveRDS(eset.ref, "Output/Rdata/suter.ref.eset_2024.05.07.RDS")

eset.ref <- readRDS("Output/Rdata/suter.ref.eset_2024.05.07.RDS")

eset.sp <- ExpressionSet(
  assayData = as.matrix(subset.obj@assays$Spatial$counts), 
  phenoData = AnnotatedDataFrame(subset.obj@meta.data)
)

genes.ref <- rownames(exprs(eset.ref))
genes.sp <- rownames(exprs(eset.sp))
common.genes <- intersect(genes.ref, genes.sp)

eset.ref.subset <- eset.ref[common.genes, ]
eset.sp.subset <- eset.sp[common.genes, ]

################################################################################

deconv <- SCDC::SCDC_prop_subcl_marker(
  eset.sp.subset,
  eset.ref.subset,
  ct.varname = "deconv",
  sample = "orig.ident",
  fl.varname = "mVC", # <- "NPvNon"
  ct.sub = as.character(unique(eset.ref.subset$deconv)),
  ct.fl.sub = as.character(unique(eset.ref.subset$mVC)), # <- $NPvNon
  weight.basis = TRUE,
  iter.max = 100
)

# # IF ABOVE DOES NOT WORK -- DECONV NP V NON ONLY
# deconv <- SCDC::SCDC_prop(
#   bulk.eset = eset.sp.subset,
#   sc.eset = eset.ref.subset,
#   ct.varname = "NPvNon",
#   ct.sub = as.character(unique(eset.ref.subset$NPvNon)),
#   weight.basis = TRUE,
#   ct.cell.size = NULL,
#   select.marker = TRUE,
#   iter.max = 100
# )

################################################################################

prop.level <- deconv$prop.est

prop.image <- plotSpatialScatterpie(
  x = subset.obj, 
  y = prop.level, 
  cell.types = colnames(prop.level), 
  img = TRUE, 
  pie.scale = 0.4
) + 
  scale_fill_manual(
    values = viridis_pal(option = "turbo")(length(colnames(prop.level))), 
    breaks = colnames(prop.level)
  )

output.dir <- "/data/mrd/spatial/Output/Figures/05_deconv_zero/"
output.file_deconv <- paste0(output.dir, "prop_image_v1_", gsub("[^A-Za-z0-9]", "_", image.name), ".pdf")
pdf(output.file_deconv)
prop.image
dev.off()



data.dir <- "/data/mrd/spatial/Output/Rdata/05_deconv_zero/"
saveRDS(
  deconv, 
  paste0(data.dir, "prop_data_v1_", gsub("[^A-Za-z0-9]", "_", image.name), ".RDS")
)

print(paste("Iteration", i, "of 18 has COMPLETED"))
