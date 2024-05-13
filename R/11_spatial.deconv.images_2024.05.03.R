library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(paletteer)
library(SPOTlight)

################################################################################
# LOAD FILES

obj <- readRDS("Output/Rdata/04_merged.joined.scored_2024.04.04.RDS")
glimpse(obj@meta.data)

markers <- readRDS("data/03_VS.markers_2024.04.05.RDS")

files <- list.files(
  path = "Output/Rdata/05_deconv", 
  pattern = "^prop_data_v1_"
)


dir = "Output/Rdata/05_deconv_zero/"
for (i in 1:length(files)) {
  
  prop <- readRDS(paste0(dir, files[i]))
  add <- prop$prop.est
  colnames(add) <- c("propVS1", "propAstro", "propEndo", "propFibro", "propMyelo", "propOligo", "propT", "propVS2", "propVS3")
  
  obj <- AddMetaData(obj, metadata = add)
  
}

################################################################################

# make neoplastic vs. non-neoplastic prop metadata
obj@meta.data <- obj@meta.data %>%
  mutate(
    propNP = propVS1 + propVS2 + propVS3
  ) 
obj@meta.data <- obj@meta.data %>%
  mutate(
    propNON = propAstro + propEndo + propFibro + propMyelo + propOligo + propT
  ) 
glimpse(obj@meta.data)

saveRDS(obj, "Output/Rdata/08_deconv_v5_2024.05.09.RDS")

# # remove spots w/ < 25% predicted neoplastic content
# # also removes entry cortex slice spots
# neo <- WhichCells(obj, expression = propNP >= 0.25) # double check that this works
# obj.sub <- subset(obj, cells = neo)
# glimpse(obj.sub@meta.data)

# make data frame of just proportions data, NP v Non-NP
props <- data.frame(row.names = rownames(obj@meta.data), Non = obj$propNON, Neo = obj$propNP)
# make data frames of props w/ all data, mixed depending on focus
props.all <- data.frame(
  row.names = rownames(obj@meta.data), 
  VS1 = obj$propVS1, 
  VS2 = obj$propVS2, 
  VS3 = obj$propVS3, 
  Astrocytes = obj$propAstro, 
  Endothelial = obj$propEndo, 
  Fibroblast = obj$propFibro, 
  Myeloid = obj$propMyelo, 
  Oligodendrocyte = obj$propOligo, 
  T.Cell = obj$propT
)
props.state <- data.frame(
  row.names = rownames(obj@meta.data), 
  VS1 = obj$propVS1, 
  VS2 = obj$propVS2, 
  VS3 = obj$propVS3, 
  Non = obj$propNON
)
props.tme <- data.frame(
  row.names = rownames(obj@meta.data), 
  Neo = obj$propNP, 
  Astrocytes = obj$propAstro, 
  Endothelial = obj$propEndo, 
  Fibroblast = obj$propFibro, 
  Myeloid = obj$propMyelo, 
  Oligodendrocyte = obj$propOligo, 
  T.Cell = obj$propT
)

################################################################################
# plots

img <- Images(obj)
run <- img[grep("_C_", img, invert = TRUE)] # remove entry cortex slices

Idents(obj) <- obj$orig.ident

dir1 <- "Output/Figures/05_deconv_zero/01.np_non"
dir.create(dir1)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(obj, idents = slice)
  pdf(paste0(dir1, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = obj, 
      y = props[spots, ],  
      cell_types = colnames(props), 
      img = TRUE, slice = run[i],
      pie_scale = 0.35, scatterpie_alpha = 0.8
    ) + 
      scale_fill_manual(
        values = paletteer_d(palette = "ggsci::lanonc_lancet", length(colnames(props))), 
        breaks = colnames(props)
      )
  )
  dev.off()
}

dir2 <- "Output/Figures/05_deconv_zero/02.all"
dir.create(dir2)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(obj, idents = slice)
  pdf(paste0(dir2, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = obj, 
      y = props.all[spots, ],  
      cell_types = colnames(props.all), 
      img = TRUE, slice = run[i],
      pie_scale = 0.35, scatterpie_alpha = 0.8
    ) + 
      scale_fill_manual(
        values = paletteer_d(palette = "ggsci::lanonc_lancet", length(colnames(props.all))), 
        breaks = colnames(props.all)
      )
  )
  dev.off()
}

dir3 <- "Output/Figures/05_deconv_zero/03.state"
dir.create(dir3)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(obj, idents = slice)
  pdf(paste0(dir3, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = obj, 
      y = props.state[spots, ],  
      cell_types = colnames(props.state), 
      img = TRUE, slice = run[i],
      pie_scale = 0.35, scatterpie_alpha = 0.9
    ) + 
      scale_fill_manual(
        values = paletteer_d(palette = "RColorBrewer::Dark2", length(colnames(props.state))), 
        breaks = colnames(props.state)
      )
  )
  dev.off()
}

dir3.1 <- "Output/Figures/05_deconv_zero/03.state.emty"
dir.create(dir3.1)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(obj, idents = slice)
  pdf(paste0(dir3.1, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = obj, 
      y = props.state[spots, ],  
      cell_types = colnames(props.state), 
      img = FALSE, slice = run[i], ####
      pie_scale = 0.40, scatterpie_alpha = 1.0
    ) + 
      scale_fill_manual(
        values = paletteer_d(palette = "RColorBrewer::Dark2", length(colnames(props.state))), 
        breaks = colnames(props.state)
      )
  )
  dev.off()
}

dir4 <- "Output/Figures/05_deconv_zero/04.tme"
dir.create(dir4)
for (i in 1:length(run)) {
  slice <- run[i]
  spots <- WhichCells(obj, idents = slice)
  pdf(paste0(dir4, "/prop_image_v5_", slice, ".pdf"))
  print(
    plotSpatialScatterpie(
      x = obj, 
      y = props.tme[spots, ],  
      cell_types = colnames(props.tme), 
      img = TRUE, slice = run[i],
      pie_scale = 0.35, scatterpie_alpha = 0.9
    ) + 
      scale_fill_manual(
        values = paletteer_d(palette = "ggsci::lanonc_lancet", length(colnames(props.tme))), 
        breaks = colnames(props.tme)
      )
  )
  dev.off()
}

################################################################################
# quantification

my.theme <- theme(
  title = element_blank(),
  plot.title = element_blank(),
  plot.subtitle = element_blank(), plot.caption = element_blank(),
  # legend.position = "none",
  axis.text = element_text(size = 11, colour = "black"),
  axis.text.x = element_text(size = 11, colour = "black"),
  axis.text.y = element_text(size = 11, colour = "black"),
  axis.title = element_text(size = 12, colour = "black"),
  axis.line = element_line(linewidth = 1, colour = "black"),
  axis.ticks.length = unit(0.1, "cm"),
  axis.ticks = element_line(linewidth = 1, colour = "black")
)

quant <- data.frame(
  row.names = rownames(obj@meta.data), 
  sliceID = obj@meta.data$orig.ident, 
  VS1 = obj@meta.data$propVS1, 
  VS2 = obj@meta.data$propVS2, 
  VS3 = obj@meta.data$propVS3, 
  Astrocytes = obj@meta.data$propAstro, 
  Endothelial = obj@meta.data$propEndo, 
  Fibroblast = obj@meta.data$propFibro, 
  Myeloid = obj@meta.data$propMyelo, 
  Oligodendrocyte = obj@meta.data$propOligo, 
  T.Cell = obj@meta.data$propT
)

# remove "C" slices ("NA" values in proportions)
quant <- quant[complete.cases(quant), ]
glimpse(quant)

quant.long <- quant %>%
  pivot_longer(cols = -sliceID, names_to = "CellType", values_to = "Proportion")
glimpse(quant.long)

quant.ag <- quant.long %>%
  aggregate(Proportion ~ sliceID + CellType, FUN = mean)
glimpse(quant.ag)
quant.ag$CellType <- factor(quant.ag$CellType, levels = c("VS1", "VS2", "VS3", "Astrocytes", "Endothelial", "Fibroblast", "Myeloid", "Oligodendrocyte", "T.Cell"))

pdf("Output/Figures/05_deconv_zero/02.all/01_prop.quant_2024.04.18.pdf", height = 7, width = 14)
quant.ag %>%
  ggplot(aes(x = sliceID, y = Proportion, fill = CellType)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(
    values = paletteer_d(
      palette = "RColorBrewer::Set1", 
      n = ifelse(length(levels(quant.ag$CellType)) <= 9, length(levels(quant.ag$CellType)), 9)
    ), 
    breaks = levels(quant.ag$CellType)
  ) + 
  theme_minimal() + 
  my.theme + 
  RotatedAxis()
dev.off()

non <- c("Astrocytes", "Endothelial", "Fibroblast", "Myeloid", "Oligodendrocyte", "T.Cell")
quant.ag <- quant.ag[-which(quant.ag$CellType %in% non), ]

quant.np <- quant.ag %>%
  pivot_wider(names_from = CellType, values_from = Proportion)
quant.np <- quant.np %>%
  mutate(
    sum = VS1+VS2+VS3
  )
quant.np <- quant.np %>%
  mutate(
    fct = 1/sum
  )
quant.np <- quant.np %>%
  mutate(
    newVS1 = VS1*fct, 
    newVS2 = VS2*fct, 
    newVS3 = VS3*fct
  )
glimpse(quant.np)
quant.np <- quant.np[, -c(2:6)]
colnames(quant.np) <- c("sliceID", "VS1", "VS2", "VS3")

quant.long <- quant.np %>%
  pivot_longer(cols = -sliceID, names_to = "CellType", values_to = "Proportion")
glimpse(quant.long)

quant.ag <- quant.long %>%
  aggregate(Proportion ~ sliceID + CellType, FUN = mean)
glimpse(quant.ag)
quant.ag$CellType <- factor(quant.ag$CellType, levels = c("VS1", "VS2", "VS3"))

pdf("Output/Figures/05_deconv_zero/02.all/02_prop.quant_np.only_2024.04.18.pdf", height = 7, width = 14)
quant.ag %>%
  ggplot(aes(x = sliceID, y = Proportion, fill = CellType)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(
    values = paletteer_d(palette = "RColorBrewer::Dark2", length(levels(quant.ag$CellType))), 
    breaks = levels(quant.ag$CellType)
  ) + 
  theme_minimal() + 
  my.theme + 
  RotatedAxis()
dev.off()


