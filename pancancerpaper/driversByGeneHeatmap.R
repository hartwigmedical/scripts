library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

load(file = "~/hmf/RData/Processed/germlineCatalog.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")
load("~/hmf/RData/processed/driverGenes.RData")
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
load(file = "~/hmf/RData/reference/simplifiedDrivers.RData")

tsgDrivers = c("Deletion","FragileDel","Indel","Missense","Multihit","Nonsense","Splice")
oncoDrivers = c("Amplification","Deletion","Fusion","Indel","Missense","Promoter")
germlineDrivers = c("Deletion", "Indel","Missense", "Multihit","Nonsense", "Splice","Synonymous")
simplifiedDrivers <- simplifiedDrivers[!simplifiedDrivers == "Germline"]

tsgDriverColours = simplifiedDriverColours[tsgDrivers]
oncoDriverColours = simplifiedDriverColours[oncoDrivers]
germlineDriverColours = simplifiedDriverColours[germlineDrivers]

germlineDriversByGene = germlineDriverCatalog %>% filter(highConfidenceGenePanel, is.na(variantLostInTumor) | !variantLostInTumor)

hpcDriversByGene = hpcDriversByGene  %>% 
  filter(driverLikelihood > 0) %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
    driver = ifelse(driver == "Inframe", "Indel", driver),
    hotspot = ifelse(driver == "Promoter", T, hotspot)
  ) %>%
  ungroup()

sortedTsgGenes = hpcDriversByGene %>% filter(type == 'TSG') %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(20, n) %>% arrange(-n)
sortedOncoGenes = hpcDriversByGene %>% filter(type == 'ONCO') %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(20, n) %>% arrange(-n)
sortedCancerTypes = hpcDriversByGene %>% group_by(cancerType) %>% summarise(n = sum(driverLikelihood)) %>% arrange(-n)
sortedGermlineGenes = germlineDriversByGene %>% group_by(gene) %>% count() %>% ungroup() %>% top_n(10, n) %>% arrange(-n)

oncoDriversByGene = hpcDriversByGene %>% filter(type == 'ONCO', gene %in% sortedOncoGenes$gene) %>% mutate(driver = factor(driver, simplifiedDrivers))
tsgDriversByGene = hpcDriversByGene %>% filter(type == 'TSG', gene %in% sortedTsgGenes$gene)  %>% mutate(driver = factor(driver, simplifiedDrivers))
germlineDriversByGene = germlineDriversByGene %>% 
  filter(gene %in% sortedGermlineGenes$gene) %>%
  mutate(
    driver = ifelse(mutationClass %in% c("Frameshift", "Inframe"), "Indel", mutationClass),
    driver = factor(driver, simplifiedDrivers)) %>% mutate(driverLikelihood = 1)

main_heatmap_data <- function(sourceData) {
  heatmapData = sourceData %>%
    group_by(cancerType, gene) %>% summarise(n = sum(driverLikelihood)) %>%
    full_join(hpcCancerTypeCounts %>% select(cancerType, cancerTypeSamples = N), by = "cancerType") %>%
    mutate(proportion = n / cancerTypeSamples) %>%
    select(-n, -cancerTypeSamples) %>%
    spread(cancerType, proportion, fill = 0) %>% 
    filter(!is.na(gene))
  
  rownames(heatmapData) <- heatmapData$gene
  heatmapData$gene <- NULL
  return (heatmapData)
}

samples_annotation <- function(sourceData) {
  annotationData1 = sourceData %>% 
    group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% 
    mutate(proportion = round(n / sum(hpcCancerTypeCounts$N), 3)) %>% select(-n) 
  rownames(annotationData1) <- annotationData1$gene
  annotationData1$gene <- NULL
  annotationData1$non <- 0
  return (annotationData1)
}

driver_annotation <- function(sourceData) {
  annotationData2 = sourceData %>% 
    group_by(gene, driver) %>% summarise(n = sum(driverLikelihood)) %>% 
    group_by(gene) %>% mutate(total = sum(n)) %>% ungroup() %>%
    mutate(proportion = n / total) %>%
    select(gene, driver, proportion) %>%
    spread(driver, proportion, fill = 0)
  rownames(annotationData2) <- annotationData2$gene
  annotationData2$gene <- NULL
  return (annotationData2)
}


biallelic_annotation <- function(sourceData) {
  biallelicData = sourceData %>% 
    mutate(biallelic = ifelse(driver %in% c("Multihit","Amp","Del"), T, biallelic)) %>%
    mutate(biallelic = ifelse(is.na(biallelic), F, biallelic)) %>%
    group_by(gene, biallelic) %>% 
    summarise(n = sum(driverLikelihood)) %>% 
    mutate(biallelic = ifelse(biallelic, "biallelic", "notBiallelic")) %>%
    spread(biallelic, n, fill = 0) %>%
    mutate(biallelicPercentage =  biallelic / (biallelic + notBiallelic) ) %>%  
    select(gene, biallelicPercentage) %>%
    mutate(nonBiallelicPercentage = 1 - biallelicPercentage)
  
  rownames(biallelicData) <- biallelicData$gene
  biallelicData$gene <- NULL
  return (biallelicData)
}

wildtypelost_annotation <- function(sourceData) {
  wildData = sourceData %>% 
    mutate(wildTypeLostInTumor = ifelse(is.na(wildTypeLostInTumor), F, wildTypeLostInTumor)) %>%
    group_by(gene, wildTypeLostInTumor) %>% 
    summarise(n = sum(driverLikelihood)) %>% 
    mutate(wildTypeLostInTumor = ifelse(wildTypeLostInTumor, "lost", "notLost")) %>%
    spread(wildTypeLostInTumor, n, fill = 0) %>%
    mutate(lostPercentage =  lost / (lost + notLost) ) %>%  
    select(gene, lostPercentage) %>%
    mutate(notLostPercentage = 1 - lostPercentage)
  
  rownames(wildData) <- wildData$gene
  wildData$gene <- NULL
  return (wildData)
}


tsgHeatmapData = main_heatmap_data(tsgDriversByGene)
tsgSamplesAnnotationData = samples_annotation(tsgDriversByGene)
tsgDriversAnnotationData = driver_annotation(tsgDriversByGene)
tsgBiallelicAnnotationData = biallelic_annotation(tsgDriversByGene)

oncoHeatmapData = main_heatmap_data(oncoDriversByGene)
oncoSamplesAnnotationData = samples_annotation(oncoDriversByGene)
oncoDriversAnnotationData = driver_annotation(oncoDriversByGene)

germlineHeatmapData = main_heatmap_data(germlineDriversByGene)
germlineSamplesAnnotationData = samples_annotation(germlineDriversByGene)
germlineDriversAnnotationData = driver_annotation(germlineDriversByGene)
germlineWildtypeLostData = wildtypelost_annotation(germlineDriversByGene)

heat_map_text <- function(value) {
  value = round(100*value, 1)
  if (value > 0) {
    return (sprintf(fmt='%.1f', value))
  }
  return ("")
}

germlineMat = (data.matrix(germlineHeatmapData)) 
germlineHeatmap = Heatmap(
  column_title = "  ",
  germlineMat, 
  row_names_gp = gpar(fontsize = 5, fontface = "italic"),
  column_names_gp = gpar(fontsize = 5),
  row_order = sortedGermlineGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.33, 0.66, 0.1), c("#f7fcf5","#bae4b3","#74c476","#238b45")), 
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (germlineMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(heat_map_text(germlineMat[i, j]), x, y, gp=gpar(fontsize=5, col = myColor))
  }
)

germlineSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(
    baseline = "min", 
    germlineSamplesAnnotationData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.01, 0.02, 0.03), labels = c("0%", "1%", "2%", "3%"),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,0.03), 
    gp = gpar(fontsize = 3, fill = "#bc80bd", col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)


germlineDriversAnnotation = rowAnnotation(
  `% Drivers` = row_anno_barplot(
    germlineDriversAnnotationData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1), labels = c("", "25%", "50%", "75%", "100%"),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,1), 
    gp = gpar(fill = germlineDriverColours, col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)


germlineWildTypeLostAnnotationColours = setNames(c("#ccebc5","white"), c("lostPercentage", "notLostPercentage"))
germlineWildTypeLostAnnotation = rowAnnotation(
  show_annotation_name = T,
  annotation_name_rot = 0,
  annotation_name_gp  = gpar(fontsize = 5),
  `% Biallelic in Tumor` = row_anno_barplot(
    germlineWildtypeLostData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1), labels = c("", "25%", "50%", "75%", "100%"),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,1), 
    gp = gpar(fill = germlineWildTypeLostAnnotationColours, col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm")
)

pGermline = germlineHeatmap + germlineSamplesAnnotation + germlineDriversAnnotation+ germlineWildTypeLostAnnotation
pGermline


oncoMat = (data.matrix(oncoHeatmapData))
oncoHeatmap = Heatmap(
  column_title = "  ",
  oncoMat, 
  row_names_gp = gpar(fontsize = 5, fontface = "italic"),
  column_names_gp = gpar(fontsize = 0),
  row_order = sortedOncoGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("#fff5eb","#fdbe85","#fd8d3c","#d94701")),
  #col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("white", "lightblue", "blue", "darkblue")),
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (oncoMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(heat_map_text(oncoMat[i, j]), x, y, gp=gpar(fontsize=5, col = myColor))
  }
)

oncoSamplesAnnotation = rowAnnotation(
  `% Samples` = 
    row_anno_barplot(
      oncoSamplesAnnotationData, 
      axis_param = list(
        at = c(0, 0.1, 0.2), labels = c("0%", "10%", "20%"),
        gp = gpar(fontsize = 5),
        labels_rot = "0",
        side = "top"),
      axis = T, 
      ylim = c(0,0.2), 
      gp = gpar(fill = "#bc80bd", col = 0), 
      border = F,
      bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = F,

  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)

oncoDriversAnnotation = rowAnnotation(
  `% Drivers` = 
    row_anno_barplot(
      oncoDriversAnnotationData, 
      axis_param = list(
        at = c(0, 0.25, 0.5, 0.75, 1), labels = c("", "25%", "50%", "75%", "100%"),
        gp = gpar(fontsize = 5),
        labels_rot = "0",
        side = "top"),
      axis = T, 
      ylim = c(0,1), 
      gp = gpar(fill = oncoDriverColours, col = 0), 
      border = F,
      bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = F,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)

oncoDriversAnnotationIndex = rowAnnotation(
  show_annotation_name = F,
  df = data.frame("Driver" = simplifiedDrivers),
  col = list(Driver = simplifiedDriverColours),
  simple_anno_size = unit(0, "cm"),
  width = unit(2, "cm"),
  annotation_legend_param = list(title = "", labels_gp = gpar(fontsize = 5), grid_height = unit(2, "mm"), grid_width = unit(2, "mm"))
)

pOnco = oncoHeatmap + oncoSamplesAnnotation + oncoDriversAnnotation + oncoDriversAnnotationIndex
pOnco

tsgMat = (data.matrix(tsgHeatmapData))
tsgHeatmap = Heatmap(
  column_title = "  ",
  tsgMat, 
  #width = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 5, fontface = "italic"),
  column_names_gp = gpar(fontsize = 0),
  row_order = sortedTsgGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("#f7fbff","#bdd7e7","#6baed6","#2171b5")),
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (tsgMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(heat_map_text(tsgMat[i, j]), x, y, gp=gpar(fontsize=5, col = myColor))
  }
)

tsgSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(
    tsgSamplesAnnotationData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.15, 0.3, 0.45, 0.6), labels = c("0%", "15%", "30%", "45%", "60%"),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,0.6), 
    gp = gpar(fill = "#bc80bd", col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = F,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)

tsgBiallelicAnnotationColours = setNames(c("#ccebc5","white"), c("biallelicPercentage", "nonBiallelicPercentage"))
tsgBiallelicAnnotation = rowAnnotation(
  show_annotation_name = F,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0,
  `% Biallelic` = row_anno_barplot(
    tsgBiallelicAnnotationData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1), labels = c("", "25%", "50%", "75%", "100%  "),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,1), 
    gp = gpar(fill = tsgBiallelicAnnotationColours, col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm")
)

tsgDriversAnnotation = rowAnnotation(
  `% Drivers` = row_anno_barplot(
    tsgDriversAnnotationData, 
    axis = T, 
    axis_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1), labels = c("", "25%", "50%", "75%", "100%  "),
      gp = gpar(fontsize = 5),
      labels_rot = "0",
      side = "top"),
    ylim = c(0,1), 
    gp = gpar(fill = tsgDriverColours, col = 0), 
    border = F,
    bar_width = 0.9), 
  width = unit(2, "cm"),
  show_annotation_name = F,
  annotation_name_gp  = gpar(fontsize = 5),
  annotation_name_rot = 0
)

pTSG = tsgHeatmap + tsgSamplesAnnotation + tsgDriversAnnotation + tsgBiallelicAnnotation
pTSG

pdf("~/hmf/RPlot/Figure 3a.pdf", width = 18.3/2.54, height = 5/2.54)
draw(pTSG, padding = unit(c(0, 5, 0, 5), "mm"), auto_adjust = F)
dev.off()

pdf("~/hmf/RPlot/Figure 3b.pdf", width = 18.3/2.54, height = 5/2.54)
draw(pOnco, padding = unit(c(0, 5, 0, 5), "mm"), auto_adjust = F)
dev.off()

p_tsg_grob = grid.grabExpr(draw(pTSG, padding = unit(c(0, 5, 0, 5), "mm"), ht_gap = unit(c(2, 2, 2), "mm"), auto_adjust = F))
p_onco_grob = grid.grabExpr(draw(pOnco, padding = unit(c(0, 5, 0, 5), "mm"), ht_gap = unit(c(2, 2, 6), "mm"), auto_adjust = F))
p_germline_gro= grid.grabExpr(draw(pGermline, padding = unit(c(0, 5, 0, 5), "mm"), ht_gap = unit(c(2, 2, 2), "mm"), auto_adjust = F))

p_final = plot_grid(p_onco_grob, p_tsg_grob, p_germline_gro, ncol = 1, rel_heights = c(1, 1, 0.95), label_size = 8, labels = "auto", label_y = c(0.82,0.82, 0.82))
#p_final
ggplot2::ggsave("~/hmf/RPlot/Figure 3.pdf", p_final, width = 183, height = 120, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/RPlot/Figure 3.png", p_final, width = 183, height = 210, units = "mm", dpi = 300)


#p_final
#save_plot("~/hmf/RPlot/Figure 3 - Driver Heatmap.png", p_final, base_width = 14, base_height = 19)



#################
### Drivers per sample by driver Type
#### TO DO: Make colours consistent
head(hpcDriversByGene)
hpcDriversByGene = hpcDriversByGene %>% left_join(hpcCancerTypeCounts %>% select(cancerType, samples = N), by = "cancerType")
hpcDriversByGene = hpcDriversByGene %>%  mutate(driverRate = driverLikelihood / samples)

tidyDriversByCancerType = hpcDriversByGene %>% group_by(cancerType, driver, type) %>% summarise(driverLikelihood = sum(driverRate))
tidyDriversByCancerTypeLevels = hpcDriversByGene %>% group_by(cancerType) %>% summarise(driverLikelihood = sum(driverRate)) %>% arrange(-driverLikelihood)
tidyDriversByCancerType$cancerType = factor(tidyDriversByCancerType$cancerType, levels= tidyDriversByCancerTypeLevels$cancerType)

ggplot(data=tidyDriversByCancerType, aes(x = cancerType, y = driverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") +
  xlab("Cancer Type") + ylab("Mean Driver Count Per Sample") +
  scale_color_manual(values=simplifiedDriverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position="bottom",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) +
  coord_flip()# + facet_wrap(~type)
