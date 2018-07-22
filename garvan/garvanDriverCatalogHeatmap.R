library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)

outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"

referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
plotDir = paste0(outputDir, "plots/")

load(paste0(processedDir, "hpcDriverCatalog.RData"))
load(paste0(processedDir, "driverGenes.RData"))
load(paste0(processedDir, "hpcCancerTypeCounts.RData"))
load(paste0(processedDir, "simplifiedDrivers.RData"))

tsgDrivers = c("Amp","Del","Indel","Missense","Multihit","Nonsense","Splice")
oncoDrivers = c("Amp","Del","Fusion","Indel","Missense","Promoter")
tsgDriverColours = simplifiedDriverColours[tsgDrivers]
oncoDriverColours = simplifiedDriverColours[oncoDrivers]

hpcDriversByGene = hpcDriversByGene  %>% 
  filter(driverLikelihood > 0, !cancerType %in% c('EPD','HEMNOS','HCC')) %>%
  #filter(driverLikelihood > 0) %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
    driver = ifelse(driver == "Inframe", "Indel", driver),
    hotspot = ifelse(driver == "Promoter", T, hotspot)
  ) %>%
  ungroup()

sortedTsgGenes = hpcDriversByGene %>% filter(type == 'TSG') %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(30, n) %>% arrange(-n)
sortedOncoGenes = hpcDriversByGene %>% filter(type == 'ONCO', !gene %in% c('CDKN2B', 'CDKN2B-AS1')) %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(30, n) %>% arrange(-n)
sortedCancerTypes = hpcDriversByGene %>% group_by(cancerType) %>% summarise(n = sum(driverLikelihood)) %>% arrange(-n)

oncoDriversByGene = hpcDriversByGene %>% filter(type == 'ONCO', gene %in% sortedOncoGenes$gene) %>% mutate(driver = factor(driver, simplifiedDrivers))
tsgDriversByGene = hpcDriversByGene %>% filter(type == 'TSG', gene %in% sortedTsgGenes$gene)  %>% mutate(driver = factor(driver, simplifiedDrivers))

main_heatmap_data <- function(sourceData, cancerTypeSamples) {
  heatmapData = sourceData %>%
    group_by(cancerType, gene) %>% summarise(n = sum(driverLikelihood)) %>%
    left_join(hpcCancerTypeCounts %>% select(cancerType, cancerTypeSamples = N), by = "cancerType") %>%
    mutate(proportion = n / cancerTypeSamples) %>%
    select(-n, -cancerTypeSamples) %>%
    spread(cancerType, proportion, fill = 0)
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
    mutate(biallelicPercentage = biallelic / (biallelic + notBiallelic) ) %>%  
    select(gene, biallelicPercentage) %>%
    mutate(nonBiallelicPercentage = 1 - biallelicPercentage)
  
  rownames(biallelicData) <- biallelicData$gene
  biallelicData$gene <- NULL
  return (biallelicData)
}


tsgHeatmapData = main_heatmap_data(tsgDriversByGene, cancerTypeSamples)
tsgSamplesAnnotationData = samples_annotation(tsgDriversByGene)
tsgDriversAnnotationData = driver_annotation(tsgDriversByGene)
tsgBiallelicAnnotationData = biallelic_annotation(tsgDriversByGene)

oncoHeatmapData = main_heatmap_data(oncoDriversByGene, cancerTypeSamples)
oncoSamplesAnnotationData = samples_annotation(oncoDriversByGene)
oncoDriversAnnotationData = driver_annotation(oncoDriversByGene)


heat_map_text <- function(value) {
  value = round(100*value, 1)
  if (value > 0) {
    return (sprintf(fmt='%.1f', value))
  }
  return ("")
}

oncoMat = (data.matrix(oncoHeatmapData))
oncoHeatmap = Heatmap(
  column_title = "  ",
  oncoMat, 
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
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
    grid.text(heat_map_text(oncoMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

oncoSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(oncoSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,0.2), gp = gpar(fill = "#bc80bd"), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_rot = 0
)

oncoDriversAnnotation = rowAnnotation(
  Drivers = row_anno_barplot(oncoDriversAnnotationData, axis = T, axis_side = "top", ylim = c(0,1), gp = gpar(fill = oncoDriverColours), border = F), 
  width = unit(3, "cm"),
  show_annotation_name = T,
  annotation_name_rot = 0
)

oncoDriversAnnotationIndex = rowAnnotation(
  df = data.frame(Driver = simplifiedDrivers),
  col = list(Driver = simplifiedDriverColours),
  width = unit(0, "cm"),
  annotation_legend_param = list(title = "")
)

pOnco = oncoHeatmap + oncoSamplesAnnotation + oncoDriversAnnotation + oncoDriversAnnotationIndex
png(file = paste0(plotDir, "OncoHeatmap.png"), width = 10, height = 7, units = "in", res = 217)
pOnco
dev.off()


tsgMat = (data.matrix(tsgHeatmapData))
tsgHeatmap = Heatmap(
  column_title = "  ",
  tsgMat, 
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  row_order = sortedTsgGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("#f7fbff","#bdd7e7","#6baed6","#2171b5")),
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (tsgMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(heat_map_text(tsgMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

tsgSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(tsgSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,0.2), gp = gpar(fill = "#bc80bd"), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_rot = 0
)

tsgBiallelicAnnotationColours = setNames(c("#ccebc5","white"), c("biallelicPercentage", "nonBiallelicPercentage"))
tsgBiallelicAnnotation = rowAnnotation(
  show_annotation_name = T,
  annotation_name_rot = 0,
  `% Biallelic` = row_anno_barplot(tsgBiallelicAnnotationData, axis = T, axis_side = "top", ylim = c(0,1), gp = gpar(fill = tsgBiallelicAnnotationColours), border = F), 
  width = unit(2, "cm")
)

tsgDriversAnnotation = rowAnnotation(
  Drivers = row_anno_barplot(tsgDriversAnnotationData, axis = T, axis_side = "top", ylim = c(0,1), gp = gpar(fill = tsgDriverColours), border = F), 
  width = unit(3, "cm"),
  show_annotation_name = T,
  annotation_name_rot = 0
)

tsgDriversAnnotationIndex = rowAnnotation(
  df = data.frame(Driver = tsgDrivers),
  col = list(Driver = tsgDriverColours),
  width = unit(0, "cm"),
  annotation_legend_param = list(title = "")
  
)

pTSG = tsgHeatmap + tsgSamplesAnnotation + tsgDriversAnnotation + tsgBiallelicAnnotation
pTSG

tsgHeatmap + tsgSamplesAnnotation + tsgBiallelicAnnotation + tsgDriversAnnotation + tsgDriversAnnotationIndex

png(file = paste0(plotDir, "TsgHeatmap.png"), width = 10, height = 7, units = "in", res = 217)
pTSG
dev.off()

