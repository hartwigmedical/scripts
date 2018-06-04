library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

load("~/hmf/RData/processed/hpcDriversByGene.RData")
load("~/hmf/RData/processed/driverGenes.RData")
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')

hpcDriversByGene = hpcDriversByGene  %>% 
  filter(driverLikelihood > 0) %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
    driver = ifelse(driver == "Inframe", "Indel", driver),
    hotspot = ifelse(driver == "Promoter", T, hotspot)
  ) %>%
  ungroup()

sortedTsgGenes = hpcDriversByGene %>% filter(type == 'TSG') %>% group_by(gene) %>% count() %>% ungroup() %>% top_n(30, n) %>% arrange(-n)
sortedOncoGenes = hpcDriversByGene %>% filter(type == 'ONCO') %>% group_by(gene) %>% count() %>% ungroup() %>% top_n(30, n) %>% arrange(-n)
sortedCancerTypes = hpcDriversByGene %>% group_by(cancerType) %>% count() %>% arrange(-n)

oncoDrivers = c("Fusion", "Del","Promoter","Missense","Indel","Amp")
oncoDrivers = factor(oncoDrivers, levels = oncoDrivers)
oncoDriverColours = c("#fdb462","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
oncoDriverColours = setNames(oncoDriverColours, oncoDrivers)

tsgDrivers = c("Fusion","Del","FragileDel","Multihit","Nonsense","Splice","Missense","Indel")
tsgDrivers = factor(tsgDrivers, levels = tsgDrivers)
tsgDriverColours = c("#fdb462","#8dd3c7","#fccde5","#b3de69", "#80b1d3","#ffffb3","#bebada","#fb8072")
tsgDriverColours = setNames(tsgDriverColours, tsgDrivers)



oncoDriversByGene = hpcDriversByGene %>% filter(type == 'ONCO', gene %in% sortedOncoGenes$gene) %>% mutate(driver = factor(driver, oncoDrivers))
tsgDriversByGene = hpcDriversByGene %>% filter(type == 'TSG', gene %in% sortedTsgGenes$gene)  %>% mutate(driver = factor(driver, tsgDrivers))

main_heatmap_data <- function(sourceData, cancerTypeSamples) {
  heatmapData = sourceData %>%
    group_by(cancerType, gene) %>% count() %>%
    left_join(hpcCancerTypeCounts %>% select(cancerType, cancerTypeSamples = N), by = "cancerType") %>%
    #mutate(proportion = round(n / cancerTypeSamples, 3)) %>%
    mutate(proportion = n / cancerTypeSamples) %>%
    select(-n, -cancerTypeSamples) %>%
    spread(cancerType, proportion, fill = 0)
  rownames(heatmapData) <- heatmapData$gene
  heatmapData$gene <- NULL
  return (heatmapData)
}

samples_annotation <- function(sourceData, max = 1) {
  annotationData1 = sourceData %>% 
    group_by(gene) %>% count() %>% 
    mutate(proportion = round(n / sum(hpcCancerTypeCounts$N), 3)) %>% select(-n) %>%
    mutate(nonProportion = max - proportion)
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
tsgSamplesAnnotationData = samples_annotation(tsgDriversByGene, max = 0.6)
tsgDriversAnnotationData = driver_annotation(tsgDriversByGene)
tsgHotspotAnnotationData = hotspot_annotation(tsgDriversByGene)
tsgBiallelicAnnotationData = biallelic_annotation(tsgDriversByGene)
tsgBiallelicAnnotationData

oncoHeatmapData = main_heatmap_data(oncoDriversByGene, cancerTypeSamples)
oncoSamplesAnnotationData = samples_annotation(oncoDriversByGene, max = 0.2)
oncoDriversAnnotationData = driver_annotation(oncoDriversByGene)

oncoMat = (data.matrix(oncoHeatmapData))
oncoHeatmap = Heatmap(
  column_title = "  ",
  oncoMat, 
  row_order = sortedOncoGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("white", "pink", "red", "darkred")),
  #col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("white", "lightblue", "blue", "darkblue")),
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (oncoMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(sprintf(fmt='%.1f', 100*oncoMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

oncoSamplesAnnotationColours = setNames(c("red","white"), c("proportion", "nonProportion"))
oncoSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(oncoSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,0.2), gp = gpar(fill = oncoSamplesAnnotationColours), border = F), 
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
  df = data.frame(Driver = oncoDrivers),
  col = list(Driver = oncoDriverColours),
  width = unit(0, "cm"),
  annotation_legend_param = list(title = "")

)

oncoHeatmap + oncoSamplesAnnotation + oncoDriversAnnotation + oncoDriversAnnotationIndex



tsgMat = (data.matrix(tsgHeatmapData))
tsgHeatmap = Heatmap(
  column_title = "  ",
  tsgMat, 
  row_order = sortedTsgGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.1, 0.2, 0.6), c("white", "lightblue", "blue", "darkblue")),
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (tsgMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(sprintf(fmt='%.1f', 100*tsgMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

tsgSamplesAnnotationColours = setNames(c("blue","white"), c("proportion", "nonProportion"))
tsgSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(tsgSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,0.6), gp = gpar(fill = tsgSamplesAnnotationColours), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_rot = 0
)

tsgBiallelicAnnotationColours = setNames(c("blue","white"), c("biallelicPercentage", "nonBiallelicPercentage"))
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

tsgHeatmap + tsgSamplesAnnotation + tsgBiallelicAnnotation + tsgDriversAnnotation + tsgDriversAnnotationIndex