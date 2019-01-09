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

tsgDrivers = c("Del","FragileDel","Indel","Missense","Multihit","Nonsense","Splice")
oncoDrivers = c("Amp","Del","Fusion","Indel","Missense","Promoter")
germlineDrivers = c("Del", "Indel","Missense", "Multihit","Nonsense", "Splice","Synonymous")
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

sortedTsgGenes = hpcDriversByGene %>% filter(type == 'TSG') %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(25, n) %>% arrange(-n)
sortedOncoGenes = hpcDriversByGene %>% filter(type == 'ONCO') %>% group_by(gene) %>% summarise(n = sum(driverLikelihood)) %>% ungroup() %>% top_n(25, n) %>% arrange(-n)
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
    mutate(proportion = 100 * round(n / sum(hpcCancerTypeCounts$N), 3)) %>% select(-n) 
  rownames(annotationData1) <- annotationData1$gene
  annotationData1$gene <- NULL
  annotationData1$non <- 0
  return (annotationData1)
}

driver_annotation <- function(sourceData) {
  annotationData2 = sourceData %>% 
    group_by(gene, driver) %>% summarise(n = sum(driverLikelihood)) %>% 
    group_by(gene) %>% mutate(total = sum(n)) %>% ungroup() %>%
    mutate(proportion = 100 * n / total) %>%
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
    mutate(biallelicPercentage = 100 * biallelic / (biallelic + notBiallelic) ) %>%  
    select(gene, biallelicPercentage) %>%
    mutate(nonBiallelicPercentage = 100 - biallelicPercentage)
  
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
    mutate(lostPercentage = 100 * lost / (lost + notLost) ) %>%  
    select(gene, lostPercentage) %>%
    mutate(notLostPercentage = 100 - lostPercentage)
  
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
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  row_order = sortedGermlineGenes$gene,
  column_order = sortedCancerTypes$cancerType,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  col = colorRamp2(c(0, 0.33, 0.66, 0.1), c("#f7fcf5","#bae4b3","#74c476","#238b45")), 
  cell_fun = function(j, i, x, y, w, h, col) {
    myColor = "black"
    if (germlineMat[i, j] > 0.15) {
      myColor = "white"
    }
    grid.text(heat_map_text(germlineMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)



germlineSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(baseline = "min", germlineSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,3), gp = gpar(fill = "#bc80bd"), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)

germlineDriversAnnotation = rowAnnotation(
  `% Drivers` = row_anno_barplot(germlineDriversAnnotationData, axis = T, axis_side = "top", ylim = c(0,100), gp = gpar(fill = germlineDriverColours), border = F), 
  width = unit(3, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)


germlineWildTypeLostAnnotationColours = setNames(c("#ccebc5","white"), c("lostPercentage", "notLostPercentage"))
germlineWildTypeLostAnnotation = rowAnnotation(
  show_annotation_name = T,
  annotation_name_rot = 0,
  annotation_name_gp  = gpar(fontsize = 7),
  `% Wild Type Lost` = row_anno_barplot(germlineWildtypeLostData, axis = T, axis_side = "top", ylim = c(0,100), gp = gpar(fill = germlineWildTypeLostAnnotationColours), border = F), 
  width = unit(2, "cm")
)


pGermline = germlineHeatmap + germlineSamplesAnnotation + germlineDriversAnnotation + germlineWildTypeLostAnnotation
pGermline

oncoMat = (data.matrix(oncoHeatmapData))
oncoHeatmap = Heatmap(
  column_title = "  ",
  oncoMat, 
  row_names_gp = gpar(fontsize = 9),
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
    grid.text(heat_map_text(oncoMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

oncoSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(oncoSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,20), gp = gpar(fill = "#bc80bd"), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)

oncoDriversAnnotation = rowAnnotation(
  `% Drivers` = row_anno_barplot(oncoDriversAnnotationData, axis = T, axis_side = "top", ylim = c(0,100), gp = gpar(fill = oncoDriverColours), border = F), 
  width = unit(3, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)

oncoDriversAnnotationIndex = rowAnnotation(
  df = data.frame(Driver = simplifiedDrivers),
  col = list(Driver = simplifiedDriverColours),
  width = unit(0, "cm"),
  annotation_legend_param = list(title = "", labels_gp = gpar(fontsize = 6.3))
)

pOnco = oncoHeatmap + oncoSamplesAnnotation + oncoDriversAnnotation + oncoDriversAnnotationIndex
pOnco

tsgMat = (data.matrix(tsgHeatmapData))
tsgHeatmap = Heatmap(
  column_title = "  ",
  tsgMat, 
  row_names_gp = gpar(fontsize = 9),
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
    grid.text(heat_map_text(tsgMat[i, j]), x, y, gp=gpar(fontsize=8, col = myColor))
  }
)

tsgSamplesAnnotation = rowAnnotation(
  `% Samples` = row_anno_barplot(tsgSamplesAnnotationData, axis = T, axis_side = "top", ylim = c(0,60), gp = gpar(fill = "#bc80bd"), border = F), 
  width = unit(2, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)

tsgBiallelicAnnotationColours = setNames(c("#ccebc5","white"), c("biallelicPercentage", "nonBiallelicPercentage"))
tsgBiallelicAnnotation = rowAnnotation(
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0,
  `% Biallelic` = row_anno_barplot(tsgBiallelicAnnotationData, axis = T, axis_side = "top", ylim = c(0,100), gp = gpar(fill = tsgBiallelicAnnotationColours), border = F), 
  width = unit(2, "cm")
)

tsgDriversAnnotation = rowAnnotation(
  `% Drivers` = row_anno_barplot(tsgDriversAnnotationData, axis = T, axis_side = "top", ylim = c(0,100), gp = gpar(fill = tsgDriverColours), border = F), 
  width = unit(3, "cm"),
  show_annotation_name = T,
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0
)

pTSG = tsgHeatmap + tsgSamplesAnnotation + tsgDriversAnnotation + tsgBiallelicAnnotation
pTSG

p_tsg_grob = grid.grabExpr(draw(pTSG))
p_onco_grob = grid.grabExpr(draw(pOnco))
p_germline_gro= grid.grabExpr(draw(pGermline))

p_final = plot_grid(p_tsg_grob, p_onco_grob, p_germline_gro, labels = "AUTO", ncol = 1, rel_heights = c(1, 1, 0.6))
p_final
save_plot("~/hmf/RPlot/Figure 3 - Driver Map.png", p_final, base_width = 14, base_height = 19)



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
