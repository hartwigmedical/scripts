library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

load('~/hmf/RData/reference/hpcCancerTypeCounts.RData')
load("~/hmf/RData/reference/CancerTypeColours.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")
driversByGene = hpcDriversByGene %>% filter(!is.na(cancerType), driverLikelihood > 0)

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")

drivers = sort(unique(driversByGene$driver))
driverColours = setNames(cosmicSignatureColours[1:length(drivers)], drivers)

tidyDriversByCancerType = driversByGene %>% 
  left_join(hpcCancerTypeCounts %>% select(cancerType, samples = N), by = "cancerType") %>% 
  mutate(driversPerSample = driverLikelihood / samples) %>%
  group_by(cancerType, driver, type) %>% 
  summarise(drivers = sum(driverLikelihood), driversPerSample = sum(driversPerSample))

tidyDriversByCancerTypeLevels = tidyDriversByCancerType %>% group_by(cancerType) %>% summarise(drivers = sum(drivers)) %>% arrange(-drivers)
tidyDriversByCancerType$cancerType = factor(tidyDriversByCancerType$cancerType, levels=tidyDriversByCancerTypeLevels$cancerType)

ggplot(data=tidyDriversByCancerType , aes(x = cancerType, y = drivers)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Driver Count by Cancer Type") + xlab("Primary Tumor Location") + ylab("Driver") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position="bottom")

tidyDriversByCancerTypeLevels = tidyDriversByCancerType %>% group_by(cancerType) %>% summarise(driversPerSample = sum(driversPerSample)) %>% arrange(-driversPerSample)
tidyDriversByCancerType$cancerType = factor(tidyDriversByCancerType$cancerType, levels=tidyDriversByCancerTypeLevels$cancerType)

ggplot(data=tidyDriversByCancerType , aes(x = cancerType, y = driversPerSample)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Average driver count by Cancer Type") + xlab("Primary Tumor Location") + ylab("Driver Per Sample") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position="bottom")

ggplot(data=tidyDriversByCancerType %>% filter(type != 'FUSION'), aes(x = cancerType, y = driversPerSample)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Average driver count by Cancer Type") + xlab("Primary Tumor Location") + ylab("Driver Per Sample") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position="bottom") + facet_grid(~type) 




tidyDriversByGene = driversByGene %>% filter(!is.na(cancerType)) %>% group_by(gene, cancerType, type, driver) %>% summarise(driverLikelihood = sum(driverLikelihood))
tidyDriversByGeneLevels = driversByGene %>% group_by(gene) %>% summarise(driverLikelihood = sum(driverLikelihood)) %>% arrange(-driverLikelihood)
tidyDriversByGeneLevels = tidyDriversByGeneLevels[1:50, ]
tidyDriversByGene = tidyDriversByGene %>% filter(gene %in% tidyDriversByGeneLevels$gene) %>% ungroup() %>%
  mutate(gene = factor(gene, levels= tidyDriversByGeneLevels$gene)) %>% group_by(gene) %>% mutate(relativeDriverLikelihood = driverLikelihood / sum(driverLikelihood)) %>%
  arrange(gene)

typeColour <- tidyDriversByGene %>% ungroup() %>% group_by(gene, type) %>% summarise(n = n()) %>% top_n(1, n) %>% filter(type != 'FUSION') %>%
  mutate(colour = ifelse(type == 'ONCO', 'red', 'blue'))

ggplot(data=tidyDriversByGene, aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Driver Count by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom") 

ggplot(data=tidyDriversByGene %>% filter(type == 'ONCO'), aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "red"), legend.position="bottom") 

ggplot(data=tidyDriversByGene %>% filter(type == 'TSG'), aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "blue"), legend.position="bottom") 


ggplot(data=tidyDriversByGene, aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Driver Count by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")


ggplot(data=tidyDriversByGene %>% filter(type == 'ONCO'), aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "red"), legend.position="bottom")

ggplot(data=tidyDriversByGene %>% filter(type == 'TSG'), aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Relative Drivers by Gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "blue"), legend.position="bottom")


topTsgGenes = driversByGene %>% filter(type == 'TSG') %>% group_by(gene)  %>% summarise(driverLikelihood = sum(driverLikelihood)) %>% top_n(50, driverLikelihood)
biallelicData = driversByGene %>% 
  mutate(biallelic = ifelse(driver %in% c("Multihit","Amp","Del"), T, biallelic)) %>%
  mutate(biallelic = ifelse(is.na(biallelic), F, biallelic)) %>%
  filter(type == 'TSG', driverLikelihood > 0, gene %in% topTsgGenes$gene) %>%
  ungroup() %>% 
  group_by(gene, biallelic) %>% 
  summarise(n = sum(driverLikelihood)) %>% 
  mutate(biallelic = ifelse(biallelic, "biallelic", "notBiallelic")) %>%
  spread(biallelic, n, fill = 0) %>%
  mutate(biallelicPercentage = biallelic / (biallelic + notBiallelic) ) %>%  
  arrange(-biallelicPercentage, -biallelic)
  #arrange(-biallelic, -biallelicPercentage)
  
biallelicDataLevels = unique(biallelicData$gene)
biallelicData$gene = factor(biallelicData$gene, biallelicDataLevels)

ggplot(data=biallelicData, aes(x = gene, y = biallelicPercentage)) +
  geom_bar(aes(fill = biallelic),  stat = "identity") + 
  ggtitle("Proportion of biallelic to non-biallelic drivers") + xlab("Gene") + ylab("Biallelic Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), legend.position="bottom")

topOncoGenes = driversByGene %>% filter(type == 'ONCO') %>% group_by(gene)  %>% summarise(driverLikelihood = sum(driverLikelihood)) %>% top_n(50, driverLikelihood)
hotspotData = driversByGene %>% 
  filter(type == 'ONCO', driverLikelihood > 0, gene %in% topOncoGenes$gene) %>%
  mutate(hotspot = ifelse(is.na(hotspot), F, hotspot)) %>%
  ungroup() %>% 
  group_by(gene, hotspot) %>% 
  summarise(n = sum(driverLikelihood)) %>% 
  mutate(hotspot = ifelse(hotspot, "hotspot", "notHotspot")) %>%
  spread(hotspot, n, fill = 0) %>%
  mutate(hotspotPercentage = hotspot / (hotspot + notHotspot) ) %>%  
  filter(hotspot > 0) %>%
  arrange(-hotspotPercentage, -hotspot)
  
hotspotDataLevels = unique(hotspotData$gene)
hotspotData$gene = factor(hotspotData$gene, hotspotDataLevels)

ggplot(data=hotspotData, aes(x = gene, y = hotspotPercentage)) +
  geom_bar(aes(fill = hotspot),  stat = "identity") + 
  ggtitle("Proportion of hotspot drivers") + xlab("Gene") + ylab("Hotspot Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), legend.position="bottom")


################################# ONCO PRINT ########################
library(ComplexHeatmap)

load("~/hmf/RData/processed/hpcDriversByGene.RData")
driversByGene = hpcDriversByGene %>% filter(!is.na(cancerType), driverLikelihood > 0)

data = sample
driverFilter = "Amp"
rm(data, driverFilter)
as_onco_print_matrix <- function(data, driverFilter) {
  result = data %>%
    mutate(field = (driver == driverFilter)) %>% 
    group_by(sampleId, gene, field) %>%
    summarise() %>% 
    spread(sampleId, field, fill = 0) %>% 
    ungroup()
  rownames(result) <- result$gene 
  result = result %>% select(-gene)
  return (data.matrix(result))
}

sampleGenes = driversByGene %>% filter(driverLikelihood > 0.5) %>% group_by(gene) %>% summarise(driverLikelihood = sum(driverLikelihood)) %>% top_n(50, driverLikelihood)
sample = driversByGene %>% filter(driverLikelihood > 0.5, gene %in% sampleGenes$gene)

mat_list = list()
for (driverType in drivers) {
  cat(driverType, "\n")
  mat_list[[driverType]] = as_onco_print_matrix(sample, driverType)
}

alter_fun = list()
for (driverType in drivers) {
  cat(driverType, "\n")
  color = driverColours[driverType]
  alter_fun[[driverType]] = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = bquote(driverColours[driverType]), col = NA))
}


generate_alter_fun <- function(drivers) {
  alter_fun = list()
  for (driverType in drivers) {
    cat(driverType, "\n")
    color = driverColours[driverType]
    alter_fun[[driverType]] = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = color, col = NA))
  }
  return (alter_fun)
}
alter_fun = generate_alter_fun(drivers)

rm(color)

jno = (driverColours[driverType])

for (driverType in drivers) {
  quote(driverColours[driverType])
  eval(driverColours[driverType])
}


#mat_list = list(snv = matrix(c(1, 0, 1, 1, 1, 0, 0, 1, 1), nrow = 3),
#                indel = matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0), nrow = 3))
#rownames(mat_list$snv) = rownames(mat_list$indel) = c("g1", "g2", "g3")
#colnames(mat_list$snv) = colnames(mat_list$indel) = c("s1", "s2", "s3")
#mat_list


oncoPrint(mat_list, alter_fun = alter_fun, col = driverColours)

drivers


driverColours["Missense"]
date()

oncoPrint(mat_list,
  alter_fun = list(
    "Fusion-Intragenic" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = driverColours["Fusion-Intragenic"], col = NA)),
    "Fusion-Coding" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = driverColours["Fusion-Coding"], col = NA)),
    "Fusion-UTR" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = driverColours["Fusion-UTR"], col = NA)),
    "Del" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Del"], col = NA)),
    "FragileDel" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["FragileDel"], col = NA)),
    "Multihit" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Multihit"], col = NA)),
    "Promoter" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Promoter"], col = NA)),
    "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Nonsense"], col = NA)),
    "Splice" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Splice"], col = NA)),
    "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Missense"], col = NA)),
    "Inframe" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Inframe"], col = NA)),
    "Indel" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = driverColours["Indel"], col = NA)),
    "Amp" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = driverColours["Amp"], col = NA))
  ), col = driverColours)

