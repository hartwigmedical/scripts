library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

load("~/hmf/RData/processed/driversByGene.RData")
load('~/hmf/RData/reference/hpcCancerTypeCounts.RData')
load("~/hmf/RData/reference/CancerTypeColours.RData")

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")

drivers = sort(unique(driversByGene$driver))
driverColours = setNames(cosmicSignatureColours[1:length(drivers)], drivers)

driversByGene = driversByGene %>% filter(!is.na(cancerType))

tidyDriversByCancerType = driversByGene %>% 
  left_join(hpcCancerTypeCounts %>% select(cancerType, samples = N), by = "cancerType") %>% 
  mutate(driversPerSample = driverLikelihood / samples) %>%
  group_by(cancerType, driver, type) %>% 
  summarise(driversPerSample = sum(driversPerSample))

tidyDriversByCancerTypeLevels = tidyDriversByCancerType %>% group_by(cancerType) %>% summarise(driversPerSample = sum(driversPerSample)) %>% arrange(-driversPerSample)
tidyDriversByCancerType$cancerType = factor(tidyDriversByCancerType$cancerType, levels=tidyDriversByCancerTypeLevels$cancerType)

ggplot(data=tidyDriversByCancerType, aes(x = cancerType, y = driversPerSample)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Drivers per sample by cancer type") + xlab("Primary Tumor Location") + ylab("Driver Per Sample") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position="bottom") 

#+ facet_wrap(~type)

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
  ggtitle("Absolute drivers per cancer type") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Relative drivers per cancer type") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Absolute driver types") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = relativeDriverLikelihood)) +
  geom_bar(aes(fill = driver), stat = "identity") + 
  ggtitle("Relative driver types") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= driverColours) +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour$colour), legend.position="bottom")


biallelicData = driversByGene %>% 
  mutate(biallelic = ifelse(driver %in% c("Multihit","Amp","Del"), T, biallelic)) %>%
  mutate(biallelic = ifelse(is.na(biallelic), F, biallelic)) %>%
  filter(type == 'TSG', !is.na(biallelic)) %>% ungroup() %>% group_by(gene, biallelic) %>% 
  summarise(driverLikelihood = sum(driverLikelihood)) %>% ungroup() %>% group_by(gene) %>%  mutate(biallelicPercentage = driverLikelihood / sum(driverLikelihood)) %>% 
  arrange(-biallelic, -biallelicPercentage)
biallelicDataLevels = unique(biallelicData$gene)
biallelicData$gene = factor(biallelicData$gene, biallelicDataLevels)


biallelicData = biallelicData %>% filter(gene %in% biallelicDataLevels[1:100])
biallelicData = biallelicData %>% filter(gene %in% tidyDriversByGene$gene)


ggplot(data=biallelicData, aes(x = gene, y = biallelicPercentage)) +
  geom_bar(aes(fill = biallelic), stat = "identity") + 
  ggtitle("Proportion of biallelic to non-biallelic drivers") + xlab("Gene") + ylab("Biallelic Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "blue"), legend.position="bottom")
