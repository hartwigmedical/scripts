library(dplyr)
library(tidyr)
library(ggplot2)

load("~/hmf/RData/output/driversByGene.RData")
load("~/hmf/RData/input/cohortByPrimaryTumorLocation.RData")
load("~/hmf/RData/input/PrimaryTumorLocationColours.RData")

driversByGene$impact = ifelse(driversByGene$multihit, "MultiHit", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$promoter), "Promoter", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$amp), "Amp", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$del), "Del", driversByGene$impact)

driversByGene = driversByGene %>% left_join(cohortByPrimaryTumorLocation %>% select(primaryTumorLocation, samples = N), by = "primaryTumorLocation")
driversByGene = driversByGene %>%  mutate(driverRate = driver / samples)

tidyDriversByCancerType = driversByGene %>% group_by(primaryTumorLocation, impact, type) %>% summarise(driver = sum(driverRate))
tidyDriversByCancerTypeLevels = driversByGene %>% group_by(primaryTumorLocation) %>% summarise(driver = sum(driverRate)) %>% arrange(-driver)
tidyDriversByCancerType$primaryTumorLocation = factor(tidyDriversByCancerType$primaryTumorLocation, levels= tidyDriversByCancerTypeLevels$primaryTumorLocation)

ggplot(data=tidyDriversByCancerType, aes(x = primaryTumorLocation, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers rate per cancer type") + xlab("Primary Tumor Location") + ylab("Driver Rate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), legend.position="bottom") + facet_wrap(~type)

tidyDriversByGene = driversByGene %>% group_by(gene, primaryTumorLocation, type, impact) %>% summarise(driver = sum(driver))
tidyDriversByGeneLevels = driversByGene %>% group_by(gene, type) %>% summarise(driver = sum(driver)) %>% arrange(-driver)
tidyDriversByGeneLevels = tidyDriversByGeneLevels[1:50, ]
tidyDriversByGene = tidyDriversByGene %>% filter(gene %in% tidyDriversByGeneLevels$gene) %>% ungroup() %>%
  mutate(gene = factor(gene, levels= tidyDriversByGeneLevels$gene)) 

typeColour <- ifelse(tidyDriversByGeneLevels$type == "ONCO", "red", "blue")
ggplot(data=tidyDriversByGene %>% filter(!is.na(primaryTumorLocation)), aes(x = gene, y = driver)) +
  geom_bar(aes(fill = primaryTumorLocation), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= primaryTumorLocationColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")

