detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_bw())



load("~/hmf/RData/processed/driverGenes.RData")
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary[is.na(highestPurityCohortSummary)] <- 0

simplifiedDrivers = c("Amp","Del","FragileDel","Fusion","Indel","Missense","Multihit","Nonsense","Promoter","Splice") 
simplifiedDriverColours = c("#fb8072","#bc80bd","#bebada", "#fdb462","#80b1d3","#8dd3c7","#b3de69","#fccde5","#ffffb3","#d9d9d9")
simplifiedDriverColours = setNames(simplifiedDriverColours, simplifiedDrivers)


########## Driver rates by type in WGD vs non WGD
wgdDriverSampleCount = highestPurityCohortSummary %>% group_by(WGD) %>% summarise(total=n())
wgdDriverRates = hpcDriversByGene %>%
  filter(driverLikelihood > 0) %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
    driver = ifelse(driver %in% c("Frameshift", "Inframe"), "Indel", driver),
    driver = ifelse(driver %in% c("FragileDel"), "Del", driver),
    driver = ifelse(driver %in% c("Missense","Nonsense","Splice"), "SNV/MNV", driver)
    ) %>%
  group_by(sampleId, driver) %>%
  summarise(n = sum(driverLikelihood)) %>%
  #count() %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, WGD), by = "sampleId") %>%
  group_by(driver, WGD) %>%
  summarise(n = sum(n)) %>%
  left_join(wgdDriverSampleCount, by = "WGD") %>%
  mutate(n = n / total) %>%
  ungroup()
wgdDriverLevels = wgdDriverRates %>% filter(WGD) %>% arrange(-n)
wgdDriverRates = mutate(wgdDriverRates, driver = factor(driver, wgdDriverLevels$driver))

p1 = ggplot(data = wgdDriverRates, aes(x = driver, y = n)) +
  geom_bar(aes(fill = WGD), stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#f1a340", "#998ec3")) +
  ggtitle("Drivers Per Sample by WGD") + 
  xlab("Driver") + ylab("Drivers") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  scale_y_continuous(expand=c(0.01, 0.01))

plot_grid(p1, labels="AUTO")

#############  Drivers vs Mutational Load
mutationalLoad = highestPurityCohortSummary %>%
  filter(msiStatus == "MSS") %>%
  mutate(
    indel = CLONAL_INDEL + SUBCLONAL_INDEL + INCONSISTENT_INDEL,
    snv = CLONAL_SNP + SUBCLONAL_SNP + INCONSISTENT_SNP + CLONAL_MNP + SUBCLONAL_MNP + INCONSISTENT_MNP,
    sv = BND + DEL + DUP + INV + INS) %>%
  select(sampleId, cancerType, indel, snv, sv) 

indelDrivers =  hpcDriversByGene %>%   
  ungroup() %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(driver %in% c("Frameshift", "Inframe"), "Indel", driver),
    driver = ifelse(grepl("Indel", impact), "Indel", driver)
  ) %>%
  filter(driverLikelihood > 0, driver == "Indel") %>%
  group_by(sampleId, cancerType) %>% summarise(drivers = sum(driverLikelihood, na.rm = T))

snvDrivers =  hpcDriversByGene %>%   
  ungroup() %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(driver %in% c("Missense", "Nonsense","Splice"), "SNV", driver),
    driver = ifelse(grepl("Missense", impact) | grepl("Nonsense", impact) | grepl("Splice", impact), "SNV", driver)
  ) %>%
  filter(driverLikelihood > 0, driver == "SNV") %>%
  group_by(sampleId, cancerType) %>% summarise(drivers = sum(driverLikelihood, na.rm = T))

cnaDrivers =  hpcDriversByGene %>%   
  ungroup() %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(driver %in% c("Amp", "Del","FragileDel"), "CNA", driver)) %>%
  filter(driverLikelihood > 0, driver == "CNA") %>%
  group_by(sampleId, cancerType) %>% summarise(drivers = sum(driverLikelihood, na.rm = T))

indelData = left_join(mutationalLoad %>% select(sampleId, cancerType, indel), indelDrivers, by = c("sampleId","cancerType"))
indelData[is.na(indelData)] <- 0
indelDataGlobal = indelData %>% ungroup() %>% summarise(load = mean(indel, na.rm = T), drivers = mean(drivers, na.rm = T)) %>% mutate(gradient = drivers/load) 
indelData = indelData %>% group_by(cancerType) %>% summarise(load = mean(indel, na.rm = T), drivers = mean(drivers, na.rm = T))

svData = left_join(mutationalLoad %>% select(sampleId, cancerType, sv), cnaDrivers, by = c("sampleId","cancerType"))
svData[is.na(svData)] <- 0
svDataGlobal = svData %>% ungroup() %>% summarise(load = mean(sv, na.rm = T), drivers = mean(drivers, na.rm = T)) %>% mutate(gradient = drivers/load) 
svData = svData %>% group_by(cancerType) %>% summarise(load = mean(sv, na.rm = T), drivers = mean(drivers, na.rm = T))

snvData = left_join(mutationalLoad %>% select(sampleId, cancerType, snv), snvDrivers, by = c("sampleId","cancerType"))
snvData[is.na(snvData)] <- 0
snvDataGlobal = snvData %>% ungroup() %>% summarise(load = mean(snv, na.rm = T), drivers = mean(drivers, na.rm = T)) %>% mutate(gradient = drivers/load) 
snvData = snvData %>% group_by(cancerType) %>% summarise(load = mean(snv, na.rm = T), drivers = mean(drivers, na.rm = T))

p1 = ggplot(svData, aes(x=load, y=drivers, color=cancerType)) + 
  geom_point() +
  geom_segment(aes(x = 0, xend = 1000, y = 0, yend = 1000 * svDataGlobal$gradient), linetype = 2, color = "black") +
  geom_point(aes(x = svDataGlobal$load, y = svDataGlobal$drivers), shape = 25, size = 2, color= "blue", fill="blue") +
  coord_cartesian(xlim= c(0,500), ylim=c(0,5))+
  ylab("CNA Drivers") + xlab("SV Mutational Load") + ggtitle("Copy Number Alteration Drivers") + 
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)

p2 = ggplot(snvData, aes(x=load, y=drivers, color=cancerType)) + 
  geom_point() +
  geom_segment(aes(x = 0, xend = 100000, y = 0, yend = 100000 * snvDataGlobal$gradient), linetype = 2, color = "black") +
  geom_point(aes(x = snvDataGlobal$load, y = snvDataGlobal$drivers), shape = 25, size = 2, color= "blue", fill="blue") +
  coord_cartesian(xlim= c(0,82000), ylim=c(0,4))+
  ylab("SNV Drivers") + xlab("SNV Mutational Load") + ggtitle("SNV/MNV Drivers") + 
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)

p3 = ggplot(indelData, aes(x=load, y=drivers, color=cancerType)) + 
  geom_point() +
  geom_segment(aes(x = 0, xend = 10000, y = 0, yend = 10000 * indelDataGlobal$gradient), linetype = 2, color = "black") +
  geom_point(aes(x = indelDataGlobal$load, y = indelDataGlobal$drivers), shape = 25, size = 2, color= "blue", fill="blue") +
  coord_cartesian(xlim= c(0,3000), ylim=c(0,1.5))+
  ylab("Indel Drivers") + xlab("Indel Mutational Load") + ggtitle("Indel Drivers") +
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)

p4 = ggplot(indelData, aes(x=load, y=drivers, color=cancerType)) + geom_point() +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 3)) +
  scale_color_manual(values=cancerTypeColours)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(p4) 

plot_grid(p1,p2,p3,legend, ncol = 1, rel_heights = c(3,3,3,2), labels = c("A","B","C"))


#############  Driver Per CancerType

driverData = hpcDriversByGene %>%
  mutate(driver = as.character(driver), 
         driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
         driver = ifelse(driver %in% c("Frameshift","Inframe"), 'Indel', driver),
         driver = factor(driver, simplifiedDrivers)
  ) %>%
  group_by(driver, cancerType) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  left_join(hpcCancerTypeCounts %>% select(cancerType, N), by = "cancerType") %>%
  mutate(driversPerSample = driverLikelihood/ N) %>%
  arrange(-driversPerSample)
driverDataLevels = driverData %>% group_by(cancerType) %>% summarise(driversPerSample = sum(driversPerSample)) %>% arrange(-driversPerSample)
driverData = driverData %>% mutate(cancerType = factor(cancerType, driverDataLevels$cancerType))

driverViolinData = hpcDriversByGene %>%
  group_by(cancerType, sampleId) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  ungroup() %>%
  mutate(cancerType = factor(cancerType, driverDataLevels$cancerType))

p1 = ggplot(driverViolinData, aes(cancerType, driverLikelihood)) + 
  geom_violin(trim = F, color = "#6baed6", fill = "#6baed6", scale = "area") +
  #stat_summary(fun.y = mean, geom = "pointrange", color = "black") + 
  #stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "red")+
  geom_boxplot(width = 0.2, outlier.shape = NA, fill  = "#deebf7") +
  #stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", color = "black", size = 0.1)+
  xlab("Cancer Type") + ylab("Drivers") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  ggtitle("") + xlab("") + ylab("Drivers per sample") + 
  scale_y_continuous(expand=c(0.01, 0.01)) +
  theme(legend.position="none") +
  coord_flip() 

p2 = ggplot(driverData, aes(cancerType, driversPerSample)) +
  geom_bar(stat = "identity", aes(fill = driver)) +
  scale_fill_manual(values = simplifiedDriverColours) +
  ggtitle("") + xlab("") + ylab("Average drivers per sample") +  
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) +
  scale_y_continuous(expand=c(0.0, 0.0)) +
  coord_flip()

plot_grid(p1,p2, ncol = 2, rel_widths =  c(1,2), labels = "AUTO")



######################################################

#Violin - Rounded total  Indel Drivers vs median Indel count
indelCount = highestPurityCohortSummary %>%
  mutate(indel = CLONAL_INDEL + SUBCLONAL_INDEL + INCONSISTENT_INDEL) %>%
  select(sampleId, indel)

#TODO: Multihit with INDEL
indelDrivers = hpcDriversByGene %>%   
  ungroup() %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(driver %in% c("Frameshift", "Inframe"), "Indel", driver)
  ) %>%
  filter(driverLikelihood > 0, driver == "Indel") %>%
  select(sampleId, driver, impact, driverLikelihood) %>%
  group_by(sampleId) %>%
  summarise(driverLikelihood = round(sum(driverLikelihood))) %>%
  mutate(driverLikelihood = ifelse(driverLikelihood >= 4, "4+", driverLikelihood))

indelData = left_join(indelCount, indelDrivers, by = "sampleId")
indelData[is.na(indelData)] <- 0
indelData = mutate(indelData, driverLikelihood = factor(driverLikelihood))
indelData$indel <- ifelse(indelData$indel > 10000, 10000, indelData$indel)

str(indelData)

ggplot(indelData, aes(driverLikelihood, indel)) + 
  geom_violin(scale = "count", draw_quantiles = c(0.5)) +
  xlab("Driver Indels") + ylab("Total Indels") 


  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  scale_colour_manual(values=cancerTypeColours, guide=FALSE) +
  ylab("Age") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_blank()) + 
  facet_grid(~cancerType)


#Violin - Rounded Total SNV Drivers vs median SNV count
#Violin - Amps & Dels vs SV count
