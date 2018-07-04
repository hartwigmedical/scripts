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
save(simplifiedDrivers, simplifiedDriverColours, file = "~/hmf/RData/reference/simplifiedDrivers.RData")


########################################### Clonality
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")

highestPurityCohortSummary$purityBucket = 
  cut(
    highestPurityCohortSummary$purity, 
    breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    labels = c("0-10%", "10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","          90-100%"))

clonalityFactor = c('SUBCLONAL','CLONAL')

clonalityLoad = highestPurityCohortSummary %>%
  filter(cancerType != 'Other') %>%  mutate(grouping = cancerType) %>% 
  mutate(total = TOTAL_INDEL + TOTAL_SNV + TOTAL_MNV, 
         subclonal = SUBCLONAL_INDEL + SUBCLONAL_SNV + SUBCLONAL_MNV) %>%
  group_by(cancerType, sampleId) %>%
  summarise(total = sum(total), SUBCLONAL = sum(subclonal), CLONAL = sum(total) - SUBCLONAL) %>%
  gather(clonality, value, CLONAL, SUBCLONAL) %>%
  mutate(percentage = value / total) %>%
  filter(clonality == 'CLONAL') %>%
  mutate(clonality = factor(clonality, clonalityFactor), type = 'All')

clonalityLoad2 = highestPurityCohortSummary %>%
  mutate(total = TOTAL_INDEL + TOTAL_SNV + TOTAL_MNV, 
         subclonal = SUBCLONAL_INDEL + SUBCLONAL_SNV + SUBCLONAL_MNV) %>%
  group_by(purityBucket, sampleId) %>%
  summarise(total = sum(total), SUBCLONAL = sum(subclonal), CLONAL = sum(total) - SUBCLONAL) %>%
  gather(clonality, value, CLONAL, SUBCLONAL) %>%
  mutate(percentage = value / total) %>%
  filter(clonality == 'CLONAL') %>%
  mutate(
    clonality = factor(clonality, clonalityFactor), 
    type = 'All')

samplesPerPurityBucket = clonalityLoad2 %>% group_by(purityBucket) %>% count()

clonalityDrivers2 = hpcDriversByGene %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, purityBucket), by = "sampleId") %>%
  filter(!is.na(clonality)) %>%
  mutate(clonality = ifelse(clonality == 'INCONSISTENT', 'CLONAL', clonality)) %>%
  group_by(purityBucket, clonality) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>% 
  spread(clonality, driverLikelihood) %>%
  mutate(subclonalPercentage = SUBCLONAL / (CLONAL + SUBCLONAL))
  
p1 = ggplot(samplesPerPurityBucket, aes(purityBucket, n)) + 
  geom_bar(fill = "#31a354", stat = "identity") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  ggtitle("") + xlab("Tumor Purity") + ylab("Count of samples") + 
  scale_y_continuous(expand=c(0.01, 0.01)) +
  theme(legend.position="none") +
  coord_flip() 

p2 = ggplot(clonalityLoad2, aes(purityBucket, percentage)) + 
  geom_violin(fill = "#31a354", scale = "width") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2) +
  ggtitle("") + xlab("") + ylab("% of variants clonal") + 
  scale_y_continuous(limits = c(0, 1), labels = percent, expand=c(0.01, 0.01)) +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  coord_flip() 

p3 = ggplot(data = clonalityDrivers2, aes(x = purityBucket, y = subclonalPercentage, width = 0.7)) +
  geom_bar(fill = "#31a354", stat = "identity") + 
  xlab("") + ylab("% of driver variants subclonal") + ggtitle("") + 
  scale_y_continuous(labels = percent, expand=c(0.01, 0.01), limits = c(0, 0.06)) +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  coord_flip()

plot_grid(p1,p2, p3, ncol = 3, labels="AUTO")


pdf(file='~/hmf/supplementaryVisualisaion.pdf', onefile=T, paper='A4') 

######### Hotspots
load(file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")

hotspotGenes = hpcDndsOncoDrivers %>% 
  filter(type == 'ONCO', !is.na(hotspot)) %>% 
  group_by(gene) %>% 
  summarise(driverLikelihood = sum(driverLikelihood)) %>% top_n(40, driverLikelihood) %>% arrange(driverLikelihood)

variantFactors = c("INDEL","MNV","SNV")
hotspotFactors = c("Hotspot","NearHotspot","NoHotspot")
hotspotData = hpcDndsOncoDrivers %>% 
  filter(type == 'ONCO', !is.na(hotspot), driverLikelihood > 0, gene %in% hotspotGenes$gene) %>%
  mutate(hotspot = ifelse(hotspot, "Hotspot", "NoHotspot"),
         hotspot = ifelse(nearHotspot, "NearHotspot", hotspot),
         hotspot = factor(hotspot, rev(hotspotFactors)),
         variant = factor(variant, rev(variantFactors)),
         impact = as.character(impact),
         impact = ifelse(impact == "Promoter", "Missense", impact),
         impact = ifelse(impact == "Inframe", "Inframe Indel", impact),
         impact = factor(impact, c("Missense", "Inframe Indel"))) %>%
  group_by(gene, impact, variant, hotspot) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  group_by(gene) %>%
  mutate(total = sum(driverLikelihood), percentage = driverLikelihood / total) %>%
  ungroup() %>%
  mutate(gene = factor(gene, hotspotGenes$gene))

p_hotspot = ggplot(data = hotspotData, aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = hotspot), stat = "identity") +
  scale_fill_manual(values = c("#fee0d2", "#fb6a4a", "#de2d26")) +
  ggtitle("Oncogene Hotspots") + 
  xlab("Gene") + ylab("Drivers")+ 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom",  strip.background = element_blank(), legend.title=element_blank()) +
  coord_flip()+ 
  facet_grid(~variant, scales = "free_x")

plot_grid(p_hotspot, labels="AUTO")


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
    indel =TOTAL_INDEL,
    snv = TOTAL_SNV + TOTAL_MNV,
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
driverData = driverData %>% 
  mutate(cancerType = factor(cancerType, driverDataLevels$cancerType)) %>%
  group_by(cancerType) %>% 
  mutate(percentage = driversPerSample / sum(driversPerSample))

driverViolinData = hpcDriversByGene %>%
  group_by(cancerType, sampleId) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  ungroup() %>%
  mutate(cancerType = factor(cancerType, driverDataLevels$cancerType))

p1 = ggplot(driverViolinData, aes(cancerType, driverLikelihood)) + 
  geom_violin( fill = "#6baed6", scale = "area") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2) +
  xlab("Cancer Type") + ylab("Drivers") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  ggtitle("") + xlab("") + ylab("Drivers per sample") + 
  scale_y_continuous(expand=c(0.01, 0.01)) +
  theme(legend.position="none") +
  coord_flip() 

p2 = ggplot(driverData, aes(cancerType, percentage)) +
  geom_bar(stat = "identity", aes(fill = driver), width=0.7) +
  scale_fill_manual(values = simplifiedDriverColours) +
  ggtitle("") + xlab("") + ylab("Average drivers per sample") +  
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) +
  scale_y_continuous(labels = percent, expand=c(0.0, 0.0), limits = c(0, 1)) +
  coord_flip()

plot_grid(p1,p2, ncol = 2, rel_widths =  c(2,3), labels = "AUTO")

dev.off()

##################################################### MSI ##################################################### 

load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
msiData = highestPurityCohortSummary %>% select(sampleId, cancerType, msiScore, msiStatus)
msiRelativeData = msiData %>% group_by(cancerType, msiStatus) %>% count() %>% group_by(cancerType) %>% mutate(percentage = n / sum(n)) 
msiRelativeData = msiData %>% filter(cancerType != 'Other') %>% group_by(cancerType, msiStatus) %>% count() %>% spread(msiStatus, n, fill = 0) %>% mutate(percentage = MSI /(MSS + MSI)) %>% arrange(percentage)
msiRelativeData = msiRelativeData %>% ungroup() %>% mutate(cancerType = factor(cancerType, msiRelativeData$cancerType))

p1 = ggplot(data = msiRelativeData, aes(y = percentage, x = cancerType)) + 
  geom_bar(stat = "identity", aes(fill = cancerType)) + ylab("% Samples MSI") + xlab("Cancer Type") + ggtitle("") +
  scale_y_continuous(limits = c(0, 0.1), labels = percent, expand=c(0.02, 0.001)) +
  scale_fill_manual(values = cancerTypeColours) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()

p2 = ggplot(data = msiData) +
  stat_ecdf(geom = "point", aes(msiScore), shape = 4) +
  annotate("text", x = 7, y = 0, label = "MSI Cutoff", size = 3, hjust = 0) +
  geom_segment(aes(x = 4, xend = 4, y = 0, yend = 1), linetype = 2) +
  coord_cartesian(xlim = c(0, 10)) +
  scale_y_continuous(labels = percent)+#, expand=c(0.02, 0.001)) +
  xlab("MSI Score") + ylab("% Samples") + ggtitle("") +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  #scale_x_log10() +
  coord_flip()


plot_grid(p1, p2, ncol = 2, labels = "AUTO")


##################################################### TMB ##################################################### 

load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
tmbData = highestPurityCohortSummary %>% select(sampleId, cancerType, TOTAL_SNV) %>% mutate(TMBScore = TOTAL_SNV / 3095, Burden = ifelse(TMBScore >= 10, "Burdened","Unburdened"))
tmbRelativeData = tmbData %>% filter(cancerType != 'Other') %>% group_by(cancerType, Burden) %>% count() %>% spread(Burden, n, fill = 0) %>% mutate(percentage = Burdened /(Burdened + Unburdened)) %>% arrange(percentage)
tmbRelativeData = tmbRelativeData %>% ungroup() %>% mutate(cancerType = factor(cancerType, tmbRelativeData$cancerType))

p3 = ggplot(data = tmbRelativeData, aes(y = percentage, x = cancerType)) + 
  geom_bar(stat = "identity", aes(fill = cancerType)) + ylab("% of samples with TMB > 10/Mb") + xlab("Cancer Type") + ggtitle("") +
  scale_y_continuous(limits = c(0, 0.7), labels = percent, expand=c(0.02, 0.001)) +
  scale_fill_manual(values = cancerTypeColours) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()

p4 = ggplot(data = tmbData) +
  stat_ecdf(geom = "point", aes(TMBScore), shape = 4) +
  annotate("text", x = 23, y = 0, label = "TMB Cutoff", size = 3, hjust = 0) +
  geom_segment(aes(x = 10, xend = 10, y = 0, yend = 1), linetype = 2) +
  scale_y_continuous(labels = percent)+#, expand=c(0.02, 0.001)) +
  xlab("SNVs Per Mb") + ylab("% Samples") + ggtitle("") +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  #scale_x_log10() +
  coord_flip()


plot_grid(p3, p4, ncol = 2, labels = "AUTO")

plot_grid(p1,p2, p3, p4, ncol = 2, nrow = 2, labels = "AUTO")

###################################################### INDELS V SNV


load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")

ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SNV, y = TOTAL_INDEL)) +
  geom_segment(aes(x = 1e2, xend=1e6, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 1e2, y = 15000, label = "MSI Threshold", size = 3, hjust = 0) +
  geom_segment(aes(y = 1e2, yend=1e6, x = 30950, xend=30950), linetype = "dashed") + annotate("text", x = 32000, y = 1.1e2, label = "TMB = 10 per Mb", size = 3, hjust = 0) +
  geom_point(aes(shape = msiStatus, color = cancerType)) + 
  scale_color_manual(values = cancerTypeColours) + 
  scale_x_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  xlab("SNVs") + ylab("Indels")


ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SNV)) +
  stat_ecdf(geom = "point", pad = FALSE, aes(shape = msiStatus)) 




highestPurityCohortSummary[highestPurityCohortSummary$msiStatus == T, ]

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
