detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_bw())

singleBlue = "#6baed6"

#hartwigRed = rgb(211,62,52, maxColorValue = 255)
#hartwigBlue = rgb(49,76,148, maxColorValue = 255)

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

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 


########################################### Figure 2 - WGD
wgdPlotData = highestPurityCohortSummary %>% 
  select(cancerType, WGD) %>%  
  group_by(cancerType, WGD) %>% count() %>%
  group_by(cancerType) %>% mutate(total = sum(n), percentage = n / total) %>%
  ungroup()

wgdPlotDataTotal = wgdPlotData %>% filter(WGD) %>% summarise(percentage = sum(n) / sum(total))
wgdPlotData = wgdPlotData %>%
  mutate(totalPercentage = wgdPlotDataTotal$percentage) %>%
  filter(cancerType != "Other")

wgdPlotLevels = wgdPlotData %>% filter(WGD) %>% arrange(-percentage)
wgdPlotData = mutate(wgdPlotData, cancerType = factor(cancerType, wgdPlotLevels$cancerType))
wgdPlotData[wgdPlotData$cancerType == "CNS", "totalPercentage"] <- NA
wgdPlotData[wgdPlotData$cancerType == "Mesothelioma", "totalPercentage"] <- NA

p1 = ggplot(data = wgdPlotData, aes(x = cancerType, y = percentage)) +
  geom_bar(aes(fill = WGD), stat = "identity") +
  geom_line(aes(x = as.numeric(cancerType), y = totalPercentage), linetype = 2) +
  annotate("text", x = 20, y = wgdPlotDataTotal$percentage, label = "Pan Cancer", size = 3) +
  annotate("text", x = 19, y = wgdPlotDataTotal$percentage, label = sprintf(fmt='(%.1f%%)', 100*wgdPlotDataTotal$percentage), size = 3) +
  #scale_fill_manual(values = c("#f1eef6", "#3182bd")) +
  scale_fill_manual(values = c("#deebf7", "#2171b5")) +
  ggtitle("") + 
  xlab("Cancer Type") + ylab("% Samples")+ 
  scale_y_continuous(labels = percent, expand=c(0.01, 0.01), limits = c(0, 1)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="none") +
  coord_flip()

wgdPDFPlotData = highestPurityCohortSummary %>% select(sampleId, WGD, ploidy)

p2 = ggplot(data=wgdPDFPlotData, aes(x=ploidy, fill = WGD)) +
  geom_histogram(position = "identity", binwidth = 0.1) + 
  scale_fill_manual(values = c(alpha("#bdd7e7", 1), alpha("#2171b5", 0.8))) +
  ggtitle("") +  xlab("Ploidy") + ylab("# Samples") +
  theme(panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(limits = c(0, 7), breaks=c(1:7))

pWGD = plot_grid(p1, p2, labels="AUTO")
save_plot("~/hmf/RPlot/Figure 2 - WGD.png", pWGD, base_width = 7, base_height = 3)

########################################### Figure 4 - Driver Per Sample
driverData = hpcDriversByGene %>%
  mutate(driver = as.character(driver), 
         driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
         driver = ifelse(driver %in% c("Frameshift","Inframe"), 'Indel', driver),
         driver = factor(driver, simplifiedDrivers)
  ) %>%
  group_by(driver, cancerType) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  left_join(hpcCancerTypeCounts %>% select(cancerType, N), by = "cancerType") %>%
  mutate(driversPerSample = driverLikelihood / N) %>%
  arrange(-driversPerSample)

driverDataLevels = driverData %>% 
  group_by(cancerType) %>% 
  summarise(driversPerSample = sum(driversPerSample)) %>% 
  arrange(-driversPerSample)

driverData = driverData %>% 
  mutate(cancerType = factor(cancerType, driverDataLevels$cancerType)) %>%
  group_by(cancerType) %>% 
  mutate(percentage = driversPerSample / sum(driversPerSample))

driverViolinData = hpcDriversByGene %>%
  group_by(cancerType, sampleId) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  ungroup() %>%
  mutate(cancerType = factor(cancerType, driverDataLevels$cancerType))
driverViolinData = merge(driverViolinData,highestPurityCohortSummary %>% select(sampleId,cancerType,cancerSubtype),by=c('sampleId','cancerType'),all=T) 
driverViolinData$driverLikelihood = driverViolinData$driverLikelihood %>% replace_na(0)

p1 = ggplot(driverViolinData, aes(cancerType, driverLikelihood)) + 
  geom_violin(fill = singleBlue, scale = "area") +
  #geom_point(data = driverDataLevels, aes(cancerType, driversPerSample)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=1) +
  xlab("Cancer Type") + ylab("Drivers") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  ggtitle("") + xlab("") + ylab("Drivers per sample") + 
  scale_y_continuous(expand=c(0.01, 0.01)) +
  scale_fill_manual(values = cancerTypeColours) +
  theme(legend.position="none") +
  coord_flip() 

p2 = ggplot(driverData, aes(cancerType, percentage)) +
  geom_bar(stat = "identity", aes(fill = driver), width=0.7) +
  scale_fill_manual(values = simplifiedDriverColours) +
  ggtitle("") + xlab("") + ylab("Drivers by variant type") +  
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), legend.title = element_blank()) +
  scale_y_continuous(labels = percent, expand=c(0.0, 0.0), limits = c(0, 1.01)) +
  coord_flip()

pDriverPerSample = plot_grid(p1,p2, ncol = 2, rel_widths =  c(2.5,3), labels = "AUTO")
save_plot("~/hmf/RPlot/Figure 4 - DriverPerSample.png", pDriverPerSample, base_width = 6, base_height = 4)

########################################### Figure 5 - Hotspots
load(file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")
load(file = "~/hmf/RData/Reference/hpcTertPromoters.Rdata")

tertData = hpcTertPromoters %>% select(sampleId, gene) %>% mutate(variant = "SNV", hotspot = T, nearHotspot = F, driverLikelihood = 1)
hotspotDataInter = hpcDndsOncoDrivers %>% select(sampleId, gene, variant,  hotspot, nearHotspot, driverLikelihood) %>%
  bind_rows(tertData)

hotspotGenes = hotspotData %>% 
  filter(!is.na(hotspot)) %>% 
  group_by(gene) %>% 
  summarise(driverLikelihood = sum(driverLikelihood)) %>% top_n(40, driverLikelihood) %>% arrange(driverLikelihood)

variantFactors = c("INDEL","MNV","SNV")
hotspotFactors = c("Hotspot","NearHotspot","NoHotspot")
hotspotColours = setNames(c("#d94701","#fd8d3c", "#fdbe85"), hotspotFactors)

hotspotData = hotspotDataInter %>% 
  filter(!is.na(hotspot), driverLikelihood > 0, gene %in% hotspotGenes$gene) %>%
  mutate(hotspot = ifelse(hotspot, "Hotspot", "NoHotspot"),
         hotspot = ifelse(nearHotspot, "NearHotspot", hotspot),
         hotspot = factor(hotspot, rev(hotspotFactors)),
         variant = factor(variant, rev(variantFactors))) %>%
  group_by(gene,  variant, hotspot) %>%
  summarise(driverLikelihood = sum(driverLikelihood)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, hotspotGenes$gene)) %>% 
  group_by(gene, hotspot) %>% spread(variant, driverLikelihood, fill = 0)  %>% gather(variant, driverLikelihood, 3:5)

p_hotspot1 = ggplot(data = hotspotData %>% filter(variant == 'SNV'), aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = hotspot), stat = "identity") +
  scale_fill_manual(values = hotspotColours) +
  ggtitle("SNV") + xlab("") + ylab("")+ 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom",  strip.background = element_blank(), legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_flip()

legend = g_legend(p_hotspot1)
p_hotspot1 = p_hotspot1 + theme(legend.position="none")

p_hotspot2 = ggplot(data = hotspotData %>% filter(variant == 'MNV'), aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = hotspot), stat = "identity") +
  scale_fill_manual(values = hotspotColours) +
  ggtitle("MNV") + xlab("") + ylab("")+ 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),  legend.position="none",  strip.background = element_blank(), legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_flip()

p_hotspot3 = ggplot(data = hotspotData  %>% filter(variant == 'INDEL'), aes(x = gene, y = driverLikelihood)) +
  geom_bar(aes(fill = hotspot), stat = "identity") +
  scale_fill_manual(values = hotspotColours) +
  ggtitle("Indel") + xlab("") + ylab("")+ 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), legend.position="none",  strip.background = element_blank(), legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_flip()


pHotspots = plot_grid(p_hotspot1, p_hotspot2, p_hotspot3, rel_widths = c(2,1,1), ncol = 3, labels=c("A"))
pHotspots2 = plot_grid(pHotspots, legend, ncol = 1, rel_heights = c(6,1)) + draw_label("Drivers", x = 0.52, y = 0.15, size = 11)
pHotspots2 
save_plot("~/hmf/RPlot/Figure 5 - Hotspots.png", pHotspots2, base_width = 5, base_height = 6)


########################################### Figure 6 - Clonality
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")

highestPurityCohortSummary$purityBucket = 
  cut(
    highestPurityCohortSummary$purity, 
    breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    labels = c("0-10%", "10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","          90-100%"))

clonalityFactor = c('SUBCLONAL','CLONAL')


clonalityLoad2 = highestPurityCohortSummary %>%
  mutate(total = TOTAL_INDEL + TOTAL_SNV + TOTAL_MNV, 
         subclonal = SUBCLONAL_INDEL + SUBCLONAL_SNV + SUBCLONAL_MNV) %>%
  group_by(purityBucket, sampleId) %>%
  summarise(total = sum(total), SUBCLONAL = sum(subclonal), CLONAL = sum(total) - SUBCLONAL) %>%
  gather(clonality, value, CLONAL, SUBCLONAL) %>%
  mutate(percentage = value / total) %>%
  filter(clonality == 'SUBCLONAL') %>%
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
  geom_bar(fill = singleBlue, stat = "identity") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  ggtitle("") + xlab("Tumor Purity") + ylab("Count of samples") + 
  scale_y_continuous(expand=c(0.01, 0.01)) +
  theme(legend.position="none") +
  coord_flip() 

p2 = ggplot(clonalityLoad2, aes(purityBucket, percentage)) + 
  geom_violin(fill = singleBlue, scale = "width") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=1.5) +
  ggtitle("") + xlab("") + ylab("% of variants subclonal") + 
  scale_y_continuous(limits = c(0, 1), labels = percent, expand=c(0.02, 0.01)) +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  coord_flip() 

p3 = ggplot(data = clonalityDrivers2, aes(x = purityBucket, y = subclonalPercentage, width = 0.7)) +
  geom_bar(fill = singleBlue, stat = "identity") + 
  xlab("") + ylab("% of driver variants subclonal") + ggtitle("") + 
  scale_y_continuous(labels = percent, expand=c(0.01, 0.01), limits = c(0, 0.06)) +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  coord_flip()

pClonality = plot_grid(p1,p2, p3, ncol = 3, labels="AUTO")
pClonality
save_plot("~/hmf/RPlot/Figure 6 - Subclonal.png", pClonality, base_width = 10, base_height = 2.5)


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
  scale_fill_manual(values = c("#fd8d3c", singleBlue)) +
  ggtitle("") + 
  xlab("Driver") + ylab("Average drivers per sample") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  scale_y_continuous(expand=c(0.01, 0.01))

plot_grid(p1, labels="AUTO")

########################################### Extended Figure 1c - coverage
load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
coverageData = highestPurityCohortSummary %>% 
  select(sampleId, tumorMeanCoverage, refMeanCoverage, cancerType) %>%
  mutate(medianTC = median(tumorMeanCoverage, na.rm = T)) %>%
  arrange(cancerType, -tumorMeanCoverage)

pCoverage = ggplot(data=coverageData)+
  stat_ecdf(aes(tumorMeanCoverage,color='Tumor (LHS)'),geom = "step", pad = FALSE) + 
  stat_ecdf(aes(refMeanCoverage*2,color='Reference (RHS)'),geom = "step", pad = FALSE) +
  scale_x_continuous(sec.axis = sec_axis(~./2, name = "Ref Mean Coverage", breaks = c(25, 50, 75, 100))) +
  coord_flip() + ggtitle("") +
  labs(x = "Tumor Mean Coverage")+
  theme(axis.title.x =  element_blank(), axis.ticks = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(labels = percent)

pCoverage = plot_grid(pCoverage, labels = "C")
pCoverage
save_plot("~/hmf/RPlot/Extended Figure 1 - Coverage.png", pCoverage, base_width = 5, base_height = 5)

########################################### Extended Figure 2 - MSI / TMB
load("~/hmf/RData/processed/hpcDriversByGene.RData")
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')

msiData = highestPurityCohortSummary %>% select(sampleId, cancerType, msiScore, msiStatus)
msiRelativeData = msiData %>% group_by(cancerType, msiStatus) %>% count() %>% group_by(cancerType) %>% mutate(percentage = n / sum(n)) 
msiRelativeData = msiData %>% filter(cancerType != 'Other') %>% group_by(cancerType, msiStatus) %>% count() %>% spread(msiStatus, n, fill = 0) %>% mutate(percentage = MSI /(MSS + MSI)) %>% arrange(percentage)
msiRelativeData = msiRelativeData %>% ungroup() %>% mutate(cancerType = factor(cancerType, msiRelativeData$cancerType))
msiGlobalData = msiRelativeData %>% summarise(MSI = sum(MSI), MSS = sum(MSS)) %>% mutate(percentage = MSI /(MSS + MSI))
msiRelativeData$globalPercentage = msiGlobalData$percentage
msiRelativeData[msiRelativeData$cancerType == "CNS", "globalPercentage"] <- NA
msiRelativeData[msiRelativeData$cancerType == "Uterus", "globalPercentage"] <- NA

tmbData = highestPurityCohortSummary %>% select(sampleId, cancerType, TOTAL_SNV) %>% mutate(TMBScore = TOTAL_SNV / 3095, Burden = ifelse(TMBScore >= 10, "Burdened","Unburdened"))
tmbRelativeData = tmbData %>% filter(cancerType != 'Other') %>% group_by(cancerType, Burden) %>% count() %>% spread(Burden, n, fill = 0) %>% mutate(percentage = Burdened /(Burdened + Unburdened)) %>% arrange(percentage)
tmbRelativeData = tmbRelativeData %>% ungroup() %>% mutate(cancerType = factor(cancerType, tmbRelativeData$cancerType))
tmbGlobalData = tmbRelativeData %>% summarise(Burdened = sum(Burdened), Unburdened = sum(Unburdened)) %>% mutate(percentage = Burdened /(Burdened + Unburdened))
tmbRelativeData$globalPercentage = tmbGlobalData$percentage
tmbRelativeData[tmbRelativeData$cancerType == "Skin", "globalPercentage"] <- NA
tmbRelativeData[tmbRelativeData$cancerType == "Lung", "globalPercentage"] <- NA

p1 = ggplot(data = msiRelativeData, aes(y = percentage, x = cancerType)) + 
  geom_bar(stat = "identity", fill = singleBlue) + 
  geom_line(aes(x = as.numeric(cancerType), y = globalPercentage), linetype = 2) +
  annotate("text", x = 20, y = msiGlobalData$percentage, label = "Pan Cancer", size = 3) +
  annotate("text", x = 19, y = msiGlobalData$percentage, label = sprintf(fmt='(%.1f%%)', 100*msiGlobalData$percentage), size = 3) +
  ylab("% Samples MSI") + xlab("Cancer Type") + ggtitle("") +
  scale_y_continuous(limits = c(0, 0.1), labels = percent, expand=c(0.02, 0.001)) +
  scale_fill_manual(values = cancerTypeColours) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()


p2 = ggplot(data = tmbRelativeData, aes(y = percentage, x = cancerType)) + 
  geom_bar(stat = "identity",fill = singleBlue) + 
  geom_line(aes(x = as.numeric(cancerType), y = globalPercentage), linetype = 2) +
  annotate("text", x = 20, y = tmbGlobalData$percentage, label = "Pan Cancer", size = 3) +
  annotate("text", x = 19, y = tmbGlobalData$percentage, label = sprintf(fmt='(%.1f%%)', 100*tmbGlobalData$percentage), size = 3) +
  ylab("% Samples TMB High") + xlab("Cancer Type") + ggtitle("") +
  scale_y_continuous(limits = c(0, 0.7), labels = percent, expand=c(0.02, 0.001)) +
  scale_fill_manual(values = cancerTypeColours) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()

p3 = ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SNV, y = TOTAL_INDEL)) +
  geom_segment(aes(x = 1e2, xend=1e6, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 1e2, y = 15000, label = "MSI Threshold", size = 3, hjust = 0) +
  geom_segment(aes(y = 1e2, yend=1e6, x = 30950, xend=30950), linetype = "dashed") + annotate("text", x = 32000, y = 1.1e2, label = "TMB High Threshold", size = 3, hjust = 0) +
  geom_point(aes(color = cancerType)) + 
  scale_color_manual(values = cancerTypeColours) + 
  scale_x_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
  xlab("SNVs") + ylab("Indels") + ggtitle("")

p12 = plot_grid(p1,p2, rel_widths = c(2,2), labels = c("A","B"))
pMsiTMB = plot_grid(p12, p3, ncol = 1, rel_heights = c(1, 2), labels = c("","C"))

########### Part 2
mutationalLoad = highestPurityCohortSummary %>%
  filter(msiStatus == "MSS") %>%
  mutate(
    indelMutLoad =TOTAL_INDEL,
    snvMutLoad = TOTAL_SNV + TOTAL_MNV,
    svMutLoad = BND + DEL + DUP + INV + INS) %>%
  select(sampleId, cancerType, indelMutLoad, snvMutLoad, svMutLoad) 

driverLoad =  hpcDriversByGene %>%   
  ungroup() %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(driver %in% c("Frameshift", "Inframe"), "indelDriverLoad", driver),
    driver = ifelse(grepl("Indel", impact), "indelDriverLoad", driver),
    driver = ifelse(driver %in% c("Missense", "Nonsense", "Splice"), "snvDriverLoad", driver),
    driver = ifelse(grepl("Missense", impact) | grepl("Nonsense", impact) | grepl("Splice", impact), "snvDriverLoad", driver),
    driver = ifelse(driver %in% c("Amp", "Del","FragileDel"), "cnaDriverLoad", driver)
  ) %>%
  filter(driver %in% c('indelDriverLoad', "snvDriverLoad", "cnaDriverLoad")) %>%
  group_by(sampleId, cancerType, driver) %>% summarise(drivers = sum(driverLikelihood, na.rm = T)) %>%
  spread(driver, drivers, fill = 0)

combinedDriverAndMutationalLoad = merge(mutationalLoad, driverLoad, by = c("sampleId","cancerType"), all.x = T)
combinedDriverAndMutationalLoad[is.na(combinedDriverAndMutationalLoad)] <- 0
combinedDriverAndMutationalLoad = combinedDriverAndMutationalLoad %>%
  group_by(cancerType) %>%
  summarise(
    indelMutLoad = mean(indelMutLoad),
    snvMutLoad = mean(snvMutLoad),
    svMutLoad = mean(svMutLoad),
    indelDriverLoad = mean(indelDriverLoad),
    snvDriverLoad = mean(snvDriverLoad),
    cnaDriverLoad = mean(cnaDriverLoad)
  )


p4 = ggplot(combinedDriverAndMutationalLoad, aes(x=svMutLoad, y=cnaDriverLoad, color=cancerType)) + 
  geom_point() +
  coord_cartesian(xlim= c(0,500), ylim=c(0,5))+
  ylab("CNA Drivers") + xlab("SV Mutational Load") + ggtitle("") + 
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)


p5 = ggplot(combinedDriverAndMutationalLoad, aes(x=snvMutLoad, y=snvDriverLoad, color=cancerType)) + 
  geom_point() +
  coord_cartesian(xlim= c(0,82000), ylim=c(0,4))+
  ylab("SNV Drivers") + xlab("SNV Mutational Load") + ggtitle("") + 
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)

p6 = ggplot(combinedDriverAndMutationalLoad, aes(x=indelMutLoad, y=indelDriverLoad, color=cancerType)) + 
  geom_point() +
  coord_cartesian(xlim= c(0,3000), ylim=c(0,1.5))+
  ylab("Indel Drivers") + xlab("Indel Mutational Load") + ggtitle("") +
  theme(legend.position="none", legend.title = element_blank()) +
  scale_color_manual(values=cancerTypeColours, guide=FALSE)


pLoads = plot_grid(p4,p5,p6, nrow = 1, labels = c("D","E","F"))

pComplete = plot_grid(pMsiTMB, pLoads, nrow = 2, rel_heights = c(3, 1))
pComplete

save_plot("~/hmf/RPlot/Extended Figure 2 - Loads.png", pComplete, base_width = 10, base_height = 14)

########################################### Extended Figure 4 - SMG tile chart

absSignificance = 0.01
load(file = "~/hmf/RData/processed/genePanel.RData")
load(file = "~/hmf/RData/processed/HmfRefCDSCv.RData")

martincorenaGenes = read.csv(file = "~/hmf/resources/Martincorena221.csv", header = T, stringsAsFactors = F) %>% select(gene_name = gene) %>% mutate(martincorena = T)
cosmicGenes = genePanel %>% filter(cosmicCurated | cosmicOncogene | cosmicTsg) %>%
  select(gene_name) %>%
  mutate(cosmic = T)
hmfGenes = HmfRefCDSCv %>% filter(qglobal_cv < absSignificance, !gene_name %in% c("POM121L12","TRIM49B","LPCAT2"), cancerType != 'BOther') %>%
  select(cancerType, gene_name, qglobal_cv) %>%
  mutate(
    gene_name = as.character(gene_name), 
    cancerType = as.character(cancerType),
    cancerType = ifelse(cancerType == 'All', "Pan Cancer", cancerType)) %>%
  left_join(martincorenaGenes, by = c("gene_name")) %>%
  left_join(cosmicGenes, by = "gene_name") %>%
  arrange(gene_name)
hmfGenes[is.na(hmfGenes)] <- F

hmfGenesLabelsAll = hmfGenes %>% filter(cancerType == 'Pan Cancer') %>% select(gene_name, qglobal_cv, martincorena, cosmic)
hmfGenesLabelsOther = hmfGenes %>% group_by(gene_name, martincorena, cosmic) %>% summarise(qglobal_cv = 1)
hmfGenesLabels = bind_rows(hmfGenesLabelsAll, hmfGenesLabelsOther) %>% 
  group_by(gene_name) %>%
  summarise(qglobal_cv = min(qglobal_cv), martincorena = any(martincorena), cosmic = any(cosmic)) %>% 
  arrange(-qglobal_cv, desc(gene_name))

hmfGenesLabels$status = "#cb181d"
hmfGenesLabels$face = "bold"
hmfGenesLabels$status = ifelse(hmfGenesLabels$martincorena, "#fd8d3c", hmfGenesLabels$status)
hmfGenesLabels$status = ifelse(hmfGenesLabels$cosmic, "black", hmfGenesLabels$status)
hmfGenesLabels$face = ifelse(hmfGenesLabels$cosmic, "plain", hmfGenesLabels$face)
hmfGenes = hmfGenes %>% mutate(gene_name = factor(gene_name, hmfGenesLabels$gene_name))

cancerTypeFactors = c('Pan Cancer', hmfGenes %>% distinct(cancerType) %>% arrange(cancerType) %>% filter(cancerType != 'Pan Cancer') %>% pull(cancerType))
hmfGenes = hmfGenes %>% mutate(cancerType = factor(cancerType, cancerTypeFactors))

pTile = ggplot(hmfGenes, aes(x = cancerType, y = gene_name))+
  geom_tile(aes(fill = qglobal_cv+ 1e-13)) + 
  scale_fill_gradient(name = "Significance", trans = "log10", low = "#2171b5", high = "#bdd7e7", guide = "colourbar",limits=c(1e-14, 1e-2), breaks=c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2), labels=c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) +
  xlab("") + ylab("") +
  theme(axis.text.y = element_text(size = 8, colour = hmfGenesLabels$status, face = hmfGenesLabels$face), axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=9)) +
  #theme(panel.border = element_blank(),axis.ticks = element_blank()) +
  #theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 20, direction = "vertical", title.position="top", title.hjust = 0, title.vjust = 0.5, nbin = 50))
pTile  

save_plot("~/hmf/RPlot/Extended Figure 4 - Tile.png", pTile, base_width = 6, base_height = 10)


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
