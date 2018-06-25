detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_grey())
theme_set(theme_bw())

#################### SETUP #################### 
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary[is.na(highestPurityCohortSummary)] <- 0
hpcCancerTypeCounts = highestPurityCohortSummary %>% 
  group_by(cancerType) %>% 
  summarise(
    N = n(), 
    medianMutationalLoad = median(CLONAL_SNP + SUBCLONAL_SNP + CLONAL_MNP + SUBCLONAL_MNP + CLONAL_INDEL + SUBCLONAL_INDEL )) %>% 
  arrange(medianMutationalLoad)
cancerTypeFactors =  factor(hpcCancerTypeCounts$cancerType, levels = hpcCancerTypeCounts$cancerType)


save(hpcCancerTypeCounts, file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
cancerTypes = sort(unique(highestPurityCohortSummary$cancerType))
cancerTypeColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
cancerTypeColours = setNames(cancerTypeColours[1:length(cancerTypes)], cancerTypes)
save(cancerTypeColours, file = "~/hmf/RData/Reference/cancerTypeColours.RData")

somaticColours = c("#a6611a","#dfc27d","#80cdc1","#018571")
somaticColours = setNames(somaticColours, c("HMF SNV","PCAWG SNV", "PCAWG MNV", "HMF MNV"))

indelSVColours = c("#d01c8b","#f1b6da","#b8e186","#4dac26")
indelSVColours = setNames(indelSVColours, c("HMF INDEL","PCAWG INDEL", "PCAWG SV", "HMF SV"))

singleSubstitutionColours = c("#14B0EF","#060809","#E00714","#BFBEBF","#90CA4B","#E9BBB8")
singleSubstitutionColours = setNames(singleSubstitutionColours, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

doubleSubstitutions = c("CC>TT","CC>AA","CC>NN","TC>NN","TT>NN","AC>NN","GC>NN","TG>NN","CT>NN","TA>NN","CG>NN","AT>NN","Other")
doubleSubstitutionColours = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#060809")
doubleSubstitutionColours = setNames(doubleSubstitutionColours, doubleSubstitutions)

indelColours = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
indelColours = setNames(indelColours, c("INS-repeat","INS-other","DEL-repeat", "DEL-other", "DEL-MH"))

svTypes = c("DUP","DEL","BND","INS","INV")
svColours = c("#33a02c","#e31a1c","#1f78b4","#ffff33","#060809","#984ea3")
svColours = setNames(svColours, svTypes)


#################### PREPARE DATA #################### 
agePlotData = highestPurityCohortSummary %>% 
  filter(cancerType != "Other") %>%
  select(sampleId, ageAtBiopsy, cancerType) %>% 
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  arrange(cancerType, -ageAtBiopsy)

cancerTypeData = highestPurityCohortSummary %>% 
  group_by(cancerType) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  mutate(percentage = round(n / n(), 1)) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  filter(cancerType != "Other")

hmfMutationalLoad = highestPurityCohortSummary %>% 
  select(sampleId, cancerType, ends_with("INDEL"), ends_with("SNP"), ends_with("MNP"), BND, DEL, INS, INV, DUP) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  mutate(
    INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL,
    MNV = INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
    SNV = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
    SV = BND + DEL + INS + INV + DUP) 

pcawgRaw = read.csv("~/hmf/resources/PCAWG_counts.txt", sep = '\t', stringsAsFactors = F)
pcawg_histology_tier2 = sort(unique(pcawgRaw$histology_tier2))

pcawgCancerTypeMapping = data.frame(histology_tier2 = pcawg_histology_tier2, cancerType = pcawg_histology_tier2, stringsAsFactors = F)
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bladder", "cancerType"] = "Urinary tract"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bone/SoftTissue", "cancerType"] = "Bone/Soft tissue"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Cervix", "cancerType"] = NA
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Head/Neck", "cancerType"] = "Head and neck"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Myeloid", "cancerType"] = "Other"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Lymphoid", "cancerType"] = "Other"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Thyroid", "cancerType"] = "Other"
pcawgCancerTypeMapping = pcawgCancerTypeMapping[!is.na(pcawgCancerTypeMapping$cancerType), ]

pcawgMutationalLoad = pcawgRaw %>% left_join(pcawgCancerTypeMapping, by = "histology_tier2") %>%
  filter(!is.na(cancerType)) %>% select(PCAWG_SNV = all.SNVs, PCAWG_INDEL = all.Indels, PCAWG_SV = SV.events, age,PCAWG_MNV = all.MNVs,  cancerType) %>%
  mutate(source = "PCAWG")

combinedMutationalLoad =  hmfMutationalLoad %>% select(sampleId, cancerType, INDEL, SNV, MNV, SV) %>%
  mutate(source = "HMF") %>%
  bind_rows(pcawgMutationalLoad) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  filter(cancerType != "Other")

combinedMutationalLoad = combinedMutationalLoad %>% 
  group_by(cancerType) %>% 
  mutate(
    medianSNV = median(SNV, na.rm = T), 
    medianMNV = median(MNV, na.rm = T), 
    medianPCAWG_SNV = median(PCAWG_SNV, na.rm = T),
    medianPCAWG_MNV = median(PCAWG_MNV, na.rm = T),
    medianINDEL = median(INDEL, na.rm = T), 
    medianSV = median(SV, na.rm = T), 
    medianPCAWG_INDEL = median(PCAWG_INDEL, na.rm = T),
    medianPCAWG_SV = median(PCAWG_SV, na.rm = T)
  ) %>% ungroup()

load(file = "~/hmf/RData/Reference/allSNPSummary.RData")
hpcSNP = allSNPSummary %>% 
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  filter(nchar(ref) == 1, nchar(alt) == 1) %>% 
  mutate(non_standard_type = paste(ref, alt, sep = '>')) %>%
  mutate(type = standard_mutation(non_standard_type)) %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad) %>%
  filter(cancerType != "Other")
hpcSNP$sampleId = factor(hpcSNP$sampleId, levels = unique(hpcSNP$sampleId))

load(file = "~/hmf/RData/Reference/allMNPSummary.RData")
hpcMNP = allMNPSummary %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  mutate(type = standard_double_mutation(paste(ref, alt, sep = '>'))) %>%
  mutate(type = factor(type, levels = doubleSubstitutions)) %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad) %>%
  filter(cancerType != "Other")
hpcMNP$sampleId = factor(hpcMNP$sampleId, levels = unique(hpcMNP$sampleId))

load(file = "~/hmf/RData/Reference/allIndelSummary.RData")
hpcINDEL = allIndelSummary %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  group_by(sampleId, type = category) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%       
  ungroup() %>%
  arrange(sampleMutationalLoad) %>%
  filter(cancerType != "Other")
hpcINDEL$sampleId = factor(hpcINDEL$sampleId, levels = unique(hpcINDEL$sampleId))      
                       
hpcSV = highestPurityCohortSummary %>% select(sampleId, cancerType, BND, DEL, DUP, INS, INV) %>%
  gather(type, n, BND, DEL, DUP, INS, INV) %>%
  group_by(sampleId, cancerType, type) %>%
  summarise(n = sum(n)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%       
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad) %>%
  filter(cancerType != "Other")
hpcSV$sampleId = factor(hpcSV$sampleId, levels = unique(hpcSV$sampleId)) 

#################### Cancer Type Summary FACETED #################### 
p1 = ggplot(data=cancerTypeData, aes(x = NA, y = n)) +
  geom_bar(aes(fill = cancerType), stat = "identity") +
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) + 
  geom_text(aes(label=paste0("(", percentage, "%)")), vjust=-0.5, size = 2) +
  #geom_text(aes(label=n), vjust=-2, size = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 5)) +  
  ylab("Samples") + 
  coord_cartesian(ylim = c(0, 600)) + facet_grid(~cancerType)

p2 = ggplot(agePlotData, aes(NA, ageAtBiopsy)) + 
  geom_violin(aes(fill=cancerType), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") + 
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  scale_colour_manual(values=cancerTypeColours, guide=FALSE) +
  ylab("Age") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_blank()) + 
  facet_grid(~cancerType)

p3 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(SNV,color='HMF SNV'), geom = "step", pad = FALSE) + geom_segment(aes(x = medianSNV, xend = medianSNV, y = 0.25, yend = 0.75, color='HMF SNV'), show.legend = F) + 
  stat_ecdf(aes(MNV,color='HMF MNV'),geom = "step", pad = FALSE) + geom_segment(aes(x = medianMNV, xend = medianMNV, y = 0.25, yend = 0.75, color='HMF MNV'), show.legend = F) + 
  stat_ecdf(aes(PCAWG_SNV,color='PCAWG SNV'), geom = "step", linetype = "dashed", pad = FALSE) + geom_segment(aes(x = medianPCAWG_SNV, xend = medianPCAWG_SNV, y = 0.25, yend = 0.75, color='PCAWG SNV'), show.legend = F) + 
  stat_ecdf(aes(PCAWG_MNV,color='PCAWG MNV'), geom = "step", linetype = "dashed", pad = FALSE) + geom_segment(aes(x = medianPCAWG_MNV, xend = medianPCAWG_MNV, y = 0.25, yend = 0.75, color='PCAWG MNV'), show.legend = F) + 
  scale_x_log10() + facet_grid(~cancerType) +
  scale_colour_manual(values=somaticColours) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), legend.position="top", legend.title = element_blank()) + 
  xlab("Somatic Variants") +
  coord_flip()


p4 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(INDEL,color='HMF INDEL'),geom = "step", pad = FALSE) + geom_segment(aes(x = medianINDEL, xend = medianINDEL, y = 0.25, yend = 0.75, color='HMF INDEL'), show.legend = F) + 
  stat_ecdf(aes(SV,color='HMF SV'),geom = "step", pad = FALSE) + geom_segment(aes(x = medianSV, xend = medianSV, y = 0.25, yend = 0.75, color='HMF SV'), show.legend = F) +
  stat_ecdf(aes(PCAWG_INDEL,color='PCAWG INDEL'),geom = "step", linetype = "dashed", pad = FALSE) + geom_segment(aes(x = medianPCAWG_INDEL, xend = medianPCAWG_INDEL, y = 0.25, yend = 0.75, color='PCAWG INDEL'), show.legend = F) +
  stat_ecdf(aes(PCAWG_SV,color='PCAWG SV'),geom = "step", linetype = "dashed", pad = FALSE) + geom_segment(aes(x = medianPCAWG_SV, xend = medianPCAWG_SV, y = 0.25, yend = 0.75, color='PCAWG SV'), show.legend = F) +
  scale_x_log10() + facet_grid(~cancerType) +
  scale_colour_manual(values=indelSVColours) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom", legend.title = element_blank()) + 
  xlab("Somatic Variants") +
  xlab("INDELs & SVs") + 
  coord_flip()

p5 = ggplot(data=hpcSNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=singleSubstitutionColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(), 
    strip.text.x = element_text(size = 5), 
    legend.position="bottom", legend.title = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) 

p6 = ggplot(data=hpcMNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=doubleSubstitutionColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x =element_blank(), 
    legend.position="bottom", legend.title = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) 

p7 = ggplot(data=hpcINDEL, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=indelColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x =element_blank(), 
    legend.position="bottom", legend.title = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) 

p8 = ggplot(data=hpcSV, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=svColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x =element_blank(), 
    legend.position="bottom", legend.title = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) 

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=1, align="v", rel_heights = c(1, 1, 3, 3, 2, 2, 2, 2), labels = c("A", "B", "C", "D", "E", "F", "G", "H"))


####################################
### COVERAGE PLOT @@@@@@@@
coverageData = highestPurityCohortSummary %>% 
  select(sampleId, tumorMeanCoverage, refMeanCoverage, cancerType) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors), medianTC = median(tumorMeanCoverage, na.rm = T)) %>%
  arrange(cancerType, -tumorMeanCoverage)

ggplot(data=coverageData)+
  stat_ecdf(aes(tumorMeanCoverage,color='Tumor (LHS)'),geom = "step", pad = FALSE) + 
  stat_ecdf(aes(refMeanCoverage*2,color='Ref (RHS'),geom = "step", pad = FALSE) +
  scale_x_continuous(sec.axis = sec_axis(~./2, name = "Ref Mean Coverage")) +
  coord_flip() + 
  labs(x = "Tumor Mean Coverage")+
  theme(axis.title.x =  element_blank()) +
  scale_y_continuous(labels = percent)

####################################
### Purity PLOT @@@@@@@@
purityData = highestPurityCohortSummary %>% 
  select(sampleId, purity,cancerType ) %>% 
  arrange(cancerType, -purity)

ggplot(data=purityData)+
  stat_ecdf(aes(purity,color='Purity'),geom = "step", pad = FALSE) + 
  coord_flip() + 
  labs(x = "Purity")+
  theme(axis.title.x =  element_blank()) +
  scale_x_continuous(labels = percent) + 
  scale_y_continuous(labels = percent)

###################################
##### Biopsy Location
#### To do: order properly (largerst to smallest with other at top)
biopsyColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
               "#6b3a9d", "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a",
               "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
               "#dea185", "#a0729d", "#8a392f")
head(highestPurityCohortSummary)

biopsyTypeCount = highestPurityCohortSummary  %>% group_by(biopsyType) %>% count() %>% mutate(`Biopsy Type`=ifelse(n<30,'Other',biopsyType)) %>% group_by(`Biopsy Type`) %>% summarise(n=sum(n))

ggplot(biopsyTypeCount, aes(x="",y=n, fill=`Biopsy Type`)) + 
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=biopsyColours) + 
  theme(axis.title.x =  element_blank(),axis.title.y =  element_blank(),axis.text.x =  element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
head(highestPurityCohortSummary)


###########@@@@@@@@@@@@@@@
####### WGD
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
  scale_fill_manual(values = c("#deebf7", "#3182bd")) +
  ggtitle("Whole Genome Duplication") + 
  xlab("Cancer Type") + ylab("% Samples")+ 
  scale_y_continuous(labels = percent, expand=c(0.01, 0.01), limits = c(0, 1)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="none") +
  coord_flip()

wgdPDFPlotData = highestPurityCohortSummary %>% 
  select(sampleId, WGD, ploidy)
  
p2 = ggplot(data=wgdPDFPlotData, aes(x=ploidy, fill = WGD)) +
  geom_histogram(position = "identity", binwidth = 0.1) + 
  scale_fill_manual(values = c(alpha("#bdd7e7", 1), alpha("#3182bd", 0.8))) +
  ggtitle("") +  xlab("Ploidy") + ylab("# Samples") +
  theme(panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(limits = c(0, 7), breaks=c(1:7))

plot_grid(p1, p2, labels="AUTO")
#deebf7

###########################
###### SNV vs INDEL
### To do:  get cancerType colours work properly
# INDEL SNV by MSI
ggplot(highestPurityCohortSummary,aes(CLONAL_SNP,CLONAL_INDEL,color=msiStatus))+geom_point() + 
  scale_x_log10() + scale_y_log10()

# INDEL SNV by CancerType
ggplot(highestPurityCohortSummary,aes(CLONAL_SNP,CLONAL_INDEL,color=cancerType))+geom_point() + 
  scale_x_log10() + scale_y_log10() + scale_fill_manual(values=cancerTypeColours)

# MNV SNV by Cancer Type
ggplot(highestPurityCohortSummary, aes(CLONAL_SNP,CLONAL_MNP,color=cancerType))+geom_point() + 
  scale_x_log10() + scale_y_log10() + scale_fill_manual(values=cancerTypeColours)

