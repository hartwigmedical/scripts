detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_grey())

#################### SETUP #################### 
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary[is.na(highestPurityCohortSummary)] <- 0
hpcCancerTypeCounts = highestPurityCohortSummary %>% 
  group_by(cancerType) %>% 
  summarise(
    N = n(), 
    medianMutationalLoad = median(CLONAL_SNP + SUBCLONAL_SNP + CLONAL_MNP + SUBCLONAL_MNP + CLONAL_INDEL + SUBCLONAL_INDEL )) %>% 
  arrange(-medianMutationalLoad)
cancerTypeFactors =  factor(hpcCancerTypeCounts$cancerType, levels = hpcCancerTypeCounts$cancerType)


save(hpcCancerTypeCounts, file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
cancerTypes = sort(unique(highestPurityCohortSummary$cancerType))
cancerTypeColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
cancerTypeColours = setNames(cancerTypeColours[1:length(cancerTypes)], cancerTypes)
save(cancerTypeColours, file = "~/hmf/RData/Reference/cancerTypeColours.RData")

singleSubstitutionColours = c("#14B0EF","#060809","#E00714","#BFBEBF","#90CA4B","#E9BBB8")
singleSubstitutionColours = setNames(singleSubstitutionColours, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

doubleSubstitutions = c("CC>AA","CC>TT","CC>Other","TC","TT","CT","TG","AC","GC","Other")
doubleSubstitutionColours = c("#374884","#E00714","#7AB154", "#060809","#BF4C86","#F1D04F","#2B4C2E","#7575A6","#864F46","#82BC9F")
doubleSubstitutionColours = setNames(doubleSubstitutionColours, doubleSubstitutions)




#################### PREPARE DATA #################### 
agePlotData = highestPurityCohortSummary %>% 
  select(sampleId, ageAtBiopsy, cancerType) %>% 
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  arrange(cancerType, -ageAtBiopsy)

cancerTypeData = highestPurityCohortSummary %>% 
  group_by(cancerType) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  mutate(percentage = round(n / n(), 1)) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors))

hmfMutationalLoad = highestPurityCohortSummary %>% 
  select(sampleId, cancerType, ends_with("INDEL"), ends_with("SNP"), ends_with("MNP"), BND, DEL, INS, INV, DUP) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  mutate(
    INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL,
    MNP = INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
    SNP = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
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
  filter(!is.na(cancerType)) %>% select(PCAWG_SNP = all.SNVs, PCAWG_INDEL = all.Indels, PCAWG_SV = SV.events, age,PCAWG_MNP = all.MNVs,  cancerType) %>%
  mutate(source = "PCAWG")

combinedMutationalLoad =  hmfMutationalLoad %>% select(sampleId, cancerType, INDEL, SNP, MNP, SV) %>%
  mutate(source = "HMF") %>%
  bind_rows(pcawgMutationalLoad) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors))

load(file = "~/hmf/RData/Reference/allSNPSummary.RData")
hpcSNP = allSNPSummary %>% 
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  filter(nchar(alt) == 1) %>% 
  mutate(type = standard_mutation(paste(ref, alt, sep = '>'))) %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad)
hpcSNP$sampleId = factor(hpcSNP$sampleId, levels = unique(hpcSNP$sampleId))


load(file = "~/hmf/RData/Reference/allMNPSummary.RData")
hpcMNP = allMNPSummary %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  filter(nchar(ref) == 2, nchar(alt) == 2) %>%
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
  arrange(sampleMutationalLoad)
hpcMNP$sampleId = factor(hpcMNP$sampleId, levels = unique(hpcMNP$sampleId))


#################### Cancer Type Summary FACETED #################### 
p1 = ggplot(data=cancerTypeData, aes(x = NA, y = n)) +
  geom_bar(aes(fill = cancerType), stat = "identity") +
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) + 
  geom_text(aes(label=paste0("(", percentage, "%)")), vjust=-0.5, size = 3) +
  geom_text(aes(label=n), vjust=-2, size = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 5)) +  
  ylab("Samples") + 
  coord_cartesian(ylim = c(0, 600)) + facet_grid(~cancerType)

p2 = ggplot(agePlotData, aes(NA, ageAtBiopsy)) + 
  geom_violin(aes(fill=cancerType), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") + 
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  xlab("Cancer Type") + ylab("Age at Biopsy") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_blank()) + 
  facet_grid(~cancerType)

p3 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(SNP,color='SNP'), geom = "step", pad = FALSE)+
  stat_ecdf(aes(INDEL,color='INDEL'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(SV,color='SV'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(MNP,color='MNP'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(PCAWG_SNP,color='SNP'),geom = "step", linetype = "dashed", pad = FALSE)+
  stat_ecdf(aes(PCAWG_INDEL,color='INDEL'),geom = "step", linetype = "dashed", pad = FALSE)+
  stat_ecdf(aes(PCAWG_SV,color='SV'),geom = "step", linetype = "dashed", pad = FALSE)+
  stat_ecdf(aes(PCAWG_MNP,color='MNP'),geom = "step", linetype = "dashed", pad = FALSE)+
  scale_x_log10() + facet_grid(~cancerType) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text = element_blank(), legend.position="top", legend.title = element_blank()) + 
  xlab("Mutational Load") + 
  coord_flip()

p4 = ggplot(data=hpcSNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=singleSubstitutionColours) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text.x =element_blank(), 
        legend.position="bottom", legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1))


p5 = ggplot(data=hpcMNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=doubleSubstitutionColours) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text.x =element_blank(), 
        legend.position="bottom", legend.title = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  facet_grid(~cancerType, scales = "free_x") + 
  guides(fill = guide_legend(nrow = 1))

plot_grid(p1, p2, p3, p4, p5, ncol=1, align="v", rel_heights = c(1, 1, 3, 1, 1), labels = c("A", "B", "C", "D", "E"))

