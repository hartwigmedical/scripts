detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)

#################### COLOURS #################### 
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
hpcCancerTypeCounts = highestPurityCohortSummary %>% group_by(cancerType) %>% summarise(N = n()) %>% arrange(-N)
save(hpcCancerTypeCounts, file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
cancerTypes = sort(unique(highestPurityCohortSummary$cancerType))
cancerTypeColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
cancerTypeColours = setNames(cancerTypeColours[1:length(cancerTypes)], cancerTypes)
save(cancerTypeColours, file = "~/hmf/RData/Reference/cancerTypeColours.RData")

#################### LOAD #################### 
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/Reference/hpcCancerTypeCounts.RData")
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
cancerTypeFactors =  factor(hpcCancerTypeCounts$cancerType, levels = hpcCancerTypeCounts$cancerType)

#################### Counts of cancer type #################### 
cancerTypeData = highestPurityCohortSummary  %>% group_by(cancerType) %>% summarise(n = n()) %>% arrange(-n)
cancerTypeData$cancerType = factor(cancerTypeData$cancerType, levels = cancerTypeFactors)

ggplot(data=cancerTypeData, aes(x = cancerType, y = n)) +
  geom_bar(aes(fill = cancerType), stat = "identity") +
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  xlab("Cancer Type") + ylab("Number of samples") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")

rm(cancerTypeData)

#################### AGE Violin Plot by Cancer Type #################### 
agePlotData = highestPurityCohortSummary %>% select(sampleId, ageAtBiopsy, cancerType) %>% arrange(cancerType, -ageAtBiopsy)
agePlotData$cancerType = factor(agePlotData$cancerType, levels = cancerTypeFactors)

ggplot(agePlotData, aes(factor(cancerType), ageAtBiopsy)) + 
  geom_violin(aes(fill=cancerType), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") + 
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  xlab("Cancer Type") + ylab("Age at Biopsy") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")

rm(agePlotData)


#################### HMF Mutational Load #################### 
hmfMutationalLoad = highestPurityCohortSummary %>% select(sampleId, cancerType, ends_with("INDEL"), ends_with("SNP"), ends_with("MNP"), BND, DEL, INS, INV, DUP)
hmfMutationalLoad[is.na(hmfMutationalLoad)] <- 0
hmfMutationalLoad = hmfMutationalLoad %>% mutate(
  INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL,
  MNP = INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
  SNP = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
  SV = BND + DEL + INS + INV + DUP)

ggplot(data=hmfMutationalLoad) +
  stat_ecdf(aes(SNP,color='SNP'), geom = "step", pad = FALSE)+
  stat_ecdf(aes(INDEL,color='INDEL'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(SV,color='SV'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(MNP,color='MNP'),geom = "step", pad = FALSE)+
  scale_x_log10() + facet_wrap(~cancerType)

#################### PCAWG Mutational Load #################### 
pcawgRaw = read.csv("~/hmf/resources/PCAWG_counts.txt", sep = '\t', stringsAsFactors = F)
pcawg_histology_tier2 = sort(unique(pcawgRaw$histology_tier2))

pcawgCancerTypeMapping = data.frame(histology_tier2 = pcawg_histology_tier2, cancerType = pcawg_histology_tier2, stringsAsFactors = F)
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bladder", "cancerType"] = "Urinary tract"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bone/SoftTissue", "cancerType"] = "Bone/Soft tissue"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Cervix", "cancerType"] = NA
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Head/Neck", "cancerType"] = "Head and neck"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Myeloid", "cancerType"] = NA
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Thyroid", "cancerType"] = NA
pcawgCancerTypeMapping = pcawgCancerTypeMapping[!is.na(pcawgCancerTypeMapping$cancerType), ]

pcawgMutationalLoad = pcawgRaw %>% left_join(pcawgCancerTypeMapping, by = "histology_tier2") %>%
  filter(!is.na(cancerType)) %>% select(PCAWG_SNP = all.SNVs, PCAWG_INDEL = all.Indels, PCAWG_SV = SV.events, age, cancerType) %>%
  mutate(source = "PCAWG")

#################### Combined Mutational Load #################### 
combinedMutationalLoad =  hmfMutationalLoad %>%
  filter(cancerType %in% pcawgCancerTypeMapping$cancerType) %>% select(sampleId, cancerType, INDEL, SNP, SV) %>%
  mutate(source = "HMF") %>%
  bind_rows(pcawgMutationalLoad)


ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(SNP,color='SNP'), geom = "step", pad = FALSE)+
  stat_ecdf(aes(INDEL,color='INDEL'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(SV,color='SV'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(PCAWG_SNP,color='SNP'),geom = "step", linetype = "dashed", pad = FALSE)+
  stat_ecdf(aes(PCAWG_INDEL,color='INDEL'),geom = "step", linetype = "dashed", pad = FALSE)+
  stat_ecdf(aes(PCAWG_SV,color='SV'),geom = "step", linetype = "dashed", pad = FALSE)+
  scale_x_log10() + facet_wrap(~cancerType)


#################### SNPS #################### 
load(file = "~/hmf/RData/Reference/allSNPSummary.RData")
allSNPSummary = allSNPSummary %>% 
  filter(nchar(alt) == 1) %>% 
  mutate(type = standard_mutation(paste(ref, alt, sep = '>'))) %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId")

allSNPSummarySampleTotals = allSNPSummary %>% ungroup() %>% group_by(sampleId) %>% summarise(sampleTotal = sum(n))
allSNPSummaryCancerTypeTotals = allSNPSummary %>% ungroup() %>% group_by(cancerType) %>% summarise(cancerTypeTotal = sum(n))

hpcSNP = allSNPSummary %>% filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  left_join(allSNPSummarySampleTotals, by = "sampleId") %>%
  left_join(allSNPSummaryCancerTypeTotals, by = "cancerType") %>%
  mutate(sampleRelativeN = n / sampleTotal) %>%
  mutate(cancerTypeRelativeN = n / cancerTypeTotal)

ggplot(data=hpcSNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity") +
#  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  xlab("Samples") + ylab("Base substitution") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position="bottom")



private static final String PINK = "233,187,184";
private static final String BLUE = "20,176,239";
private static final String BLACK = "6,8,9";
private static final String RED = "224,7,20";
private static final String GREY = "191,190,191";
private static final String GREEN = "144,202,75";


private static String color(PurityAdjustedSomaticVariant variant) {
  
  if (signature("C", "A", variant)) return BLUE;
  if (signature("C", "G", variant)) return BLACK;
  if (signature("C", "T", variant)) return RED;
  if (signature("T", "A", variant)) return GREY;
  if (signature("T", "C", variant)) return GREEN;
  if (signature("T", "G", variant)) return PINK;
  
  return "purple";
}


########################### RUBBISH ###########################



#load(file = "~/hmf/RData/reference/highestPuritySomaticSummary.RData")
#hmfSomatics = highestPuritySomaticSummary
#hmfSomatics[is.na(hmfSomatics)] = 0
#hmfSomatics = hmfSomatics %>% mutate(INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL, SNP = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP)


cohortSummary[is.na(cohortSummary)] <- 0
cohortSummary$age <- ifelse(is.na(cohortSummary$biopsyDate), NA, as.numeric(substr(cohortSummary$biopsyDate,0,4)) - cohortSummary$birthYear)
cohortSummary$primaryTumorLocation <- ifelse(cohortSummary$primaryTumorLocation == 0, NA, cohortSummary$primaryTumorLocation)

hmfMutationalLoad = cohortSummary %>% mutate(
  INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL,
  SNP = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP + INCONSISTENT_MNP + SUBCLONAL_MNP + CLONAL_MNP,
  SV = BND + DEL + INS + INV + DUP)

clonalitySummary = cohortSummary %>%
  select(sampleId, primaryTumorLocation, INCONSISTENT_INDEL, SUBCLONAL_INDEL, CLONAL_INDEL, INCONSISTENT_SNP, SUBCLONAL_SNP, CLONAL_SNP) %>%
  gather(type, value, INCONSISTENT_INDEL, SUBCLONAL_INDEL, CLONAL_INDEL, INCONSISTENT_SNP, SUBCLONAL_SNP, CLONAL_SNP) %>%
  separate(type, c("clonality", "type")) %>% group_by(primaryTumorLocation, clonality, type) %>% summarise(value = sum(value)) %>%
  left_join(cohortByPrimaryTumorLocation, by = "primaryTumorLocation") %>% mutate(value = value / N)




jonSNPLevels = jon %>% filter(type == 'SNP') %>% group_by(primaryTumorLocation) %>% summarise(n = sum(value)) %>% arrange(-n)
jon$primaryTumorLocation = factor(jon$primaryTumorLocation, levels=jonSNPLevels$primaryTumorLocation)

ggplot(data=jon %>% filter(type == 'SNP'), aes(x = primaryTumorLocation, y = value)) +
  geom_bar(aes(fill = clonality), stat = "identity") +
  ggtitle("SNPs per cancer type") + xlab("Primary Tumor Location") + ylab("SNPs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), legend.position="bottom")

ggplot(data=jon %>% filter(type == 'SNP'), aes(x = primaryTumorLocation, y = value)) +
  geom_bar(aes(color = clonality), stat = "identity") +
  ggtitle(paste0("TSG - " ,"Missense")) + xlab("Codons") + ylab("")


ggplot(data=hmfMutationalLoad) +
  stat_ecdf(aes(SNP,color='SNP'), geom = "step", pad = FALSE)+
  stat_ecdf(aes(INDEL,color='INDEL'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(SV,color='SV'),geom = "step", pad = FALSE)+
  scale_x_log10() + facet_wrap(~primaryTumorLocation)









ggplot(data=cohortSummary ) +
  stat_ecdf(aes(SNP,colour='SNP'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(INDEL,colour='INDEL'),geom = "step", pad = FALSE)+
  stat_ecdf(aes(SV,colour='SV'),geom = "step", pad = FALSE)+
  scale_x_log10() + facet_wrap(~primaryTumorLocation)




pcawgSomatics = pcawgRaw[!is.na(pcawgRaw$cancerType), c("all.SNVs", "all.Indels", "SV.events", "age", "cancerType")]
pcawgSomatics$source <- "PCAWG"
colnames(pcawgSomatics) <- c("snvs", "indels", "svs", "age", "cancerType", "source")

#### VISUALISATION
load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
cohortByPrimaryTumorLocation = highestPurityCohort %>% group_by(primaryTumorLocation) %>% summarise(N = n())
save(cohortByPrimaryTumorLocation, file = '~/hmf/RData/reference/cohortByPrimaryTumorLocation.RData')
primaryTumorLocations = unique(highestPurityCohort$primaryTumorLocation)
primaryTumorLocations= primaryTumorLocations[!is.na(primaryTumorLocations)]

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
primaryTumorLocationColours = setNames(cosmicSignatureColours[1:length(primaryTumorLocations)], primaryTumorLocations)
save(primaryTumorLocationColours, file = "~/hmf/RData/reference/primaryTumorLocationColours.RData")



load(file = "~/hmf/RData/reference/primaryTumorLocationColours.RData")


cohortSummary$age=as.numeric(str_sub(cohortSummary$samplingDate,1,4))-cohortSummary$birthYear   # Need to make this more accurate
cohortSummary = cohortSummary %>% mutate(SV=BND+DUP+DEL+INV)
PCAWGCounts = read.table('~/Dropbox/HMF Australia team folder/PCAWG Data/PCAWG_counts.txt',sep='\t',header = TRUE)
View(PCAWGCounts)

ggplot(data=cohortSummary,aes(SNP,INDEL))+geom_point(aes(colour=primaryTumorLocation))+scale_x_log10()+scale_y_log10() +
  scale_colour_manual( values= primaryTumorLocationColours) +
  theme(legend.position="bottom")

ggplot(data=hmfMutationalLoad,aes(age,SNP))+geom_point(aes(colour=primaryTumorLocation))+scale_y_log10() +
  scale_colour_manual( values= primaryTumorLocationColours) +
  theme(legend.position="bottom") + facet_wrap(~primaryTumorLocation)






















############# Load PURPLE DATA
load("~/hmf/purple.RData")
cohort$age <- ifelse(is.na(cohort$biopsyDate), NA, as.numeric(substr(cohort$biopsyDate,0,4)) - cohort$birthYear)
purpleRaw = cohort
rm(allSamples)
rm(cohort)
rm(backupSamples)

# TRANSFORM PURPLE into usable
purpleVariants = purpleRaw[!is.na(purpleRaw$cancerType), c("snv", "indels", "svTotal", "age", "cancerType")]
purpleVariants$source <- "PURPLE"
colnames(purpleVariants) <- c("snvs", "indels", "svs", "age", "cancerType", "source")

# PURPLE CancerTypes
unique(purpleVariants$cancerType)

#############Load PCAWG Data
pcawgRaw = read.csv("/Users/jon/hmf/pcawg/PCAWG_counts.txt", sep = '\t')
pcawgRaw$cancerType <- NA

# Mapping helpers
#unique(purpleVariants$cancerType)
#unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('organ_system')])
#unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('tumour_type_abbr')])

### Map PCAWG cancerTypes to PURPLE cancerType
pcawgRaw$cancerType <- ifelse(pcawgRaw$tumour_type_abbr == "Melanoma", "Melanoma", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$tumour_type_abbr == "Sarcoma_soft_tissue", "Sarcoma", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "BREAST", "Breast", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "LIVER", "Liver", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "KIDNEY", "Kidney", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "OVARY", "Ovary", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "PROSTATE GLAND", "Prostate", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "STOMACH", "Stomach", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "PANCREAS", "Pancreas", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "UTERUS, NOS", "Uterus", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "BRAIN, & CRANIAL NERVES, & SPINAL CORD, (EXCL. VENTRICLE, CEREBELLUM)", "Brain", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "LUNG & BRONCHUS", "Lung", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "GALLBLADDER & EXTRAHEPATIC BILE DUCTS", "Gall bladder", pcawgRaw$cancerType)

# Mapped successfully
#dim(pcawgRaw[!is.na(pcawgRaw$cancerType),])

# Missing
#pcawgRaw[is.na(pcawgRaw$cancerType),]

# TRANSFORM PCAWG into usable
pcawgSomatics = pcawgRaw[!is.na(pcawgRaw$cancerType), c("all.SNVs", "all.Indels", "SV.events", "age", "cancerType")]
pcawgSomatics$source <- "PCAWG"
colnames(pcawgSomatics) <- c("snvs", "indels", "svs", "age", "cancerType", "source")



#############  Combine into MEGA data frame
purpleIndels = purpleVariants[, c("indels", "cancerType", "age", "source")]
colnames(purpleIndels) <- c("value", "cancerType", "age", "source")
purpleIndels$type <- "Purple Indels"

purpleSnvs = purpleVariants[, c("snvs", "cancerType", "age", "source")]
colnames(purpleSnvs) <- c("value", "cancerType", "age", "source")
purpleSnvs$type <- "Purple SNVs"

purpleSvs = purpleVariants[, c("svs", "cancerType", "age", "source")]
colnames(purpleSvs) <- c("value", "cancerType", "age", "source")
purpleSvs$type <- "Purple SVs"

pcawgIndels = pcawgSomatics[, c("indels", "cancerType", "age", "source")]
colnames(pcawgIndels) <- c("value", "cancerType", "age", "source")
pcawgIndels$type <- "PCAWG Indels"

pcawgSnvs = pcawgSomatics[, c("snvs", "cancerType", "age", "source")]
colnames(pcawgSnvs) <- c("value", "cancerType", "age", "source")
pcawgSnvs$type <- "PCAWG SNVs"

pcawgSvs = pcawgSomatics[, c("svs", "cancerType", "age", "source")]
colnames(pcawgSvs) <- c("value", "cancerType", "age", "source")
pcawgSvs$type <- "PCAWG SVs"

combined = rbind(purpleSvs, purpleIndels, purpleSnvs, pcawgSvs, pcawgIndels, pcawgSnvs)

rm(purpleIndels)
rm(purpleSnvs)
rm(purpleSvs)
rm(pcawgIndels)
rm(pcawgSnvs)
rm(pcawgSvs)

#############  Indel, SNV and SV CDFs
mappedCancerTypes = unique(pcawgSomatics[!is.na(pcawgSomatics$cancerType), c("cancerType")])

ggplot(data=subset(combined, cancerType %in% mappedCancerTypes), aes(x=value, color = type)) +
  stat_ecdf(geom = "step", pad = FALSE) +
  facet_grid(cancerType~., scales = "free") +
  scale_x_log10()

#############  Indels by age
ggplot(data=subset(combined, cancerType %in% mappedCancerTypes & type %in% c('PCAWG Indels', "Purple Indels")),  aes(x=age, y=value, color = type)) +
  geom_point() +
  facet_grid(cancerType~.) +
  scale_y_log10()  +
  labs(title="Indel count by age")

#############  SNVs by age
ggplot(data=subset(combined, cancerType %in% mappedCancerTypes & type %in% c('PCAWG SNVs', "Purple SNVs")),  aes(x=age, y=value, color = type)) +
  geom_point() +
  facet_grid(cancerType~.) +
  scale_y_log10()  +
  labs(title="SNVs by age")

#ggplot(data=subset(combined, cancerType %in% mappedCancerTypes),
#       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source))  + facet_grid(cancerType ~ type, scales = "free") + scale_x_log10()


#############  EXPERIMENTAL RUBBISH BELOW

#ggplot(data=subset(combined, cancerType %in% c('Breast', 'Melanoma')),
#       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast")  + facet_wrap( ~type, scales = "free")


#############  PLOT INDELS
#combinedSomatics = rbind(pcawgSomatics, purpleVariants)

#ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')),
#      aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast Indels")
#
#ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')),
#       aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma Indels")



#ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')),
#       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast SNVs")

#ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')),
#       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma SNVs")

#aes(chCopyNumber)
#ggplot(aes(indels), data=allSomatics)+ stat_ecdf(geom = "step", pad = FALSE) + facet_wrap( ~cancerType )


#purpleBreast = allSomatics[allSomatics$cancerType == 'Breast',]
#purpleLung = allSomatics[allSomatics$cancerType == 'Lung',]

#ggplot(aes(indels), data=purpleBreast) + stat_ecdf(geom = "step", pad = FALSE)
#ggplot(aes(indels), data=purpleLung) + stat_ecdf(geom = "step", pad = FALSE)



#gplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')),
#       aes(x=indels, color=cancerType))+ stat_ecdf(geom = "step", pad = FALSE)

#ggplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')),
#       aes(x=indels))+ stat_ecdf(geom = "step", pad = FALSE, aes(color = cancerType))


