library(RMySQL)
library(ggplot2)
library(dplyr)

extractCancerTypeCounts <- function(cancerTypes, copyNumberTargets) {
  cancerTypeColumns = gsub(" ", ".", cancerTypes)
  cancerTypeColumns = ifelse(is.na(cancerTypeColumns), "NA.", cancerTypeColumns)
  availableCancerColumns = intersect(cancerTypeColumns, names(copyNumberTargets))
  result = copyNumberTargets %>% select(availableCancerColumns)
  result[is.na(result)] <- 0
  return (result)
}

tidyTargets <- function(copyNumberTargets, maxTargets = 30) {
  names = names(copyNumberTargets)
  namesInd = setNames(1:ncol(copyNumberTargets), names)

  topTargets = copyNumberTargets %>% group_by(target) %>% summarise(N = sum(N)) %>% top_n(maxRecords, N) %>% arrange(-N)
  tidyCopyNumberTargets = copyNumberTargets %>% 
    ungroup() %>%
    filter(target %in% topTargets$target) %>%
    mutate(gene = factor(target, levels=topTargets$target)) %>% 
    select(gene, c(namesInd["Biliary"]:namesInd["Uterus"])) %>%
    gather(cancerType, N, -gene) %>%
    filter(!is.na(N)) %>%
    group_by(gene, cancerType) %>%
    summarise(N = sum(N))

  return (tidyCopyNumberTargets)
}

###################### PART 1 ###################### 

load("~/hmf/RData/reference/cancerTypeColours.RData")
load("~/hmf/RData/processed/geneCopyNumberDeleteTargets.RData")
load("~/hmf/RData/processed/geneCopyNumberAmplificationTargets.RData")


significantDels = tidyTargets(geneCopyNumberDeleteTargets)
significantAmps = tidyTargets(geneCopyNumberAmplificationTargets, maxTargets = 40)

ampsPlot = ggplot(data=significantAmps, aes(gene, N)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Top Gene Amplifications") + xlab("Genes") + ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual( values= cancerTypeColours)

delsPlot = ggplot(data=significantDels, aes(gene, N)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Top Gene Deletions") + xlab("Genes") + ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual( values= cancerTypeColours) 


purple::multiplot(ampsPlot, delsPlot, cols = 1)

#cohortByCancerType = cohort %>% group_by(cancerType) %>% summarise(N = n())
#save(cohortByCancerType, file = '~/hmf/RData/cohortByCancerType.RData')

###################### PART 4 ###################### 
rm(list=ls())

containsCandidate <- function(target, category) {
  splitCategory = strsplit(category, ",")
  result = rep(F, length(target))
  for (i in 1: length(result)) {
    result[i] = target[i] %in% splitCategory[[i]]
  }
  
  return (result)
}

categoryStatus <- function(copyNumberTargets, del) {
  copyNumberTargets$status <- NA
  if (del) {
    copyNumberTargets$status <- ifelse(containsCandidate(copyNumberTargets$target, copyNumberTargets$fragile), "Fragile", NA)
    copyNumberTargets$status <- ifelse(is.na(copyNumberTargets$status) & containsCandidate(copyNumberTargets$target, copyNumberTargets$cosmicTsg), "Known", copyNumberTargets$status)
  } else {
    copyNumberTargets$status <- ifelse(is.na(copyNumberTargets$status) & containsCandidate(copyNumberTargets$target, copyNumberTargets$cosmicOncogene), "Known", copyNumberTargets$status)
  }
  copyNumberTargets$status <- ifelse(is.na(copyNumberTargets$status) & copyNumberTargets$method == "panel", "Significant", copyNumberTargets$status)
  copyNumberTargets$status <- ifelse(is.na(copyNumberTargets$status), "Unknown", copyNumberTargets$status)
}


countByCancerType <- function(cohortByCancerType, copyNumberTargets, del) {
  colnames(cohortByCancerType) <- c("cancerType", "cancerTypeSamples")
  
  cancerCounts = extractCancerTypeCounts(cohortByCancerType$cancerType, copyNumberTargets)
  status = categoryStatus(copyNumberTargets, del)
  copyNumberTargets = copyNumberTargets %>% select(gene, target, N)
  copyNumberTargets = cbind(copyNumberTargets, status) 
  copyNumberTargets = cbind(copyNumberTargets, cancerCounts) 
  
  tidyCopyNumberTargets = copyNumberTargets %>% gather(cancerType, N, c(5:ncol(copyNumberTargets))) %>% group_by(cancerType, status) %>% summarise(N = sum(N))
  tidyCopyNumberTargets$cancerType = gsub("\\.", " ", tidyCopyNumberTargets$cancerType)
  tidyCopyNumberTargets$cancerType = ifelse(tidyCopyNumberTargets$cancerType == "NA ", NA, tidyCopyNumberTargets$cancerType)
  
  tidyCopyNumberTargets = dplyr::inner_join(tidyCopyNumberTargets, cohortByCancerType, by = "cancerType")
  tidyCopyNumberTargets$rate = tidyCopyNumberTargets$N / tidyCopyNumberTargets$cancerTypeSamples
  
  tidyCopyNumberCancerTypeLevels = tidyCopyNumberTargets %>% group_by(cancerType) %>% summarise(rate = sum(rate)) %>% arrange(-rate) %>% select(cancerType)
  tidyCopyNumberTargets$cancerType = factor(tidyCopyNumberTargets$cancerType, levels = tidyCopyNumberCancerTypeLevels$cancerType)
  
  return (tidyCopyNumberTargets)
}

#load("~/hmf/RData/cancerTypes.RData")
#load("~/hmf/RData/geneCopyNumberDeleteTargets.RData")
#load("~/hmf/RData/geneCopyNumberAmplificationTargets.RData")
load('~/hmf/RData/cohortByCancerType.RData')
tidyDeletesByCancerType = countByCancerType(cohortByCancerType, geneCopyNumberDeleteTargets %>% filter(N >= 5), del = T)
#tidyDeletesByCancerType = dplyr::inner_join(tidyDeletesByCancerType, cohortByCancerType, by = "cancerType")
#tidyDeletesByCancerType$rate = tidyDeletesByCancerType$N / tidyDeletesByCancerType$cancerTypeSamples

p1 = ggplot(data=tidyDeletesByCancerType, aes(cancerType, rate)) +
  geom_bar(aes(fill = status), stat = "identity") + 
  ggtitle("Gene deletions") + xlab("Cancer Type") + ylab("Rate per sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

tidyAmpsByCancerType = countByCancerType(cohortByCancerType, geneCopyNumberAmplificationTargets %>% filter(N > 15), del = F)
#tidyAmpsByCancerType = dplyr::inner_join(tidyAmpsByCancerType, cohortByCancerType, by = "cancerType")
#tidyAmpsByCancerType$rate = tidyAmpsByCancerType$N / tidyAmpsByCancerType$cancerTypeSamples

p2 = ggplot(data=tidyAmpsByCancerType, aes(cancerType, rate)) +
  geom_bar(aes(fill = status), stat = "identity") + 
  ggtitle("Gene amplifications") + xlab("Cancer Type") + ylab("Rate per sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

multiplot(p2, p1, cols = 1)
  
###################### PART 2 ###################### 
rm(list = ls())
load("~/hmf/RData/amps.RData")
load("~/hmf/RData/dels.RData")

ampsSharedWithPcawg = amps %>% filter(pcawgGeneCount > 0) %>% select(gene_name, hmfCount, pcawgGeneCount) %>% 
  group_by(gene_name)  %>% summarise(hmfCount = sum(hmfCount), pcawgGeneCount = max(pcawgGeneCount)) %>% 
  mutate(total = hmfCount + pcawgGeneCount) %>% arrange(-total)
ampsSharedWithPcawg$hmfType <- "HMF"
ampsSharedWithPcawg$pcawgType <- "PCAWG"
ampsSharedWithPcawg$gene_name <- factor(ampsSharedWithPcawg$gene_name, levels = ampsSharedWithPcawg$gene_name)

tidySharedAmps = rbind(ampsSharedWithPcawg %>% select(gene_name, count = pcawgGeneCount, type = pcawgType), ampsSharedWithPcawg %>% select(gene_name, count = hmfCount, type = hmfType))
p1 = ggplot(data=tidySharedAmps, aes(gene_name, count)) +
  geom_bar(aes(fill = type), stat = "identity") +
  ggtitle("PCAWG Amplifications") + xlab("Genes") + ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

delsSharedWithPcawg = dels %>% filter(pcawgGeneCount > 0) %>% select(gene_name, hmfCount, pcawgGeneCount) %>% 
  group_by(gene_name)  %>% summarise(hmfCount = sum(hmfCount), pcawgGeneCount = max(pcawgGeneCount)) %>% 
  mutate(total = hmfCount + pcawgGeneCount) %>% arrange(-total)
delsSharedWithPcawg$hmfType <- "HMF"
delsSharedWithPcawg$pcawgType <- "PCAWG"
delsSharedWithPcawg$gene_name <- factor(delsSharedWithPcawg$gene_name, levels = delsSharedWithPcawg$gene_name)


tidySharedDels = rbind(delsSharedWithPcawg %>% select(gene_name, count = pcawgGeneCount, type = pcawgType), delsSharedWithPcawg %>% select(gene_name, count = hmfCount, type = hmfType))
p2 = ggplot(data=tidySharedDels, aes(gene_name, count)) +
  geom_bar(aes(fill = type), stat = "identity") +
  ggtitle("PCAWG Deletions") + xlab("Genes") + ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

multiplot(p1, p2, cols = 1)

ampsSharedWithPcawg$cnv <- "Amplifications"
delsSharedWithPcawg$cnv <- "Deletions"
sharedWithPcawg = rbind(ampsSharedWithPcawg, delsSharedWithPcawg)
ggplot(data=sharedWithPcawg,aes(hmfCount, pcawgGeneCount,label=gene_name, colour=factor(cnv)))+
  geom_point() +
  geom_text(size=2,hjust = 0, nudge_x = 0.03)+labs(title='Gene Amplifications') +
  scale_x_log10() +
  scale_y_log10() 



###################### PART 4 ###################### 




