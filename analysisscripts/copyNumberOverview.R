#install.packages("devtools")
#devtools::install_github("hadley/multidplyr")

detach("package:purple", unload=TRUE)
library(purple)
library(tidyr)
library(dplyr)
library(multidplyr)
library(doParallel)
library(GenomicRanges)

isMaleSexChromosome <- function(gender, chromosome) {
  return (gender != "FEMALE" & chromosome %in% c('X','Y'))
}

isLoh <- function(minAllelePloidy, chromosome, gender) {
  result <- ifelse(isMaleSexChromosome(gender, chromosome), FALSE, minAllelePloidy < 0.5)
  return (result)
}

isDel <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
  result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber < cutoff * ploidy / 2, copyNumber < cutoff * ploidy)
  return (result)
}

isAmp <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
  result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber > cutoff * ploidy / 2,  copyNumber > cutoff * ploidy)
  return (result)
}


load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
highestPurityCohort$gender <- ifelse(substr(highestPurityCohort$gender, 1,4) == 'MALE', 'MALE', highestPurityCohort$gender)

load(file = "~/hmf/RData/Reference/hpcCopyNumbers.RData")
highestPurityCopyNumbers = hpcCopyNumbers %>%
  left_join(highestPurityCohort[, c("sampleId", "ploidy", "cancerType", "gender")], by = "sampleId") %>%
  mutate(minAllelePloidy = pmax(0, (1-actualBaf) * copyNumber),
    loh = isLoh(minAllelePloidy, chromosome, gender),
    absDel = isDel(copyNumber, chromosome, gender, 1, 0.5),
    relDel = loh & isDel(copyNumber, chromosome, gender, ploidy, 0.6),
    amp1_4 = isAmp(copyNumber, chromosome, gender, ploidy, 1.4),
    amp2_0 = isAmp(copyNumber, chromosome, gender, ploidy, 2),
    amp3_0 =isAmp(copyNumber, chromosome, gender, ploidy, 3))

highestPurityCopyNumbers$region = GRanges(paste0("chr", highestPurityCopyNumbers$chromosome), ranges = IRanges(start = highestPurityCopyNumbers$start, end = highestPurityCopyNumbers$end))

library(BSgenome.Hsapiens.UCSC.hg19)
bins <- tileGenome(seqinfo(Hsapiens), tilewidth=300000, cut.last.tile.in.chrom=TRUE)
length(bins)

ol = as.matrix(findOverlaps(bins, highestPurityCopyNumbers$region, type = "any"))
olCopyNumbers = cbind(bin = ol[, 1], highestPurityCopyNumbers[ol[, 2], c("sampleId", "cancerType", "loh", "absDel", "relDel", "amp1_4", "amp2_0", "amp3_0")])

no_cores <- 7
cl<-makeCluster(no_cores, type="FORK"); date()
binSampleSummary = olCopyNumbers %>% partition(bin, sampleId, cluster = cl) %>%
  summarise(cancerType = dplyr::first(cancerType),
            loh = any(loh), absDel = any(absDel), relDel = any(relDel), amp1_4 = any(amp1_4), amp2_0 = any(amp2_0), amp3_0 = any(amp3_0))  %>%
  collect() %>%
  as_tibble()
date()
binCancerTypeSummary = binSampleSummary %>% partition(bin, cancerType, cluster = cl) %>%
  summarise(loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0)) %>%
  collect() %>%
  as_tibble()
date()
stopCluster(cl)

binSummary = binSampleSummary %>% group_by(bin)  %>%  summarise(cancerType = 'All', loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0))
binSummary = bind_rows(binSummary, binCancerTypeSummary)

primaryTumorLocationCounts = highestPurityCohort %>% group_by(cancerType) %>% summarise(nTotal = as.numeric(n()))
primaryTumorLocationMaleCounts = highestPurityCohort %>% filter(gender == 'MALE') %>% group_by(cancerType) %>% summarise(nMale = as.numeric(n()))
primaryTumorLocationCounts = merge(primaryTumorLocationCounts, primaryTumorLocationMaleCounts, by = 'cancerType', all = T)
primaryTumorLocationCounts$nMale <- ifelse(is.na(primaryTumorLocationCounts$nMale), 0, primaryTumorLocationCounts$nMale)
allCounts = highestPurityCohort %>% summarise(cancerType = "All", nTotal = as.numeric(n()))
allCountsMale = highestPurityCohort %>% filter(gender == 'MALE') %>% summarise(cancerType = "All", nMale = as.numeric(n()))
allCounts = merge(allCounts, allCountsMale, by = 'cancerType', all = T)
primaryTumorLocationCounts = rbind(primaryTumorLocationCounts, allCounts)
rm(primaryTumorLocationMaleCounts, allCounts, allCountsMale)


highestPurityCopyNumberSummary = left_join(binSummary, primaryTumorLocationCounts, by = "cancerType")
highestPurityCopyNumberSummary$region = bins[highestPurityCopyNumberSummary$bin]
highestPurityCopyNumberSummary$chromosome = substring(as.character(seqnames(highestPurityCopyNumberSummary$region)), 4)
highestPurityCopyNumberSummary$start = start(highestPurityCopyNumberSummary$region)
highestPurityCopyNumberSummary$end = end(highestPurityCopyNumberSummary$region)
highestPurityCopyNumberSummary$region <- NULL
highestPurityCopyNumberSummary = highestPurityCopyNumberSummary %>%
  mutate(n = ifelse(chromosome == 'Y', nMale, nTotal)) %>%
  mutate(lohPercentage = loh /n, absDelPercentage = absDel / n, relDelPercentage = relDel / n, amp1_4Percentage = amp1_4 / n,  amp2_0Percentage = amp2_0 / n, amp3_0Percentage = amp3_0 / n) %>%
  ungroup() %>%
  select(chromosome, start, end, cancerType, lohPercentage, absDelPercentage, relDelPercentage, amp1_4Percentage, amp2_0Percentage, amp3_0Percentage)
save(highestPurityCopyNumberSummary, file = "~/hmf/RData/processed/hpcCopyNumberSummary.RData")




load(file = "~/hmf/RData/processed/hpcCopyNumberSummary.RData")
primaryTumorLocations = unique(highestPurityCopyNumberSummary$cancerType)
for (location in primaryTumorLocations) {

  locationString = gsub(" ", "", location, fixed = TRUE)
  locationString = gsub("/", "", locationString, fixed = TRUE)

  primaryTumorLocationSummary = highestPurityCopyNumberSummary %>%
    filter(cancerType == location) %>%
    mutate(chromosome = paste0("hs", chromosome))

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, lohPercentage) %>% mutate(lohPercentage = -lohPercentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".loh.circos"), sep = "\t", row.names = F, col.names = F, quote = F)

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, absDelPercentage) %>% mutate(absDelPercentage = -absDelPercentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".absDel.circos"), sep = "\t", row.names = F, col.names = F, quote = F)

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, relDelPercentage) %>% mutate(relDelPercentage = -relDelPercentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".relDel.circos"), sep = "\t", row.names = F, col.names = F, quote = F)

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp1_4Percentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".amp1_4.circos"), sep = "\t", row.names = F, col.names = F, quote = F)

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp2_0Percentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".amp2_0.circos"), sep = "\t", row.names = F, col.names = F, quote = F)

  write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp3_0Percentage),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".amp3_0.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
}

tidyTargets <- function(copyNumberTargets) {
  names = names(copyNumberTargets)
  namesInd = setNames(1:ncol(copyNumberTargets), names)
  
  topTargets = copyNumberTargets %>% group_by(target) %>% summarise(N = sum(N)) %>% arrange(-N)
  tidyCopyNumberTargets = copyNumberTargets %>% 
    ungroup() %>%
    filter(target %in% topTargets$target) %>%
    mutate(gene = target) %>%
    select(gene, chromosome, start, end, c(namesInd["Biliary"]:namesInd["Uterus"])) %>%
    gather(cancerType, N, -gene, -chromosome, -start, -end) %>%
    filter(!is.na(N)) %>%
    group_by(gene, chromosome, cancerType) %>%
    summarise(N = sum(N), start = min(start), end = max(end))
  
  return (tidyCopyNumberTargets)
}


load(file = "~/hmf/RData/reference/canonicalTranscripts.RData")
load(file = "~/hmf/RData/processed/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/processed/geneCopyNumberDeleteTargets.RData")

geneCopyNumberDeleteTargets$telomere <- gsub("_telomere","", geneCopyNumberDeleteTargets$telomere)
geneCopyNumberDeleteTargets$centromere <- gsub("_centromere","", geneCopyNumberDeleteTargets$centromere)

geneCopyNumberDeleteTargets = geneCopyNumberDeleteTargets %>% mutate(target =  coalesce(telomere, centromere, target))
tidyDelTargets = tidyTargets(geneCopyNumberDeleteTargets) %>% ungroup() %>% mutate(chromosome = paste0("hs", chromosome), info = "color=vdred")
tidyDelTargetsAll = tidyDelTargets %>% group_by(gene, chromosome, start, end, info) %>% summarise(N = sum(N)) %>% mutate(cancerType = 'All')
tidyDelTargets = bind_rows(tidyDelTargetsAll, tidyDelTargets)

tidyAmpTargets = tidyTargets(geneCopyNumberAmplificationTargets) %>% ungroup() %>% mutate(chromosome = paste0("hs", chromosome), info = "color=black")
tidyAmpTargetsAll = tidyAmpTargets %>% group_by(gene, chromosome, start, end, info) %>% summarise(N = sum(N)) %>% mutate(cancerType = 'All')
tidyAmpTargets = bind_rows(tidyAmpTargets, tidyAmpTargetsAll)

label_size <- function(genes) {
  label_size_single <- function(gene) {
    if (nchar(gene) < 7) {
      return ("label_size=40p")
    }
    
    if (nchar(gene) < 10) {
      return ("label_size=30p")
    }
    
    return ("label_size=24p")
  }
  
  sapply(genes,label_size_single)
  
}


for (location in primaryTumorLocations) {
  locationString = gsub(" ", "", location, fixed = TRUE)
  locationString = gsub("/", "", locationString, fixed = TRUE)
  
  targetDel = tidyDelTargets %>% ungroup() %>% filter(N >= 3, cancerType == location) %>% top_n(10, N)
  targetAmp = tidyAmpTargets %>% ungroup() %>% filter(N >= 3, cancerType == location) %>% top_n(10, N)
  genes = bind_rows(targetDel, targetAmp)
  genes = genes %>% mutate(size = label_size(gene)) %>%
    mutate(info = paste(info, size, sep = ","))
  
  
  write.table(genes %>% select(chromosome, start, end, gene, info),
              file = paste0("~/hmf/analysis/copyNumberSummary/", locationString, ".genes.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
}


for (location in primaryTumorLocations) {
  locationString = gsub(" ", "", location, fixed = TRUE)
  locationString = gsub("/", "", locationString, fixed = TRUE)

  cmd = paste0("sed 's/CANCER/", locationString,"/g' circos.template > ",locationString,".conf\n")
  #cmd = paste0("/Users/jon/hmf/tools/circos-0.69-6/bin/circos -nosvg -conf /Users/jon/hmf/analysis/copyNumberSummary/",locationString, ".conf -outputdir /Users/jon/hmf/analysis/copyNumberSummary -outputfile ",locationString,".png\n")
  cat(cmd)
}

