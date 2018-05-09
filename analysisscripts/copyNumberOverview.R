#install.packages("devtools")
#devtools::install_github("hadley/multidplyr")

detach("package:purple", unload=TRUE)
library(purple)
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

load(file = "~/hmf/RData/reference/highestPurityCopyNumbers.RData")
highestPurityCopyNumbers = highestPurityCopyNumbers %>%
  left_join(highestPurityCohort[, c("sampleId", "ploidy", "primaryTumorLocation", "gender")], by = "sampleId") %>%
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
olCopyNumbers = cbind(bin = ol[, 1], highestPurityCopyNumbers[ol[, 2], c("sampleId", "primaryTumorLocation", "loh", "absDel", "relDel", "amp1_4", "amp2_0", "amp3_0")])

no_cores <- 7
cl<-makeCluster(no_cores, type="FORK"); date()
binSampleSummary = olCopyNumbers %>% partition(bin, sampleId, cluster = cl) %>%
  summarise(primaryTumorLocation = first(primaryTumorLocation),
            loh = any(loh), absDel = any(absDel), relDel = any(relDel), amp1_4 = any(amp1_4), amp2_0 = any(amp2_0), amp3_0 = any(amp3_0))  %>%
  collect() %>%
  as_tibble()
date()
binCancerTypeSummary = binSampleSummary %>% partition(bin, primaryTumorLocation, cluster = cl) %>%
  summarise(loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0)) %>%
  collect() %>%
  as_tibble()
date()
stopCluster(cl)

binSummary = binSampleSummary %>% group_by(bin)  %>%  summarise(primaryTumorLocation = 'All', loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0))
binSummary = bind_rows(binSummary, binCancerTypeSummary)

primaryTumorLocationCounts = highestPurityCohort %>% group_by(primaryTumorLocation) %>% summarise(nTotal = as.numeric(n()))
primaryTumorLocationMaleCounts = highestPurityCohort %>% filter(gender == 'MALE') %>% group_by(primaryTumorLocation) %>% summarise(nMale = as.numeric(n()))
primaryTumorLocationCounts = merge(primaryTumorLocationCounts, primaryTumorLocationMaleCounts, by = 'primaryTumorLocation', all = T)
primaryTumorLocationCounts$nMale <- ifelse(is.na(primaryTumorLocationCounts$nMale), 0, primaryTumorLocationCounts$nMale)
allCounts = highestPurityCohort %>% summarise(primaryTumorLocation = "All", nTotal = as.numeric(n()))
allCountsMale = highestPurityCohort %>% filter(gender == 'MALE') %>% summarise(primaryTumorLocation = "All", nMale = as.numeric(n()))
allCounts = merge(allCounts, allCountsMale, by = 'primaryTumorLocation', all = T)
primaryTumorLocationCounts = rbind(primaryTumorLocationCounts, allCounts)
rm(primaryTumorLocationMaleCounts, allCounts, allCountsMale)


highestPurityCopyNumberSummary = left_join(binSummary, primaryTumorLocationCounts, by = "primaryTumorLocation")
highestPurityCopyNumberSummary$region = bins[highestPurityCopyNumberSummary$bin]
highestPurityCopyNumberSummary$chromosome = substring(as.character(seqnames(highestPurityCopyNumberSummary$region)), 4)
highestPurityCopyNumberSummary$start = start(highestPurityCopyNumberSummary$region)
highestPurityCopyNumberSummary$end = end(highestPurityCopyNumberSummary$region)
highestPurityCopyNumberSummary$region <- NULL
highestPurityCopyNumberSummary = highestPurityCopyNumberSummary %>%
  mutate(n = ifelse(chromosome == 'Y', nMale, nTotal)) %>%
  mutate(lohPercentage = loh /n, absDelPercentage = absDel / n, relDelPercentage = relDel / n, amp1_4Percentage = amp1_4 / n,  amp2_0Percentage = amp2_0 / n, amp3_0Percentage = amp3_0 / n) %>%
  ungroup() %>%
  select(chromosome, start, end, primaryTumorLocation, lohPercentage, absDelPercentage, relDelPercentage, amp1_4Percentage, amp2_0Percentage, amp3_0Percentage)
save(highestPurityCopyNumberSummary, file = "~/hmf/RData/processed/highestPurityCopyNumberSummary.RData")




load(file = "~/hmf/RData/processed/highestPurityCopyNumberSummary.RData")
primaryTumorLocations = unique(highestPurityCopyNumberSummary$primaryTumorLocation)
primaryTumorLocations = primaryTumorLocations[!is.na(primaryTumorLocations)]
primaryTumorLocations
for (location in primaryTumorLocations) {

  locationString = gsub(" ", "", location, fixed = TRUE)
  locationString = gsub("/", "", locationString, fixed = TRUE)

  primaryTumorLocationSummary = highestPurityCopyNumberSummary %>%
    filter(primaryTumorLocation == location) %>%
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



for (location in primaryTumorLocations) {
  locationString = gsub(" ", "", location, fixed = TRUE)
  locationString = gsub("/", "", locationString, fixed = TRUE)

  cmd = paste0("sed 's/CANCER/", locationString,"/g' circos.template > ",locationString,".conf\n")
  cat(cmd)

  #cmd = paste0("/Users/jon/hmf/tools/circos-0.69-6/bin/circos -nosvg -conf /Users/jon/hmf/analysis/copyNumberSummary/",locationString, ".conf -outputdir /Users/jon/hmf/analysis/copyNumberSummary -outputfile ",locationString,".png\n")
  #cat(cmd)
}

