library(dplyr)
library(tidyr)
library(RMySQL)
library(ggplot2)

troughIndex <- function(y) {
  result = c()
  for (i in c(2:(length(y) - 1))) {
    if (y[i-1] > y[i] & y[i] < y[i+1]) {
      result= c(result, i)
    }
  }
  return (result)
}

peakIndex <- function(y) {
  result = c()
  for (i in c(2:(length(y) - 1))) {
    if (y[i-1] < y[i] & y[i] > y[i+1]) {
      result= c(result, i)
    }
  }
  return (result)
}

peaks <- function(sampleId, x, y) {
  
  result = data.frame()
  
  peakHeight = NA
  peakPloidy = NA
  troughHeight = NA
  troughPloidy = NA
  
  for (i in c(2:(length(y) - 1))) {
    if (y[i-1] < y[i] & y[i] > y[i+1]) {
      #cat("Found peak at index", i, "\n")
      peakHeight = round(y[i], 3)
      peakPloidy = round(x[i], 3)
    } 
    
    if (y[i-1] > y[i] & y[i] <= y[i+1]) {
      prevTroughHeight = troughHeight;
      prevTroughPloidy = troughPloidy
      troughHeight = round(y[i],3)
      troughPloidy = round(x[i],3)
      
      if (!is.na(peakPloidy)) {
        sampleResult = data.frame(sampleId = c(sampleId), peakPloidy = c(peakPloidy), peakHeight = c(peakHeight), 
                                  prevTroughPloidy = c(prevTroughPloidy), prevTroughHeight = c(prevTroughHeight),
                                  nextTroughPloidy = c(troughPloidy), nextTroughHeight = c(troughHeight), stringsAsFactors = F)
        result = bind_rows(result, sampleResult)
      }
    } 
    
  }
  if (!is.na(peakPloidy)) {
    sampleResult = data.frame(sampleId = c(sampleId), peakPloidy = c(peakPloidy), peakHeight = c(peakHeight), 
                              prevTroughPloidy = c(troughPloidy), prevTroughHeight = c(troughHeight),
                              nextTroughPloidy = c(NA), nextTroughHeight = c(NA), stringsAsFactors = F)
    result = bind_rows(result, sampleResult)
  }
  
  result = result %>% mutate(order = row_number()) %>% arrange(-peakHeight) %>% 
    mutate(rankHeight = row_number())  %>% select(sampleId, order, rankHeight, everything()) %>% arrange(order)
  
  return (result)
}

somatic_peaks2 <- function(somatics) {
  result = data.frame()
  for (sample in unique(somatics$sampleId)) {
    densityData = somatics %>% 
      filter(sampleId == sample) %>% 
      mutate(ploidy = pmax(0, adjustedCopyNumber * adjustedVaf)) %>%
      filter(ploidy < 4) 
    
    d <- density(densityData$ploidy, bw = 0.04, n = 100)
    result = bind_rows(result, peaks(sample, d$x, d$y))
  }
  
  return (result)
}

somatic_peaks_by_chromosome <- function(somatics) {
  result = data.frame()
  for (sample in unique(somatics$sampleId)) {
    cat(sample, "\n")
    
    densityData = somatics %>% 
      filter(sampleId == sample, !chromosome %in% c('Y','MT')) %>% 
      mutate(ploidy = pmax(0, adjustedCopyNumber * adjustedVaf)) %>%
      filter(ploidy < 4) 
    

    for (chrom in unique(densityData$chromosome)) {
      dd2 = densityData %>% filter(chromosome == chrom) 
      d <- density(dd2$ploidy, bw = 0.04, n = 100)
      chromResult = peaks(sample, d$x, d$y)
      if (nrow(chromResult) > 0) {
        chromResult$chromosome = chrom
        result = bind_rows(result, chromResult)
      }
    }
    
  }
  
  return (result)
}

allSomaticPeaks_p2 = somatic_peaks2(allSomatics_p2)
save(allSomaticPeaks_p2, file = "~/hmf/RData/allSomaticPeaks_p2.RData")

allSomaticPeaksByChromosome_p2 = somatic_peaks_by_chromosome(allSomatics_p2)

d = density(densityData$ploidy, bw = 0.04, n = 100)
peaks("CPCT02060182T", d$x, d$y)
peakIndex(d$y)
d$x
d$y

View(allSomaticPeaks_p2 %>% filter(rankHeight == 1))
d = density(densityData$ploidy, bw = 0.04, n = 100)
jon = peaks("CPCT02180057T", d$x, d$y)

collapsed_values <- function(indexes, values) {
  return (paste0(sapply(indexes, function (x) {round(values[x], 3)}), collapse = ","))
}

somatic_peaks <- function(somatics) {
  result = data.frame()
  
  for (sample in unique(somatics$sampleId)[1:10]) {
    densityData = somatics %>% 
      filter(sampleId == sample) %>% 
      mutate(ploidy = pmax(0, adjustedCopyNumber * adjustedVaf)) %>%
      filter(ploidy < 4) %>%
      pull(ploidy)
    
    d <- density(densityData, bw = 0.04, n = 100)
    troughIndexes = troughIndex(d$y)
    troughPloidies = collapsed_values(troughIndexes, d$x)
    troughHeights = collapsed_values(troughIndexes, d$y)
    
    peakIndexes = peakIndex(d$y)
    peakPloidies = collapsed_values(peakIndexes, d$x)
    peakHeights = collapsed_values(peakIndexes, d$y)
    
    sampleResult = data.frame(sampleId = sample, peakPloidy = peakPloidies, peakHeight = peakHeights, troughPloidy = troughPloidies, troughHeight = troughHeights)
    result = bind_rows(sampleResult, result)
   }
  
  return (result)
}

rm(troughIndexes, troughPloidies, troughHeights, peakIndexes, peakPloidies, peakHeights, d, allSomaticPeaks_p1, allSomaticPeaks_p2, densityData)
rm(jon, jon1, jon2, sampleResult, sample, somatics)

load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
allSomaticPeaks_p1 = somatic_peaks(allSomatics_p1)
save(allSomaticPeaks_p1, file = "~/hmf/RData/allSomaticPeaks_p1.RData")
rm(allSomatics_p1)

load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
allSomaticPeaks_p2 = somatic_peaks(allSomatics_p2)
save(allSomaticPeaks_p2, file = "~/hmf/RData/allSomaticPeaks_p2.RData")
rm(allSomatics_p2)





####### RUBBISH
unique(allSomatics_p2$sampleId)
sample = "CPCT02180041T"
sample = "CPCT02190023T"
sample = "CPCT02170022T"
sample = "CPCT02170014T"
sample = "CPCT02170016T"

sample = "CPCT02120090T"
sample = "CPCT02120091T"
sample = "CPCT02120099T"
sample = "CPCT02010450TII"


1   CPCT02010618T 3.425147   0.39
2   DRUP01010054T 3.019199   0.39
3   DRUP01010078T 2.779752   0.39
4 CPCT02010450TII 3.641782   0.38
5   CPCT02040117T 3.158549   0.38
6   CPCT02040198T 2.594004   0.37
7   CPCT02040174T 1.735438   0.35
8   CPCT02190031T 1.923636   0.32

sample = "CPCT02060150T"
sample = "CPCT02010679T"
#load("~/hmf/RData/Reference/allPurity.RData")
#highestPurityCohort = highestPurityCohort %>% filter(sampleId %in% allSomatics_p2$sampleId)
#View(highestPurityCohort %>% filter(ploidy>2.5))

#highestPurityCohort %>% filter(cancerType == "Lung", purity < 0.4) %>% arrange(-purity) %>% select(sampleId, ploidy, purity)
densityData = allSomatics_p2 %>% 
  filter(sampleId == sample) %>% 
  mutate(
    ploidy = pmax(0, adjustedCopyNumber * adjustedVaf),
    majorAllelePloidy = adjustedCopyNumber - minorAllelePloidy,
    minorAllelePloidyDistance = abs(minorAllelePloidy - round(minorAllelePloidy, 0)),
    majorAllelePloidyDistance = abs(majorAllelePloidy - round(majorAllelePloidy, 0))
    ) %>%
  filter(!chromosome %in%  c('Y','MT'), ploidy < 4)
  #filter(ploidy < 4, abs(adjustedCopyNumber - round(adjustedCopyNumber, 0) ) < 0.2)

ggplot(data = densityData, aes(minorAllelePloidyDistance)) + 
  geom_histogram(aes(y =..density..), bins = 200, fill="blue", col="blue",  alpha = .4) + 
  geom_histogram(aes(x = majorAllelePloidyDistance, y =..density..), bins = 200, fill="red", col="red",  alpha = .4)+
  xlim(-0.05, 0.5) + facet_wrap(~chromosome)

ggplot(data = densityData, aes(minorAllelePloidy)) + 
  geom_histogram(aes(y =..density..), bins = 200, fill="blue", col="blue",  alpha = .4) + 
  geom_histogram(aes(x = adjustedCopyNumber - minorAllelePloidy, y =..density..), bins = 200, fill="red", col="red",  alpha = .4)+
  xlim(-1, 6) #+ facet_wrap(~chromosome)

#highestPurityCohort %>% filter(sampleId == sample) %>% select(purity)

ggplot(data = densityData, aes(ploidy)) + 
  geom_histogram(aes(y =..density..), bins = 200, fill="blue", col="blue",  alpha = .4) + 
  geom_density(bw = 0.04, n = 200) #+ facet_wrap(~chromosome)


ggplot(highestPurityCohortSummary, aes(purity, SUBCLONAL_SNV)) + geom_point() + facet_wrap(~cancerType)
  scale_y_log10()
