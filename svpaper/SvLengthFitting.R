library(RMySQL)
library(purple)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

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


allSvs = read.csv(file = "/Users/jon/hmf/analysis/svPaper/SVA_SVS.csv")
load(file = "/Users/jon/hmf/analysis/svPaper/cohort.RData")
hpcSvs = allSvs %>% filter(SampleId %in% highestPurityCohort$sampleId)

#rm(list=setdiff(ls(), c("hpcSvs", "peaks")))

hpcDels = hpcSvs %>% filter(ResolvedType == 'SimpleSV', Type == 'DEL') %>% 
  mutate(length = PosEnd - PosStart + 1,  log2Length = log2(length)) %>% 
  group_by(SampleId)

samples = hpcDels %>% group_by(SampleId) %>% count() %>% arrange(-n)
  
result = data.frame(stringsAsFactors = F)
for (sample in samples[1:100, ]$SampleId) {
  cat(sample, "\n")
  
  sampleDels = hpcDels %>% filter(SampleId == sample)
  sampleDensity = density(sampleDels$log2Length, bw = 0.6, n = 200)
  samplePeaks = peaks(sample, sampleDensity$x, sampleDensity$y)
  samplePlot = ggplot(data = sampleDels, aes(log2Length)) + 
    geom_histogram(aes(y =..density..), bins = 200, fill="blue", col="blue",  alpha = .4) + 
    xlim(1, 25) +
    geom_density(bw = 0.6, n = 200)
  
  save_plot(paste0("/Users/jon/hmf/analysis/svPaper/svLength/", sample, ".del.png"), samplePlot, base_width = 3, base_height = 3)
  
  result = bind_rows(result, samplePeaks)
}


delPeakSummary = result %>% 
  filter(peakHeight >= 0.01) %>%
  mutate(rankHeight = order) %>%
  select(sampleId, order, rankHeight, peakLength = peakPloidy, peakHeight) %>%
  group_by(sampleId) %>% mutate(order = paste0("Peak", order), rankHeight = paste0("Peak", rankHeight, "Height")) %>% spread(order, peakLength, fill = 0) %>% 
  spread(rankHeight, peakHeight, fill = 0) %>%
  group_by(sampleId) %>%
  summarise(
    Peak1 = sum(Peak1), Peak2 = sum(Peak2), Peak3 = sum(Peak3), Peak4 = sum(Peak4),
    Peak1Height = sum(Peak1Height), Peak2Height = sum(Peak2Height), Peak3Height = sum(Peak3Height), Peak4Height = sum(Peak4Height) 
    )


delPeak = result
save(delPeak, delPeakSummary, file = "/Users/jon/hmf/analysis/svPaper/svLength/delPeaks.RData")

#### INDIVIDUAL SAMPLE
sample = "CPCT02010790T"
sample = "CPCT02050135T"
sampleDels = hpcDels %>% filter(SampleId == sample)
sampleDensity = density(sampleDels$log2Length, bw = 0.6, n = 50)
peaks(sample, sampleDensity$x, sampleDensity$y)
ggplot(data = sampleDels, aes(log2Length)) + 
  geom_histogram(aes(y =..density..), bins = 50, fill="blue", col="blue",  alpha = .4) + 
  xlim(1, 25) +
  geom_density(bw = 0.6, n = 50) #+ facet_wrap(~ChrStart)



### PLOT ALL 
ggplot(data = dels, aes(length1)) + 
  geom_histogram(aes(y =..density..), bins = 200, fill="blue", col="blue",  alpha = .4) + 
  xlim(1, 25) +
  geom_density(bw = 0.3, n = 200) #+ facet_wrap(~ChrStart)



