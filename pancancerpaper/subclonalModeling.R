
### FUNCTIONS
add_bucket <- function(somatics, binwidth, offset) {
  domain = seq(0 - binwidth / 2 + offset , 10 + binwidth / 2 + offset, binwidth)
  somatics = somatics %>%
    mutate(
      bucketFactor = cut(somaticPloidy, breaks = domain),
      bucket = as.character((as.numeric(bucketFactor) - 1) * binwidth + offset)
    ) %>%
    select(-bucketFactor)
  return (somatics)
}

#peakBinwidth = 0.01;width = 10
estimate_somatic_peak <- function(somatics, peakBinwidth, width) {
  histo = add_bucket(somatics, peakBinwidth, 0) %>%
    group_by(bucket) %>%
    summarise(weight = sum(weight)) %>%
    arrange(bucket)
  
  for (i in c(1:width)) {
    histo[, paste0("lag", i)] <- lag(histo$weight, i, default = 0)
    histo[, paste0("lead", i)] <- lead(histo$weight, i, default = 0)
  }
  histo$total = rowSums(histo[, -1])
  histo = histo %>% filter(weight > 0)
  
  peak_y = max(histo$total)
  peak_x = first(histo$bucket[histo$total == peak_y])
  return (as.numeric(peak_x))
}

model_peak <- function(peak, somatics, binwidth) {
  peakOffset = round(peak - round(peak / binwidth) * binwidth, 2)
  
  peakSomatics =  somatics %>% 
    filter(somaticPloidy > peak - binwidth/2, somaticPloidy < peak + binwidth/2) %>% 
    select(somaticPloidy, alleleReadCount, totalReadCount, weight)
  peakSomaticsWeight = peakSomatics$weight
  
  domain = seq(peakOffset, 10 + peakOffset, binwidth)
  for (ploidy in domain) {
    peakSomatics = peakSomatics %>% mutate(
      lowerBoundAlleleReadCount = pmax(0, ploidy - binwidth/2.0) / somaticPloidy * alleleReadCount,
      lowerBoundAlleleReadCountRounded = round(lowerBoundAlleleReadCount, 0),
      lowerBoundAddition = lowerBoundAlleleReadCountRounded + 0.5 - lowerBoundAlleleReadCount,
      
      upperBoundAlleleReadCount = pmax(0, ploidy + binwidth/2.0) / somaticPloidy * alleleReadCount,
      upperBoundAlleleReadCountRounded = round(upperBoundAlleleReadCount, 0),
      upperBoundSubtraction = upperBoundAlleleReadCountRounded + 0.5 - upperBoundAlleleReadCount,
      
      peakSliceLikelihood = 
        pbinom(upperBoundAlleleReadCountRounded, totalReadCount, alleleReadCount/totalReadCount) - 
        pbinom(lowerBoundAlleleReadCountRounded, totalReadCount, alleleReadCount/totalReadCount) +
        lowerBoundAddition * dbinom(lowerBoundAlleleReadCountRounded, totalReadCount, alleleReadCount/totalReadCount) -
        upperBoundSubtraction * dbinom(upperBoundAlleleReadCountRounded, totalReadCount, alleleReadCount/totalReadCount)
    )
    peakSomatics[, paste(ploidy)] <-  round(peakSomatics$peakSliceLikelihood, 2)
  }

  peakBucket = peakSomatics[[paste(peak)]]
  peakSomatics = peakSomatics[, c(-1:-11)]
  normalisePeakSomatics = peakSomaticsWeight * peakSomatics / peakBucket 
  
  peakDistribution = colSums(normalisePeakSomatics)
  peakDistribution = data.frame(bucket = (names(peakDistribution)), count = round(peakDistribution,2), stringsAsFactors = F)
  peakDistributionNormalised = peakDistribution %>% mutate(bucket = round(as.numeric(bucket) - peakOffset, 3))
  
  #ggplot() +
  #geom_histogram(data = peakSomatics, aes(x = somaticPloidy), binwidth = 0.1) +
  #geom_line(data = peakDistributionNormalised, aes(x = bucket, y = count, group = 1)) + xlim(0, 5)
  
  peakHeight = max(peakDistribution$count)
  peakArea =  sum(peakDistribution$count) 
  peakAvgWeight = mean(peakSomaticsWeight)

  peakSummary = data.frame(peak = peak, height = peakHeight, area = peakArea, avgWeight = peakAvgWeight, stringsAsFactors = F)
  return (list(peak = peak, offset = peakOffset, binwidth = binwidth, summary = peakSummary, distribution = peakDistribution, normalisedDistribution = peakDistributionNormalised))
}


remove_peak_from_somatics <- function(peakModel, somatics) {
  
  offset = peakModel[["offset"]]
  binwidth = peakModel[["binwidth"]]
  modelDistribution = peakModel[["distribution"]]
  
  somaticsWithDistribution = add_bucket(somatics, binwidth, offset) %>%
    group_by(bucket) %>%
    mutate(bucketWeight = sum(weight)) %>% 
    ungroup() %>%
    left_join(modelDistribution, by = "bucket") %>%
    ungroup()
  somaticsWithDistribution = somaticsWithDistribution %>%
    mutate(weight = ifelse(bucketWeight == 0, 0, weight - abs(count / bucketWeight))) %>%
    select(-count, -bucketWeight)
  
  return (somaticsWithDistribution)
}

#minAverageWeight = 0.15; minClonalCutoff = 0.9;
#enrichedSomatics = filteredSomatics
model_somatics <- function(enrichedSomatics, binwidth, minAverageWeight = 0.4, minClonalCutoff = 0.85, estimateWidth) {
  peaks = list()
  sample = first(enrichedSomatics$sampleId)

  somatics = enrichedSomatics %>% 
    filter(
      somaticPloidy > 0, 
      somaticPloidy < 10, 
      !chromosome %in% c('X','Y')) %>%
    select(somaticPloidy, alleleReadCount, totalReadCount) %>% 
    mutate(weight = 1) 

  initialWeight = sum(somatics$weight)
  unexplained = 1
  i = 1
  peakDistribution = data.frame()
  
  while (unexplained > 0.01 & i < 10 & any(somatics$weight > 0)) {
  #for (i in c(1:3)) {
      
    peakPloidy = estimate_somatic_peak(somatics, 0.01, estimateWidth)
    peakModel = model_peak(peakPloidy, somatics, binwidth)

    somatics = remove_peak_from_somatics(peakModel, somatics)
    unexplained = sum(pmax(0,somatics$weight)) / initialWeight
    cat("Unexplained: ", unexplained, sum(pmax(0,somatics$weight)), initialWeight, "\n")

    summary = peakModel[["summary"]] %>% 
      mutate(
        include = avgWeight > minAverageWeight,
        subclonal = peak < minClonalCutoff,
        sampleId = sample) %>% 
      select(sampleId, everything())

    jon = peakModel[["distribution"]] 
    jon2 = peakModel[["normalisedDistribution"]] 
    
    if (summary$include) {
      normalisedDistribution = peakModel[["normalisedDistribution"]] %>% mutate(sampleId = sample, peak = as.character(peakModel[["peak"]]), bucket = as.numeric(bucket), subclonal = summary$subclonal)
      peakDistribution = bind_rows(peakDistribution, normalisedDistribution)
    } 

    i = i + 1
  }
  
  totalCount = sum(peakDistribution$count)
  peakDistribution$count = peakDistribution$count / totalCount * initialWeight
  
  return (peakDistribution)
}

somatic_enrichment <- function(somatics) {

  somatics = somatics %>%
    filter(filter == 'PASS') %>%
    mutate(
      majorAllelePloidy = adjustedCopyNumber - minorAllelePloidy,
      somaticPloidy = adjustedVaf * adjustedCopyNumber,
      CN = round(adjustedCopyNumber),
      CN = ifelse(CN > 5, "6+", CN),
      CN = factor(CN, levels = c("6+", rev(-1:5)), ordered = T))
  
  return (somatics)
}

subclonal_regions <- function(samplePeakDistribution, subclonalCutoff = 0.85) {
  peakClonality = samplePeakDistribution %>% 
    mutate(clonality = ifelse(peak < subclonalCutoff, "subclonal", "clonal")) %>%
    group_by(sampleId, bucket, clonality) %>% 
    summarise(count = sum(count)) %>% spread(clonality, count, fill = 0)
  
  if (!"subclonal" %in% names(peakClonality)) {
    peakClonality$subclonal <- 0
  }
  
  if (!"clonal" %in% names(peakClonality)) {
    peakClonality$clonal <- 0
  }
  
  peakClonality = peakClonality %>%
    mutate(clonal = round(clonal, 2), subclonal = round(subclonal, 2)) %>%
    mutate(
      total = clonal + subclonal, 
      subclonal = ifelse(total == 0, 0, round(subclonal / total, 2)), 
      clonal =  1 - subclonal) %>%
    select(-total) %>%
    filter(subclonal > 0)  %>%
    ungroup() %>%
    mutate(bucket = as.character(bucket)) %>%
    select(sampleId, bucket, subclonalLikelihood = subclonal)
  
  return (peakClonality)
}

annotate_somatics_with_clonality <- function(enrichedSomatics, subclonalRegions) {
  somaticsWithBuckets = add_bucket(enrichedSomatics, binwidth, 0) %>%
    left_join(subclonalRegions, by = c("sampleId", "bucket"))
  somaticsWithBuckets$subclonalLikelihood <- ifelse(is.na(somaticsWithBuckets$subclonalLikelihood), 0, somaticsWithBuckets$subclonalLikelihood)
  somaticsWithBuckets$clonalLikelihood <- 1 - somaticsWithBuckets$subclonalLikelihood
  return (somaticsWithBuckets)
}





#subclonalRegions = subclonal_regions(samplePeakDistribution)
#str(subclonalRegions)
#initialSomatics = somatic_enrichment(somatics)
#model_somatics(enrichedSomatics, binwidth)


### TEST DATA
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())


##### FROM RDATA FILES
#load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
#load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
samples = unique(allSomatics_p2$sampleId)
peakDistributions = data.frame()
modifiedSomatics = data.frame()

binwidth = 0.05
estimateWidth = 10
i = 0
for (sample in samples) {

  i = i + 1
  cat("Processing:", sample, " (",i,")", "\n")
  somatics = allSomatics_p2 %>% filter(sampleId == sample) %>% mutate(somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber)) %>% 
    select(-clonalLikelihood, -subclonalLikelihood)
  filteredSomatics = somatics %>% filter(!chromosome %in% c('X','Y'))
  
  samplePeakDistribution = model_somatics(filteredSomatics,binwidth, estimateWidth = estimateWidth)
  peakDistributions = bind_rows(samplePeakDistribution, peakDistributions)
  
  #samplePeakDistributionTotal = samplePeakDistribution %>% group_by(bucket) %>% summarise(count = sum(count))
  #p = ggplot() +
  #  geom_histogram(data=filteredSomatics, aes(x = somaticPloidy), binwidth = binwidth, alpha = 0.8) +
  #  geom_line(data=samplePeakDistributionTotal, aes(x = bucket, y = count), position = "identity", alpha = 0.8) +
  #  geom_line(data=samplePeakDistribution, aes(x = bucket, y = count, color = peak), position = "identity") +
  #  geom_area(data=samplePeakDistribution %>% filter(subclonal), aes(x = bucket, y = count, color = peak, fill = subclonal), position = "identity",  alpha = 0.3) +
  #  ggtitle(paste0(sample, ", binwidth:", binwidth, ", estimate width: ", estimateWidth)) + xlab("Ploidy") + ylab("Count") +
  #  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
  #  xlim(0, 4)
  #save_plot(paste0("/Users/jon/hmf/analysis/somatics/", sample, ".png"), p, base_height = 6, base_width = 8)  
  
  subclonalRegions = subclonal_regions(samplePeakDistribution)
  annotatedSomatics = annotate_somatics_with_clonality(somatics, subclonalRegions) %>% select(-bucket, -somaticPloidy) 
  modifiedSomatics = bind_rows(annotatedSomatics, modifiedSomatics)
}

rm(allSomatics_p2)
allSomatics_p2 = modifiedSomatics
save(allSomatics_p2, file = "/Users/jon/hmf/analysis/somatics/RData/allSomatics_p2.RData")

peakDistributions_p2 = peakDistributions
save(peakDistributions_p2, file = "/Users/jon/hmf/analysis/somatics/RData/peakDistributions_p2.RData")

head(allSomatics_p1)

### FIX UP ODL
rm(list=setdiff(ls(), "allSomatics_p1"))
load(file = "/Users/jon/hmf/analysis/somatics/RData/modifiedSomatics_p2.RData")
allSomatics_p2 = modifiedSomatics
rm(modifiedSomatics)
save(allSomatics_p2, file = "/Users/jon/hmf/analysis/somatics/RData/allSomatics_p2.RData")


##### APPLY PEAK_DISTRIBUTIONS
library(multidplyr)
rm(list=setdiff(ls(), c("allSomatics_p2", "peakDistributions_p2") ))


allSomatics_p2 = allSomatics_p2 %>% select(-subclonalLikelihood, -clonalLikelihood)
head(allSomatics_p2)
head(peakDistributions_p2)

subclonalRegions_p2 = subclonal_regions(peakDistributions_p2, 0.85)
save(subclonalRegions_p2, file = "/Users/jon/hmf/analysis/somatics/RData/subclonalRegions_p2.RData")


allSomatics_p2 = allSomatics_p2 %>% mutate(somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber))
date()
allSomatics_p2 = annotate_somatics_with_clonality(allSomatics_p2, subclonalRegions_p2)
allSomatics_p2 = allSomatics_p2 %>% select(-bucket, -somaticPloidy) 
date()
save(allSomatics_p2, file = "/Users/jon/hmf/RData/Reference/allSomatics_p2.RData")



load("/Users/jon/hmf/analysis/somatics/RData/peakDistributions_p1.RData")
subclonalRegions_p1 = subclonal_regions(peakDistributions_p1, 0.85)
save(subclonalRegions_p1, file = "/Users/jon/hmf/analysis/somatics/RData/subclonalRegions_p1.RData")


load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
load(file = "~/hmf/analysis/somatics/RData/subclonalRegions_p1.RData")

allSomatics_p1 = allSomatics_p1  %>% 
  select(-subclonalLikelihood, -clonalLikelihood) %>% 
  mutate(somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber))

allSomatics_p1 = annotate_somatics_with_clonality(allSomatics_p1, subclonalRegions_p1) %>% 
  select(-bucket, -somaticPloidy) 

save(allSomatics_p1, file = "/Users/jon/hmf/RData/Reference/allSomatics_p1.RData")


################# Methods Figure 1 - Subclonal.png 
library(cowplot)
library(scales)
theme_set(theme_bw())

singleBlue = "#6baed6"
singleRed = "#d94701"

binwidth = 0.05
estimateWidth = 10
query = "select * from somaticVariant where sampleId = 'xxxxx'"
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
somatics = dbGetQuery(dbProd, query) %>% mutate(somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber)) 
dbDisconnect(dbProd)
rm(dbProd)

filteredSomatics = somatics %>% filter(!chromosome %in% c('X','Y'))
samplePeakDistribution = model_somatics(filteredSomatics,binwidth, estimateWidth = estimateWidth)
samplePeakDistribution2 = samplePeakDistribution %>% filter(!peak %in% c("2.7","5.97") )

samplePeakDistributionTotal = samplePeakDistribution2 %>% group_by(bucket) %>% summarise(count = sum(count))

pDistribution = ggplot() +
  geom_histogram(data=filteredSomatics, aes(x = somaticPloidy), binwidth = binwidth, fill=singleBlue, col=singleBlue,  alpha = .4) +
  geom_line(data=samplePeakDistributionTotal, aes(x = bucket, y = count), position = "identity", alpha = 0.8) +
  geom_line(data=samplePeakDistribution2, aes(x = bucket, y = count, color = peak), position = "identity") +
  geom_area(data=samplePeakDistribution2 %>% filter(subclonal), aes(x = bucket, y = count, fill = subclonal), position = "identity",  alpha = 0.3, color = singleRed) +
  geom_segment(aes(x = 0.85, xend = 0.85, y = 0, yend = 1900), linetype = "dashed") +
  geom_text(aes(label = "Subclonal Peaks", x = 0.8, y = 1900), size = 3.5, hjust = 1) +
  geom_text(aes(label = "Clonal Peaks", x = 0.9, y = 1900), size = 3.5, hjust = 0) +
  ggtitle("") + xlab("Ploidy") + ylab("") +
  scale_y_continuous(expand=c(0.02, 0.02)) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position="none") +
  scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 


jon = samplePeakDistribution %>% group_by(bucket, subclonal) %>% summarise(count = sum(count))%>% group_by(bucket) %>% mutate(percent = count / sum(count))

pLikelihood = ggplot(data = jon %>% filter(subclonal)) +
  geom_bar(width = 0.05, aes(x = bucket, y = percent), stat = "identity", fill=singleRed, col=singleRed,  alpha = 0.3) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
  xlab("") + ylab("") +
  scale_y_continuous(labels = percent, expand=c(0.02, 0.02), limits = c(0, 1)) +
  scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 

pSubclonal = plot_grid(pDistribution, pLikelihood, nrow = 2, rel_heights = c(5, 1), labels = "AUTO")
pSubclonal

save_plot("~/hmf/RPlot/Methods Figure 1 - Subclonal.png", pSubclonal, base_width = 9, base_height = 7)  





