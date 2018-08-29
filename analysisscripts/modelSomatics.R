
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

save(peak, binwidth, peakSomatics, file = "~/hmf/analysis/somatics/CPCT02010646T.RData")

#somatics = initialSomatics %>% mutate(weight = 1)
#peak = 0.99
model_peak <- function(peak, somatics, binwidth) {
  peakOffset = round(peak - round(peak / binwidth) * binwidth, 2)
  
  peakSomatics =  somatics %>% 
    filter(somaticPloidy > peak - binwidth/2, somaticPloidy < peak + binwidth/2) %>% 
    select(somaticPloidy, alleleReadCount, totalReadCount, weight)
  peakSomaticsWeight = peakSomatics$weight
  
  domain = seq(peakOffset, 10 + peakOffset, binwidth)
  for (ploidy in domain) {
    peakSomatics = peakSomatics %>% mutate(
      lowerBoundAlleleReadCount = floor(pmax(0, ploidy - binwidth/2) / somaticPloidy * alleleReadCount),
      upperBoundAlleleReadCount = floor((ploidy + binwidth/2) / somaticPloidy * alleleReadCount),
      peakSliceLikelihood = pbinom(upperBoundAlleleReadCount, totalReadCount, alleleReadCount/totalReadCount) - pbinom(lowerBoundAlleleReadCount, totalReadCount, alleleReadCount/totalReadCount)
    )
    peakSomatics[, paste(ploidy)] <-  round(peakSomatics$peakSliceLikelihood, 2)
  }
 
  peakBucket = peakSomatics[[paste(peak)]]
  peakSomatics = peakSomatics[, c(-1:-7)]
  normalisePeakSomatics = peakSomaticsWeight * peakSomatics / peakBucket 
  
  peakDistribution = col62Sums(normalisePeakSomatics)
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
model_somatics <- function(enrichedSomatics, binwidth, minAverageWeight = 0.4, minClonalCutoff = 0.9) {
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
      
    peakPloidy = estimate_somatic_peak(somatics, 0.01, 15)
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

subclonal_regions <- function(peakDistribution) {
  peakClonality = peakDistribution %>% 
    mutate(clonality = ifelse(subclonal, "subclonal", "clonal")) %>%
    group_by(sampleId, bucket, clonality) %>% 
    summarise(count = sum(count)) %>% spread(clonality, count, fill = 0) %>%
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

subclonalRegions = subclonal_regions(samplePeakDistribution)
str(subclonalRegions)
#initialSomatics = somatic_enrichment(somatics)
#model_somatics(enrichedSomatics, binwidth)



### TEST DATA
#save(sampleDetails, enrichedProdSomatics, enrichedPilotSomatics, file = "~/hmf/analysis/fit/CPCT02110040T.RData")
#load(file = "~/hmf/analysis/fit/CPCT02110040T.RData")
#model_somatics(enrichedProdSomatics, 0.1)
#sample = "CPCT02020197T"

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())


##### FROM RDATA FILES
#load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
samples = unique(allSomatics_p2$sampleId)
peakDistributions = data.frame()
binwidth = 0.05

subclonalFill = c("red", "black")
subclonalFill = setNames(subclonalFill, c("TRUE", "FALSE"))

for (sample in samples) {
  sample = "CPCT02010646T"
  cat("Processing:", sample, "\n")
  somatics = allSomatics_p2 %>% filter(sampleId == sample) %>% mutate(somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber))
  filteredSomatics = somatics %>% filter(!chromosome %in% c('X','Y'))
  
  samplePeakDistribution = model_somatics(filteredSomatics,binwidth)
  samplePeakDistributionTotal = samplePeakDistribution %>% group_by(bucket) %>% summarise(count = sum(count))
  #sum(samplePeakDistributionTotal$count)
  #nrow(initialSomatics %>% filter(somaticPloidy > 0, somaticPloidy < 10))
  
  p = ggplot() +
    geom_histogram(data=filteredSomatics, aes(x = somaticPloidy), binwidth = binwidth, alpha = 0.8) +
    geom_line(data=samplePeakDistributionTotal, aes(x = bucket, y = count), position = "identity", alpha = 0.8) +
    geom_line(data=samplePeakDistribution, aes(x = bucket, y = count, color = peak), position = "identity") +
    geom_area(data=samplePeakDistribution %>% filter(subclonal), aes(x = bucket, y = count, color = peak, fill = subclonal), position = "identity",  alpha = 0.3) +
    ggtitle(paste0(sample, ", binwidth:", binwidth)) + xlab("Ploidy") + ylab("Count") +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    xlim(0, 4)
  
p
  #subclonalRegions = subclonal_regions(samplePeakDistribution)
  #annotatedSomatics = annotate_somatics_with_clonality(somatics, subclonalRegions)
  #sum(annotatedSomatics$clonalLikelihood) + sum(annotatedSomatics$subclonalLikelihood)
  
    save_plot(paste0("/Users/jon/hmf/analysis/somatics/", sample, ".png"), p, base_height = 6, base_width = 8)  
  
  peakDistributions = bind_rows(samplePeakDistribution, peakDistributions)
}

rm(allSomatics_p1)

peakDistributions_p2 = peakDistributions
save(peakDistributions_p2, file = "/Users/jon/hmf/analysis/somatics/peakDistributions_p2.RData")


save(somatics, file = "/Users/jon/hmf/analysis/somatics/CPCT02020689T.RData")

#### FROM DATABASE
binwidth = 0.1
dbPaper = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
cohort = dbGetQuery(dbPaper, "SELECT * FROM purity where qcStatus = 'PASS'")
distributions = data.frame()
for (sample in cohort[1:1, "sampleId"]) {
  sample = "CPCT02010615T"
  cat("Processing:", sample, "\n")
  initialSomatics = (dbGetQuery(dbPaper, paste0("SELECT * FROM somaticVariant where filter = 'PASS' and sampleId = '", sample, "'")))
  enrichedSomatics = somatic_enrichment(initialSomatics)
  
  distributions = bind_rows(model_somatics(enrichedSomatics, binwidth), distributions)
}

sum(distributions$count)

dbDisconnect(dbPaper)
rm(dbPaper)



sublconality = distributions %>% group_by(sampleId, bucket, subclonal) %>% summarise(count = sum(count)) %>% spread(subclonal, count)









##### OTHER OPTIONS FOR FINDING PEAKS.. KEEP FOR LATER JUST IN CSE


#Option 0
#peakPloidy = estimate_somatic_peak(somatics, 0.01, 5)
#peakModel = model_peak(peakPloidy, somatics, binwidth)


#Option 2
#seedPloidy = estimate_somatic_peak(somatics)
#impact = data.frame()
#for (j in c(-10:10)) {
#  currentPloidy = seedPloidy + j * 0.01;
#  currentModel = model_peak(currentPloidy, somatics, binwidth)
#  currenImpact = change_in_weights(currentModel, somatics)
#  impact = bind_rows(impact, currenImpact)
#}
#peakPloidy = first(impact %>% mutate(imbalance = abs(imbalance)) %>% arrange(imbalance) %>% pull(peak))
#peakModel = model_peak(peakPloidy, somatics, binwidth)

#Option 3
#peakPloidy = estimate_somatic_peak(somatics)
#peakModel = model_peak(peakPloidy, somatics, binwidth)
#peakImpact = change_in_weights(currentModel, somatics)
#repeat {
#  if (peakImpact$right > peakImpact$left) {
#    nextPloidy = peakPloidy + 0.01
#  } else {
#    nextPloidy = peakPloidy - 0.01
#  }
#  
#  nextModel = model_peak(nextPloidy, somatics, binwidth)
#  nextImpact = change_in_weights(nextModel, somatics)
#  if (abs(nextImpact$imbalance) > abs(peakImpact$imbalance)) {
#    break;
#  } else {
#    peakPloidy = nextPloidy;
#    peakModel = nextModel;
#    peakImpact = nextImpact;
#  }
#}

change_in_weights <- function(peakModel, somatics) {
  peak = peakModel[["peak"]]
  offset = peakModel[["offset"]]
  binwidth = peakModel[["binwidth"]]
  
  modelDistribution = peakModel[["distribution"]] %>% filter(peakWeight > 0) %>% select(bucket, weightToRemove = peakWeight)
  actualDistribution = add_bucket(somatics, binwidth, offset) %>%
    group_by(bucket) %>%
    summarise(currentWeight = sum(pmax(0, weight))) %>% 
    filter(bucket %in% modelDistribution$bucket)
  
  combinedDistribtuion = left_join(modelDistribution, actualDistribution, by = "bucket")
  combinedDistribtuion[is.na(combinedDistribtuion)] <- 0
  change = combinedDistribtuion %>% 
    mutate(
      after = currentWeight - weightToRemove,
      type = ifelse(as.numeric(bucket) < peak, "left", "right")) %>% 
    group_by(type) %>% 
    summarise(value = sum(after)) %>% 
    spread(type, value) %>%
    mutate(
      peak = peak,
      imbalance = left - right
    ) %>%
    select(peak, everything())
  
  return (change)
}



totalDistribution = totalDistribution %>% group_by(bucket) %>% summarise(peakWeight = sum(peakWeight)) %>% mutate(peak = "Total")
peakDistribution = bind_rows(peakDistribution, totalDistribution)
peakDistribution = peakDistribution %>%  mutate(peak = factor(peak, levels = unique(peakDistribution$peak)))

p = ggplot() + 
  geom_histogram(data = enrichedSomatics, aes(x = somaticPloidy, fill = CN), binwidth = binwidth, alpha = 0.6) +
  geom_line(data = peakDistribution, aes(x = bucket, y = peakWeight, color = peak), stat = "identity") +
  ggtitle(paste0(sample, "", "")) +
  xlim(0, 5)

jonDistribution = totalDistribution %>% group_by(bucket, clonality) %>% summarise(peakWeight = sum(peakWeight))
clonalityProportion = jonDistribution %>% group_by(bucket) %>% spread(clonality, peakWeight) %>% 
  mutate(TOTAL = CLONAL + SUBCLONAL, CLONAL = ifelse(TOTAL == 0, 0, CLONAL / TOTAL), SUBCLONAL = ifelse(TOTAL == 0, 0, SUBCLONAL / TOTAL)) %>% select(-TOTAL) %>%
  gather(clonality, value, CLONAL, SUBCLONAL)

ggplot()+
  geom_line(data = jonDistribution, aes(x = bucket, y = peakWeight, color = clonality), position = "stack")

ggplot()+
  geom_bar(data = jonDistribution, aes(x = bucket, y = peakWeight, fill = clonality), position = "stack", stat = "identity") + xlim(0, 5) 

ggplot() +
  geom_bar(data = clonalityProportion, width = 0.1,  aes(x = bucket, y = value, fill = clonality), stat = "identity") + xlim(0, 5) 


dev.off()
p

save_plot(paste0("/Users/jon/hmf/analysis/somatics/", sample, ".png"), p, base_height = 6, base_width = 7)  
