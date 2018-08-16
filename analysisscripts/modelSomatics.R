
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

find_somatic_peak <- function(somatics, peakBinwidth = 0.01, width = 5) {
  histo = add_bucket(somatics  %>% filter(weight > 0), peakBinwidth, 0) %>%
    group_by(bucket) %>%
    summarise(weight = sum(weight)) %>%
    arrange(bucket)
  
  for (i in c(1:width)) {
    histo[, paste0("lag", i)] <- lag(histo$weight, i, default = 0)
    histo[, paste0("lead", i)] <- lead(histo$weight, i, default = 0)
  }
  histo$total = rowSums(histo[, -1])
  peak_y = max(histo$total)
  peak_x = first(histo$bucket[histo$total == peak_y])
  return (as.numeric(peak_x))
}

somatics = initialSomatics %>% mutate(weight = 1)
model_peak <- function(somatics, binwidth) {

  peak = find_somatic_peak(somatics)
  peakOffset = peak - round(peak / binwidth) * binwidth
  
  peakSomatics =  somatics %>% 
    filter(somaticPloidy > peak - binwidth/2, somaticPloidy < peak + binwidth/2) %>% 
    select(somaticPloidy, alleleReadCount, totalReadCount, weight)
  peakSomaticsWeight = pmax(0, peakSomatics[["weight"]])
  
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
  
  peakDistribution = colSums(normalisePeakSomatics)
  peakDistribution = data.frame(bucket = (names(peakDistribution)), peakWeight = round(peakDistribution,2), stringsAsFactors = F)
  peakDistributionNormalised = peakDistribution %>% mutate(bucket = round(as.numeric(bucket) - peakOffset, 3))
  
  peakHeight = max(peakDistribution$peakWeight)
  peakArea =  sum(peakDistribution$peakWeight) 
  peakAvgWeight = mean(peakSomaticsWeight)

  peakSummary = data.frame(peak = peak, height = peakHeight, area = peakArea, weight = peakAvgWeight, stringsAsFactors = F)
  return (list(peak = peak, offset = peakOffset, binwidth = binwidth, summary = peakSummary, distribution = peakDistribution, normalisedDistribution = peakDistributionNormalised))
}


remove_peak_from_somatics <- function(peak, somatics) {
  
  offset = peak[["offset"]]
  binwidth = peak[["binwidth"]]
  modelDistribution = peak[["distribution"]]
  
  somatics = add_bucket(somatics, binwidth, offset) %>%
    group_by(bucket) %>%
    mutate(bucketWeight = sum(weight)) %>% 
    ungroup() %>%
    left_join(modelDistribution, by = "bucket") %>%
    ungroup()
  somatics = somatics %>%
    mutate(weight = weight - abs(peakWeight / bucketWeight)) %>%
    select(-peakWeight, -bucketWeight)
  
  return (somatics)
}

initialSomatics = enrichedProdSomatics


model_somatics <- function(initialSomatics, binwidth) {
  peaks = list()
  sample = first(initialSomatics$sampleId)

  somatics = initialSomatics %>% 
    select(somaticPloidy, alleleReadCount, totalReadCount) %>% mutate(weight = 1) %>%
    filter(somaticPloidy >= 0, somaticPloidy <= 10)
  
  peakSummary = data.frame()
  peakDistribution = data.frame()
  initialWeight = sum(somatics$weight)
  unexplained = 1
  i = 1
  
  while (unexplained > 0.05) {
    peak = model_peak(somatics, binwidth)
    somatics = remove_peak_from_somatics(peak, somatics)
    unexplained = sum(pmax(0,somatics$weight)) / initialWeight
    if (is.nan(unexplained)) {
      unexplained = 0
    }
    
    cat("Unexplained: ", unexplained, sum(pmax(0,somatics$weight)), initialWeight, "\n")
    
    summary = peak[["summary"]] %>% mutate(unexplained = unexplained, sampleId = sample) %>% select(sampleId, everything())
    peakSummary = bind_rows(peakSummary, summary)
    
    distribution = peak[["normalisedDistribution"]] %>% mutate(peak = as.character(peak[["peak"]]) , bucket = as.numeric(bucket))
    peakDistribution = bind_rows(peakDistribution, distribution)
    
    peaks[[i]] = peak
    i = i + 1
  }
  
  peakSummary = peakSummary %>% mutate(rank = row_number())
  peakDistributionTotal = peakDistribution %>% group_by(bucket) %>% summarise(peakWeight = sum(peakWeight)) %>% mutate(peak = "Total")
  peakDistribution = bind_rows(peakDistribution, peakDistributionTotal)
  peakDistribution = peakDistribution %>%  mutate(peak = factor(peak, levels = unique(peakDistribution$peak)))
  
  p = ggplot() + 
    geom_histogram(data = initialSomatics, aes(x = somaticPloidy, fill = CN), binwidth = binwidth, alpha = 0.6) +
    geom_line(data = peakDistribution, aes(x = bucket, y = peakWeight, color = peak), stat = "identity") +
    ggtitle(sample) +
    xlim(0, 5)

  save_plot(paste0("/Users/jon/hmf/analysis/somatics/", sample, ".png"), p, base_height = 6, base_width = 7)  
  
  return (peakSummary)
}

somatic_enrichment <- function(somatics) {

  somatics = somatics %>%
    filter(!chromosome %in% c('X','Y'), filter == 'PASS') %>%
    mutate(
      majorAllelePloidy = adjustedCopyNumber - minorAllelePloidy,
      somaticPloidy = adjustedVaf * adjustedCopyNumber,
      CN = round(adjustedCopyNumber),
      CN = ifelse(CN > 5, "6+", CN),
      CN = factor(CN, levels = c("6+", rev(-1:5)), ordered = T))
  
  return (somatics)
}

initialSomatics = somatic_enrichment(somatics)
model_somatics(enrichedSomatics, binwidth)



### TEST DATA
#save(sampleDetails, enrichedProdSomatics, enrichedPilotSomatics, file = "~/hmf/analysis/fit/CPCT02110040T.RData")
#load(file = "~/hmf/analysis/fit/CPCT02110040T.RData")
#model_somatics(enrichedProdSomatics, 0.1)


#### ALGO

dbPaper = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
cohort = dbGetQuery(dbPaper, "SELECT * FROM purity where qcStatus = 'PASS'")
peakSummary = data.frame()
for (sample in cohort[1:40, "sampleId"]) {
  sample = "CPCT02030224T"
  cat("Processing:", sample, "\n")
  somatics = dbGetQuery(dbPaper, paste0("SELECT * FROM somaticVariant where filter = 'PASS' and sampleId = '", sample, "'"))
  peakSummary = bind_rows(model_somatics(somatic_enrichment(somatics), binwidth), peakSummary)
}



dbDisconnect(dbPaper)
rm(dbPaper)



save(peakSummary, file = "/Users/jon/hmf/analysis/somatics/peakSummary.RData")

0/3.0



