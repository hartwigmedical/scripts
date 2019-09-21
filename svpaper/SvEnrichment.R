library(RMySQL)
library(purple)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)

########################## BINS
library(BSgenome.Hsapiens.UCSC.hg19)
bins_100k <- tileGenome(seqinfo(Hsapiens), tilewidth=100000, cut.last.tile.in.chrom=TRUE)
bins_1M <- tileGenome(seqinfo(Hsapiens), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
bins_10M <- tileGenome(seqinfo(Hsapiens), tilewidth=10000000, cut.last.tile.in.chrom=TRUE)

#getCentromeres(Hsapiens) TODO

########################## Average Mappability
mappability <- function(bins, genome) {
  df = data.frame(chromosome = seqnames(bins), start = start(bins), end = end(bins)) %>%
    mutate(chromosome = as.character(chromosome), row = row_number()) %>%
    filter(nchar(chromosome) <= 5, chromosome != 'chrM') %>%
    group_by(row) %>%
    mutate(N = countPattern("N", getSeq(genome, chromosome, start, end))) %>%
    ungroup() %>%
    mutate(chromosome = gsub("chr", "", chromosome)) %>%
    select(-row)
  return (df)
}

mappability_100k = mappability(bins_100k, BSgenome.Hsapiens.UCSC.hg19)
mappability_1M = mappability(bins_1M, BSgenome.Hsapiens.UCSC.hg19)
mappability_10M = mappability(bins_10M, BSgenome.Hsapiens.UCSC.hg19)

save(mappability_100k, mappability_1M, mappability_10M, file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")


########################## Average Copy Number 
averageCopyNumber <- function(bins, hpcCopyNumbers) {
  hpcCopyNumberRegions = GRanges(paste0("chr", hpcCopyNumbers$chromosome), ranges = IRanges(start = hpcCopyNumbers$start, end = hpcCopyNumbers$end))
  ol = as.matrix(findOverlaps(bins, hpcCopyNumberRegions, type = "any"))
  olCopyNumbers = cbind(bin = ol[, 1], binChr = seqnames(bins)[ol[, 1]], binStart = start(bins)[ol[, 1]], binEnd = end(bins)[ol[, 1]],  hpcCopyNumbers[ol[, 2], ])

  averageCopyNumbers = olCopyNumbers %>% group_by(bin, binChr, binStart, binEnd) %>%
    mutate(averageStart = pmax(start, binStart), averageEnd = pmin(end, binEnd), weight = (averageEnd - averageStart) / 1000000,  weightedCopyNumber = copyNumber * weight) %>%
    summarise(averageCopyNumber = sum(copyNumber * weight) / sum(weight))

  averageCopyNumbers = averageCopyNumbers %>% ungroup() %>% mutate(chromosome = gsub("chr", "", binChr)) %>% select(chromosome, start = binStart, end = binEnd, averageCopyNumber)
  return (averageCopyNumbers)
}

load(file = "/Users/jon/hmf/analysis/cohort/hpcCopyNumbers.RData")
averageCopyNumber_100k = averageCopyNumber(bins_100k, hpcCopyNumbers)
averageCopyNumber_1M = averageCopyNumber(bins_1M, hpcCopyNumbers)
averageCopyNumber_10M = averageCopyNumber(bins_10M, hpcCopyNumbers)
save(averageCopyNumber_100k, averageCopyNumber_1M, averageCopyNumber_10M, file = "/Users/jon/hmf/analysis/svPaper/averageCopyNumbers.RData")

########################## Average Replication
replication = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>%
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)

averageReplication <- function(bins, replication) {
  result = averageCopyNumber(bins, replication %>% select(chromosome, start, end, copyNumber = replication))
  result = result %>%
    mutate(averageReplication = averageCopyNumber / 100) %>%
    select(chromosome, start, end, averageReplication)
  return (result)
}


averageReplication_100k = averageReplication(bins_100k, replication)
averageReplication_1M = averageReplication(bins_1M, replication)
averageReplication_10M = averageReplication(bins_10M, replication)
save(averageReplication_100k, averageReplication_1M, averageReplication_10M, file = "/Users/jon/hmf/analysis/svPaper/averageReplication.RData")

########################## Average GC
gc = read.table("~/hmf/resources/GC_profile.1000bp.cnp", sep = '\t', header = F, stringsAsFactors = F) %>%
  mutate(chromosome = V1, start = V2 + 1, end = start + 1000, gc = V3) %>%
  select(chromosome, start, end, gc) %>%
  filter(gc > -1)

averageGC <- function(bins, gc) {
  result = averageCopyNumber(bins, gc %>% select(chromosome, start, end, copyNumber = gc))
  result = result %>%
    mutate(averageGC = averageCopyNumber) %>%
    select(chromosome, start, end, averageGC)
  return (result)
}

averageGC_100k = averageGC(bins_100k, gc)
averageGC_1M = averageGC(bins_1M, gc)
averageGC_10M = averageGC(bins_10M, gc)
save(averageGC_100k, averageGC_1M, averageGC_10M, file = "/Users/jon/hmf/analysis/svPaper/averageGC.RData")


########################## Uniquely Mappable
mappability_150 = read.table(file = "/Users/jon/hmf/resources/unique.mappability.bed", sep = "\t")
colnames(mappability_150) <- c("chromosome", "start0", "end", "id", "mappability")
mappability_150 = mappability_150 %>%
  filter(mappability >= 1.00, chromosome != 'MT') %>%
  mutate(start = start0 + 1) %>%
  select(chromosome, start, end, mappability)

uniquelyMappable <- function(bins, mappability) {
  mappableRegions = GRanges(paste0("chr",mappability$chromosome), ranges = IRanges(start = mappability$start, end = mappability$end))
  ol = as.matrix(findOverlaps(mappableRegions, bins, type = "any"))

  df = mappability[ol[, 1], ]
  df$binStart <- start(bins)[ol[, 2]]
  df$binEnd <- end(bins)[ol[, 2]]
  df$actualStart = pmax(df$start, df$binStart)
  df$actualEnd = pmin(df$end, df$binEnd)
  df$size = df$actualEnd - df$actualStart + 1
  sum(df$size)

  result = df %>%
    group_by(chromosome, binStart, binEnd) %>%
    summarise(uniquelyMappableBases = sum(size)) %>%
    mutate(uniquelyMappablePercentage = uniquelyMappableBases / (binEnd - binStart + 1))
}

uniquelyMappable_100k = uniquelyMappable(bins_100k, mappability_150)
uniquelyMappable_1M = uniquelyMappable(bins_1M, mappability_150)
uniquelyMappable_10M = uniquelyMappable(bins_10M, mappability_150)
save(uniquelyMappable_100k, uniquelyMappable_1M, uniquelyMappable_10M, file = "/Users/jon/hmf/analysis/svPaper/uniquelyMappable.RData")


############################################ Functions ############################################ 
replication_bucket <- function(replication) {
  cut(replication, breaks = c(-0.1, seq(0.062, 0.78, 0.02), 1.1))
}

gc_bucket <- function(x) {
  cut(x, breaks = c(-0.1, seq(0.26, 0.64, 0.02), 1.1))
}


enrich <- function(criteria, mappability, averageCopyNumber, uniquelyMappable, averageReplication, averageGC) {
  if (nrow(criteria) == 0) {
    return (criteria)
  }
  
  svsRegions = GRanges(criteria$Chr, ranges = IRanges(start = criteria$Pos, end = criteria$Pos))
  binRegions = GRanges(mappability$chromosome, ranges = IRanges(start = mappability$start, end = mappability$end))
  uniquelyMappableRegions = GRanges(as.character(uniquelyMappable$chromosome), ranges = IRanges(start = uniquelyMappable$binStart, end = uniquelyMappable$binEnd))
  replicationRegions = GRanges(as.character(averageReplication$chromosome), ranges = IRanges(start = averageReplication$start, end = averageReplication$end))
  gcRegions = GRanges(as.character(averageGC$chromosome), ranges = IRanges(start = averageGC$start, end = averageGC$end))
  
  
  ol = as.matrix(findOverlaps(svsRegions, binRegions, type = "within"))
  ol2 = as.matrix(findOverlaps(svsRegions, uniquelyMappableRegions, type = "within"))
  ol3 = as.matrix(findOverlaps(svsRegions, replicationRegions, type = "within"))
  ol4 = as.matrix(findOverlaps(svsRegions, gcRegions, type = "within"))
    
  criteria[ol[, 1], "N"] <- mappability[ol[, 2], ]$N
  criteria[ol[, 1], "averageCopyNumber"] <- averageCopyNumber[ol[, 2], ]$averageCopyNumber
  criteria[ol[, 1], "binChromosome"] <- mappability[ol[, 2], ]$chromosome
  criteria[ol[, 1], "binStart"] <- mappability[ol[, 2], ]$start
  criteria[ol[, 1], "binEnd"] <- mappability[ol[, 2], ]$end
  
  criteria[ol2[, 1], "uniquelyMappablePercentage"] <- uniquelyMappable[ol2[, 2], ]$uniquelyMappablePercentage
  criteria$uniquelyMappablePercentage <- ifelse(is.na(criteria$uniquelyMappablePercentage), 0, criteria$uniquelyMappablePercentage)
  
  criteria[ol3[, 1], "averageReplication"] <- averageReplication[ol3[, 2], ]$averageReplication
  criteria$averageReplication <- ifelse(is.na(criteria$averageReplication), 0.00001, criteria$averageReplication)
  criteria$averageReplication <- ifelse(criteria$averageReplication < 0.00001, 0.00001, criteria$averageReplication)
  criteria$replicationBucket = replication_bucket(criteria$averageReplication)
  
  criteria[ol4[, 1], "averageGC"] <- averageGC[ol4[, 2], ]$averageGC
  criteria$averageGC <- ifelse(is.na(criteria$averageGC), 0.00001, criteria$averageGC)
  criteria$averageGC <- ifelse(criteria$averageGC < 0.00001, 0.00001, criteria$averageGC)
  criteria$gcBucket = gc_bucket(criteria$averageGC)
  
  criteria = criteria %>% 
    left_join(replicationFactor, by = c("feature", "replicationBucket")) %>%
    mutate(replicationFactor = ifelse(is.na(replicationFactor), 1, replicationFactor)) %>% 
    left_join(gcFactor, by = c("feature", "gcBucket")) %>%
    mutate(gcFactor = ifelse(is.na(gcFactor), 1, gcFactor))
  
  return (criteria)
}


normalise <- function(svs, buckets) {
  normalised_svs = svs %>%
    group_by(feature) %>%
    mutate(featureRows=n()) %>%
    group_by(binChromosome, binStart, binEnd, N, featureRows, averageCopyNumber, uniquelyMappablePercentage, averageReplication, feature, gcFactor, replicationFactor) %>%
    summarise(
      unnormalisedBucketCount = n()) %>%
    group_by(feature) %>%
    mutate(
      unnormalisedBucketCount = unnormalisedBucketCount,
      #expectedBucketCount = nrows / buckets * replicationFactor * averageCopyNumber / meanCopyNumber,
      expectedBucketCount = featureRows * replicationFactor * averageCopyNumber * gcFactor /sum(replicationFactor * averageCopyNumber * gcFactor),
      p=ppois(unnormalisedBucketCount, expectedBucketCount, FALSE),
      q=p.adjust(p,"BH", buckets),
      buckets = buckets)
  
    normalised_svs %>% group_by(feature) %>% count()
    
  return (normalised_svs)
}

enrich_and_normalise <- function(criteria, qFilter = 0.01) {
  
  criteria = featuredBreakends
  
  enriched_100k = enrich(criteria, mappability_100k, averageCopyNumber_100k, uniquelyMappable_100k, averageReplication_100k, averageGC_100k)
  enriched_1M = enrich(criteria, mappability_1M, averageCopyNumber_1M, uniquelyMappable_1M, averageReplication_1M, averageGC_1M)
  enriched_10M = enrich(criteria, mappability_10M, averageCopyNumber_10M, uniquelyMappable_10M, averageReplication_10M, averageGC_10M)
  
  normalised_100k = normalise(enriched_100k, bins_100k) #%>% filter(q < qFilter)
  normalised_1M = normalise(enriched_1M, bins_1M) # %>% filter(q < qFilter)
  normalised_10M = normalise(enriched_10M, bins_10M) #%>% filter(q < qFilter)
  
  criteriaResult = data.frame()
  criteriaResult = bind_rows(criteriaResult, normalised_100k)
  criteriaResult = bind_rows(criteriaResult, normalised_1M)
  criteriaResult = bind_rows(criteriaResult, normalised_10M)
  
  return (criteriaResult)
}



############################################ BINS ############################################ 
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
allowedBins = mappability_100k %>% mutate(mappable = (100000 - N) / 100000) %>% filter(mappable >= 0.8)
binRegions = GRanges(allowedBins$chromosome, ranges = IRanges(start = allowedBins$start, end = allowedBins$end))
regions_1M = GRanges(mappability_1M$chromosome, ranges = IRanges(start = mappability_1M$start, end = mappability_1M$end))
regions_10M = GRanges(mappability_10M$chromosome, ranges = IRanges(start = mappability_10M$start, end = mappability_10M$end))

bins_100k = nrow(allowedBins)
bins_1M = length(unique(as.matrix(findOverlaps(regions_1M, binRegions, type = "any"))[, 1]))
bins_10M = length(unique(as.matrix(findOverlaps(regions_10M, binRegions, type = "any"))[, 1]))
save(bins_100k, bins_1M, bins_10M, file = "/Users/jon/hmf/analysis/svPaper/binSizes.RData")

############################################ RUN FROM HERE ############################################ 
load(file = "/Users/jon/hmf/analysis/svPaper/binSizes.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/uniquelyMappable.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/averageReplication.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/averageCopyNumbers.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/averageGC.RData")
load(file = "~/hmf/analysis/svPaper/replicationFactor.RData")
load(file = "~/hmf/analysis/svPaper/gcFactor.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/featuredBreakends.RData")
featuredBreakends = featuredBreakends %>% mutate(IsFS = F)

buckets = enrich_and_normalise(featuredBreakends, 100000) 
save(featuredBreakends, buckets, replicationFactor, file = "/Users/jon/hmf/analysis/svPaper/buckets.RData")



############################################ Plot ############################################ 
plot_enriched_locations <- function(dataf, selectedBucket, subtypes, file) {
  pdf(file=file,width=10, height = 6)
  for (selectedSubtype in subtypes) {
    df = dataf %>%
      filter(buckets == selectedBucket) %>%
      filter(feature == selectedSubtype) %>%
      mutate(
        significant = q < 0.01,
        size = ifelse(significant, 2, 0.3),
        binChromosome = factor(binChromosome, levels = c(1:22, 'X', 'Y'), ordered = T))
    
    myPlot = ggplot(df) +
      geom_point(aes(x = binStart, y = unnormalisedBucketCount, color = significant, size = significant), stat = "identity", position = "stack") +
      facet_wrap(~binChromosome) + ggtitle(unique(df$feature)) +
      scale_y_log10() +
      scale_size_manual(values= c(0.2, 0.7)) +
      xlab("location") + ylab("un-normalised bucket count")
    
    print(myPlot)
    
  }
  dev.off()
}

buckets[is.na(buckets)] <- 1

uniqueResolvedTypes = as.character(unique(featuredBreakends$feature))
plot_enriched_locations(buckets, 313, uniqueResolvedTypes, "~/hmf/analysis/svPaper/plot/ResolvedType10M.pdf")
plot_enriched_locations(buckets, 2915, uniqueResolvedTypes, "~/hmf/analysis/svPaper/plot/ResolvedType1M.pdf")
plot_enriched_locations(buckets, 28475, uniqueResolvedTypes, "~/hmf/analysis/svPaper/plot/ResolvedType100k.pdf")
