library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

singleBlue = "#6baed6"
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")

filter_unmappable <- function(x, mappability) {
  allowedBins = mappability %>% mutate(mappable = (100000 - N) / 100000) %>% filter(mappable >= 0.8)
  binRegions = GRanges(allowedBins$chromosome, ranges = IRanges(start = allowedBins$start, end = allowedBins$end)) 
  xRegions = GRanges(x$chromosome, ranges = IRanges(start = x$start, end = x$end)) 
  ol = as.matrix(findOverlaps(xRegions, binRegions, type = "any"))
  result = x[ol[, 1], ]
} 


hmfTheme = theme_bw() +
  theme(
    axis.title = element_text(size=7),
    axis.text = element_text(size=5),
    axis.ticks = element_blank(),
    legend.title = element_text(size=5), 
    legend.text = element_text(size=5),
    legend.key.size = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size=6)
  )
theme_set(hmfTheme)  


to_gc_bucket <- function(x) {
  cut(x, breaks = c(-0.1, seq(0.26, 0.64, 0.02), 1.1))
}

from_gc_bucket <- function(x) {
  (as.numeric(x) - 1) * 0.02 + 0.25
}

to_rep_bucket <- function(x) {
  cut(x, breaks = c(-0.1, seq(0.062, 0.78, 0.02), 1.1))
}

from_rep_bucket <- function(x) {
  (as.numeric(x) - 1) * 0.02 + 0.061
}


normalise_value <- function(df, averageReplication, to_bucket) {

  dfHist = df %>%
    filter(value > 0) %>%
    mutate(bucket = to_bucket(value)) %>%
    group_by(feature, bucket) %>%
    summarise(count = n()) %>%
    filter(!is.na(bucket)) %>% 
    group_by(feature) %>% 
    mutate(
      actualDensity = count / sum(count))
  
  refHist = averageReplication %>% 
    #mutate(replication = replication / 100) %>% 
    #filter(replication >= 0.06, replication <= 0.8) %>%
    mutate(bucket = to_bucket(value)) %>%
    group_by(bucket) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(expectedDensity = count / sum(count)) %>%
    select(-count)
  
  # Verify all 1s
  dfHist %>% group_by(feature) %>% summarise(totalDensity = sum(actualDensity))
  sum(refHist$expectedDensity)
  
  normalisedHist = dfHist %>% inner_join(refHist, by = "bucket") %>%
    group_by(feature) %>%
    mutate(
      ratio = actualDensity / expectedDensity,
      sumProduct = sum(ratio * actualDensity),
      factor = ratio / sumProduct)
  
  return (normalisedHist)
}

plot_replication <- function(df) {
  ggplot(df) +
    geom_col(aes(x = bucketCentre, y = factor), width = 0.02, fill = singleBlue, position = "identity") +
    facet_wrap( ~feature, nrow = 4) + 
    xlab("Replication") + ylab("Factor")
}

load(file = "/Users/jon/hmf/analysis/svPaper/featuredBreakends.RData")  

averageGC_1k = read.table("~/hmf/resources/GC_profile.1000bp.cnp", sep = '\t', header = F, stringsAsFactors = F) %>%
  mutate(chromosome = V1, start = V2 + 1, end = start + 1000, gc = V3) %>%
  select(chromosome, start, end, gc) %>%
  filter(gc > -1)
averageGC_1k = filter_unmappable(averageGC_1k, mappability_100k)

averageReplication_1k = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>% 
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)
averageReplication_1k = filter_unmappable(averageReplication_1k, mappability_100k)

replicationBreakends = featuredBreakends %>% filter(replication >= 0.06, replication <= 0.8) %>% mutate(value = replication)
averageReplication_1k = averageReplication_1k %>% mutate(value = replication / 100.0) %>% filter(value >= 0.06, value <= 0.8)
replicationBreakendNormalised = normalise_value(replicationBreakends, averageReplication_1k, to_rep_bucket)
replicationBreakendNormalised$bucketCentre = from_rep_bucket(replicationBreakendNormalised$bucket)
plot_replication(replicationBreakendNormalised)

replicationFactor = replicationBreakendNormalised %>% select(feature, replicationBucket = bucket, replicationFactor = factor)
save(replicationFactor, file = "~/hmf/analysis/svPaper/replicationFactor.RData")


########## GC

gcBreakends = featuredBreakends %>% mutate(value = gc)
averageGC_1k = averageGC_1k %>% mutate(value = gc)
gcBreakendNormalised = normalise_value(gcBreakends, averageGC_1k, to_gc_bucket)
gcBreakendNormalised$bucketCentre = from_gc_bucket(gcBreakendNormalised$bucket)
plot_replication(gcBreakendNormalised) + xlab("GC")

gcFactor = gcBreakendNormalised %>% select(feature, gcBucket = bucket, gcFactor = factor)
save(gcFactor, file = "~/hmf/analysis/svPaper/gcFactor.RData")


########## GC NORMALISED
gcNormalisedBreakends = featuredBreakends %>%
  mutate(gcBucket = to_gc_bucket(gc)) %>%
  left_join(gcFactor, by = c("feature","gcBucket"))


replicationBreakends = gcNormalisedBreakends %>% filter(replication >= 0.06, replication <= 0.8) %>% mutate(value = replication)
averageReplication_1k = averageReplication_1k %>% mutate(value = replication / 100.0) %>% filter(value >= 0.06, value <= 0.8)
replicationBreakendNormalised = normalise_value(replicationBreakends, averageReplication_1k, to_rep_bucket)
replicationBreakendNormalised$bucketCentre = from_rep_bucket(replicationBreakendNormalised$bucket)
plot_replication(replicationBreakendNormalised) + xlab("GC")




