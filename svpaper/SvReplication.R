library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

singleBlue = "#6baed6"

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


normalise_replication <- function(df, averageReplication) {
  replication_bucket <- function(replication) {
    cut(replication, breaks = c(0, seq(0.062, 0.78, 0.02), 1))
  }
  
  dfHist = df %>%
    filter(replication > 0) %>%
    mutate(bucket = replication_bucket(replication)) %>%
    group_by(feature, bucket) %>%
    summarise(count = n()) %>%
    filter(!is.na(bucket)) %>% 
    group_by(feature) %>% 
    mutate(
      actualDensity = count / sum(count))
  
  refHist = averageReplication %>% 
    mutate(replication = replication / 100) %>% 
    filter(replication >= 0.06, replication <= 0.8) %>%
    mutate(bucket = replication_bucket(replication)) %>%
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
    geom_bar(aes(x = bucket, y = factor), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
    scale_x_discrete(breaks = unique(df$bucket)[seq(1,37, 9)], labels = (seq(1,37, 9) * 0.02 + 0.05)) +
    facet_wrap( ~feature, nrow = 4) + 
    xlab("Replication") + ylab("Factor")
}

averageReplication_1k = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>% 
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)

load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
allowedBins = mappability_100k %>% mutate(mappable = (100000 - N) / 100000) %>% filter(mappable >= 0.8)
binRegions = GRanges(allowedBins$chromosome, ranges = IRanges(start = allowedBins$start, end = allowedBins$end)) 
replicationRegions = GRanges(averageReplication_1k$chromosome, ranges = IRanges(start = averageReplication_1k$start, end = averageReplication_1k$end)) 
ol = as.matrix(findOverlaps(replicationRegions, binRegions, type = "any"))
averageReplication_1k = averageReplication_1k[ol[, 1], ]


load(file = "/Users/jon/hmf/analysis/svPaper/featuredBreakends.RData")  
replicationBreakends = featuredBreakends %>% filter(replication >= 0.06, replication <= 0.8)
replicationBreakendNormalised = normalise_replication(replicationBreakends, averageReplication_1k)
plot_replication(replicationBreakendNormalised)

replicationFactor = replicationBreakendNormalised %>% select(feature, replicationBucket = bucket, replicationFactor = factor)
save(replicationFactor, file = "~/hmf/analysis/svPaper/replicationFactor.RData")
