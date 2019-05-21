library(tidyr)
library(dplyr)
library(ggplot2)

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

averageReplication_1k = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>% 
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)

load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
allowedBins = mappability_100k %>% mutate(mappable = (100000 - N) / 100000) %>% filter(mappable >= 0.8)
binRegions = GRanges(allowedBins$chromosome, ranges = IRanges(start = allowedBins$start, end = allowedBins$end)) 
replicationRegions = GRanges(averageReplication_1k$chromosome, ranges = IRanges(start = averageReplication_1k$start, end = averageReplication_1k$end)) 
ol = as.matrix(findOverlaps(replicationRegions, binRegions, type = "any"))
averageReplication_1k = averageReplication_1k[ol[, 1], ]
rm(list=setdiff(ls(), c("averageReplication_1k")))

load(file = "/Users/jon/hmf/analysis/svPaper/hpcDelsDups.RData")

df = hpcDelsDups %>% 
  filter(ResolvedType=='SimpleSV', Type %in% c('DEL', 'DUP')) %>%
  select(subtype, RepOriginStart, RepOriginEnd) %>% gather(type, replication, c(2,3)) %>%
  filter(replication >= 0.06, replication <= 0.8)


replication_bucket <- function(replication) {
  cut(replication, breaks = c(0, seq(0.062, 0.78, 0.02), 1))
}

dfHist = df %>%
  filter(replication > 0) %>%
  mutate(bucket = replication_bucket(replication)) %>%
  group_by(subtype, bucket) %>%
  summarise(count = n()) %>%
  filter(!is.na(bucket)) %>% 
  group_by(subtype) %>% 
  mutate(
    actualDensity = count / sum(count))

refHist = averageReplication_1k %>% 
  mutate(replication = replication / 100) %>% 
  filter(replication >= 0.06, replication <= 0.8) %>%
  mutate(bucket = replication_bucket(replication)) %>%
  group_by(bucket) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(expectedDensity = count / sum(count)) %>%
  select(-count)
  
# Verify all 1s
dfHist %>% group_by(subtype) %>% summarise(totalDensity = sum(actualDensity))
sum(refHist$expectedDensity)

normalisedHist = dfHist %>% inner_join(refHist, by = "bucket") %>%
  group_by(subtype) %>%
  mutate(
    ratio = actualDensity / expectedDensity,
    sumProduct = sum(ratio * actualDensity),
    factor = ratio / sumProduct) 
    
verification = normalisedHist %>% group_by(subtype) %>% summarise(jon = sum(factor * count / sum(count)))
replicationFactor = normalisedHist %>% select(subtype, replicationBucket = bucket, replicationFactor = factor)
save(replicationFactor, file = "~/hmf/analysis/svPaper/replicationFactor.RData")


ggplot(normalisedHist) +
  geom_bar(aes(x = bucket, y = factor), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = unique(normalisedHist$bucket)[seq(1,37, 9)], labels = (seq(1,37, 9) * 0.02 + 0.05)) +
  facet_wrap( ~subtype, nrow = 3) + xlab("Replication") + ylab("Factor")

ggplot(dfHist) +
  geom_bar(aes(x = bucket, y = density), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = seq(0, 1, 0.1), labels = as.character(seq(0, 1, 0.1))) +
  facet_wrap( ~subtype, nrow = 3) + xlab("Replication") + ylab("Density")# + xlim(0.3, 0.9)

ggplot(refHist) +
  geom_bar(aes(x = bucket, y = density), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = seq(0, 1, 0.1), labels = as.character(seq(0, 1, 0.1))) +
  xlab("Replication") + ylab("Density")# + xlim(0.3, 0.9)


