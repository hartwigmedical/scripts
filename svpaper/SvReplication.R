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


normalise_replication <- function(df, averageReplication) {
  replication_bucket <- function(replication) {
    cut(replication, breaks = c(0, seq(0.062, 0.78, 0.02), 1))
  }
  
  dfHist = df %>%
    filter(replication > 0) %>%
    mutate(bucket = replication_bucket(replication)) %>%
    group_by(ResolvedType, bucket) %>%
    summarise(count = n()) %>%
    filter(!is.na(bucket)) %>% 
    group_by(ResolvedType) %>% 
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
  dfHist %>% group_by(ResolvedType) %>% summarise(totalDensity = sum(actualDensity))
  sum(refHist$expectedDensity)
  
  normalisedHist = dfHist %>% inner_join(refHist, by = "bucket") %>%
    group_by(ResolvedType) %>%
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
    facet_wrap( ~ResolvedType, nrow = 3) + 
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


#allData = svData %>% 
#  select(ResolvedType, RepOriginStart, RepOriginEnd, PosStart, PosEnd, ClusterCount, FoldbackLnkStart, FoldbackLnkEnd, LocTopTypeStart, LocTopTypeEnd, LEStart, LEEnd) %>% 
#  gather(type, replication, c(2,3)) %>%
#  filter(replication >= 0.06, replication <= 0.8)

load(file = "/Users/jon/hmf/analysis/svPaper/svData.RData")
beData = rbind(
  svData %>% mutate(IsFoldback = FoldbackLnkStart == FoldbackLnkEnd, FoldbackLnk=FoldbackLnkStart,LocTopType=LocTopTypeStart,Chr=ChrStart,Pos=PosStart,Orient=OrientStart,IsStart=T,LE=LEStart, replication = RepOriginStart, Len=PosEnd - PosStart + 1),
  svData %>% mutate(IsFoldback = FoldbackLnkStart == FoldbackLnkEnd, FoldbackLnk=FoldbackLnkEnd,LocTopType=LocTopTypeEnd,Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,IsStart=F, LE=LEEnd, replication = RepOriginEnd, Len=PosEnd - PosStart + 1)) %>%
  select(SampleId,Id,IsStart,ClusterId,Type,ResolvedType,FoldbackLnk,LocTopType,Chr,Pos,Orient,LE,ClusterCount,Len,replication, IsFoldback) %>%
  filter(replication >= 0.06, replication <= 0.8)

unique(allData$ResolvedType)

resolveTypeDF = beData %>% mutate(ResolvedType = ifelse(ResolvedType %in% c('DEL','DUP','LINE','COMPLEX'),as.character(ResolvedType), "Other"))
resolveTypeNormalised = normalise_replication(resolveTypeDF, averageReplication_1k)
plot_replication(resolveTypeNormalised)

 

delDF = beData %>% filter(ResolvedType == 'DEL', ClusterCount == 1) %>%
  mutate(
    length = Len,
    label = cut(length, c(0, 500, 10000, 500000, 1e100), labels = c("Short", "Medium", "Long", "VeryLong")),
    ResolvedType = paste0(label, "Del"),
    ResolvedType = factor(ResolvedType, levels = c('ShortDel', "MediumDel", "LongDel", "VeryLongDel"), ordered = T) )
delNormalised = normalise_replication(delDF, averageReplication_1k)
plot_replication(delNormalised)

dupDF = beData %>% filter(ResolvedType == 'DUP', ClusterCount == 1) %>%
  mutate(
    length = Len,
    label = cut(length, c(0, 500, 10000, 500000, 1e100), labels = c("Short", "Medium", "Long", "VeryLong")),
    ResolvedType = paste0(label, "Dup"),
    ResolvedType = factor(ResolvedType, levels = c('ShortDup', "MediumDup", "LongDup", "VeryLongDup"), ordered = T) )
dupNormalised = normalise_replication(dupDF, averageReplication_1k)
plot_replication(dupNormalised)


featuresDF = beData %>% mutate(
  ResolvedType = as.character(ResolvedType),
  ResolvedType = ifelse(ResolvedType == 'INV' & IsFoldback, "Foldback", ResolvedType),
  ResolvedType = ifelse(LocTopType=='TI_ONLY', "TISource", ResolvedType),
  ResolvedType = ifelse(ResolvedType == 'LINE' & LE == "None", "LineInsertion", ResolvedType)) %>%
  filter(ResolvedType %in% c("Foldback", "TISource","LineInsertion"))
featuresNormalised = normalise_replication(featuresDF, averageReplication_1k)
plot_replication(featuresNormalised)

pdf("~/hmf/analysis/svPaper/plot/replication.pdf")
plot_replication(resolveTypeNormalised)
plot_replication(delNormalised)
plot_replication(dupNormalised)
plot_replication(featuresNormalised)
dev.off()

########## JUNK?

verification = normalisedHist %>% group_by(ResolvedType) %>% summarise(total = sum(factor * count / sum(count)))
replicationFactor = normalisedHist %>% select(ResolvedType, replicationBucket = bucket, replicationFactor = factor)
save(replicationFactor, file = "~/hmf/analysis/svPaper/replicationFactor.RData")


types = unique(svData$ResolvedType)
for (type in types) {
  cat("Processing", type, "\n")
  replicationPlot = plot_replication(normalisedHist %>% filter(grepl(type, ResolvedType)))
  filename = paste0("~/hmf/analysis/svPaper/plot/replication", type, ".png")
  ggsave(filename, replicationPlot, width = 183, height = 161, units = "mm", dpi = 300)
}

all = plot_replication(normalisedHist)
all

delReplicationPlot = plot_replication(normalisedHist %>% filter(grepl("Del", subtype)))
delReplicationPlot
ggsave("~/hmf/analysis/svPaper/plot/delReplication.png", delReplicationPlot, width = 183, height = 120, units = "mm", dpi = 300)

dupReplicationPlot = plot_replication(normalisedHist %>% filter(grepl("Dup", subtype)))
ggsave("~/hmf/analysis/svPaper/plot/dupReplication.png", dupReplicationPlot, width = 183, height = 120, units = "mm", dpi = 300)


simplePlot = plot_replication(normalisedHist %>% filter(grepl("Dup", subtype) | grepl("Del", subtype)))
ggsave("~/hmf/analysis/svPaper/plot/simpleReplication.png", simplePlot, width = 183, height = 161, units = "mm", dpi = 300)
ggsave("~/hmf/analysis/svPaper/plot/simpleReplication.pdf", simplePlot, width = 183, height = 161, units = "mm", dpi = 300)

featurePlot = plot_replication(normalisedHist %>% filter(grepl("Fold", subtype) | grepl("ShortTI", subtype)))
ggsave("~/hmf/analysis/svPaper/plot/featureReplication.png", featurePlot, width = 183, height = 161, units = "mm", dpi = 300)
ggsave("~/hmf/analysis/svPaper/plot/featureReplication.pdf", featurePlot, width = 183, height = 161, units = "mm", dpi = 300)


load(file = "/Users/jon/hmf/analysis/svPaper/svData.RData")
df = svData %>% filter(isShortTI, !ResolvedType %in% c("COMPLEX", "LINE")) %>% filter(xor(LocTopTypeStart=='TI_ONLY', LocTopTypeEnd=='TI_ONLY')) %>%
  mutate(
    TIReplication = ifelse(LocTopTypeStart=='TI_ONLY', RepOriginStart, RepOriginEnd), 
    NonTIReplication = ifelse(LocTopTypeStart=='TI_ONLY', RepOriginEnd, RepOriginStart))

ggplot(df) +
    geom_point(aes(x = TIReplication, y = NonTIReplication, color = ResolvedType), size = 0.3) + ggtitle("ShortTI TI Vs Non TI Replication" ) 

