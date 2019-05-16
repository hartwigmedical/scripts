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


#load(file = "/Users/jon/hmf/analysis/cohort/cohort.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/hpcSvs.RData")
averageReplication_1k = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>% 
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)

df = hpcSvs %>% filter(ResolvedType=='SimpleSV',Type %in% c('DEL', 'DUP')) %>%
  mutate(
    subtype = "Unknown",
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DUP' & PosEnd-PosStart<=500, "ShortDup", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DUP' & PosEnd-PosStart>500 & PosEnd-PosStart<=80000, "MediumDup", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DUP' & PosEnd-PosStart>80000 & PosEnd-PosStart<=1e6, "LongDup", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DUP' & PosEnd-PosStart>1e6 & PosEnd-PosStart<=5e6, "VeryLongDup", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DUP' & PosEnd-PosStart>5e6, "SuperLongDup", subtype),
    
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DEL' & PosEnd-PosStart<=500, "ShortDel", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DEL' & PosEnd-PosStart>500 & PosEnd-PosStart<=10000, "MediumDel", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DEL' & PosEnd-PosStart>10000 & PosEnd-PosStart<=1e6, "LongDel", subtype),
    subtype = ifelse(ResolvedType=='SimpleSV' & Type=='DEL' & PosEnd-PosStart>1e6, "VeryLongDel", subtype)
    ) %>% 
  select(subtype, RepOriginStart, RepOriginEnd) %>% gather(type, replication, c(2,3))

replication_bucket <- function(x) {
  
  
}

dfHist = df %>%
  mutate(
    bucketFactor = cut(replication, breaks = seq(0, 1, 0.02)),
    bucket = (as.numeric(bucketFactor) - 1) * 0.02 + 0.01
  ) %>%
  group_by(subtype, bucket) %>%
  summarise(count = n()) %>%
  filter(!is.na(bucket)) %>% 
  group_by(subtype) %>% 
  mutate(
    actualDensity = count / sum(count), 
    bucket = as.character(bucket))


refHist = averageReplication_1k %>%
  mutate(
    bucketFactor = cut(replication / 100, breaks = seq(0, 1, 0.02)),
    bucket = (as.numeric(bucketFactor) - 1) * 0.02 + 0.01) %>%
  group_by(bucket) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(!is.na(bucket), bucket >= 0.06, bucket <= 0.8) %>%
  mutate(expectedDensity = count / sum(count), bucket = as.character(bucket)) %>%
  select(-count)
  
# Verify all 1s
dfHist %>% group_by(subtype) %>% summarise(totalDensity = sum(actualDensity))
sum(refHist$expectedDensity)

normalisedHist = dfHist %>% inner_join(refHist, by = "bucket") %>%
  mutate(factor = actualDensity / expectedDensity, bucket = factor(bucket)) %>%
  group_by(subtype) %>%
  mutate(normalisedDensity = factor / sum(factor))

str(normalisedHist)

ggplot(normalisedHist) +
  geom_bar(aes(x = factor(bucket), y = factor), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = seq(0.07, 0.79, 0.09)) +
  facet_wrap( ~subtype, nrow = 3) + xlab("Replication") + ylab("Factor")# + xlim(0.3, 0.9)

ggplot(dfHist) +
  geom_bar(aes(x = bucket, y = density), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = seq(0, 1, 0.1), labels = as.character(seq(0, 1, 0.1))) +
  facet_wrap( ~subtype, nrow = 3) + xlab("Replication") + ylab("Density")# + xlim(0.3, 0.9)

ggplot(refHist) +
  geom_bar(aes(x = bucket, y = density), stat = "identity", width = 1, fill = singleBlue, position = "stack") +
  scale_x_discrete(breaks = seq(0, 1, 0.1), labels = as.character(seq(0, 1, 0.1))) +
  xlab("Replication") + ylab("Density")# + xlim(0.3, 0.9)




MediumDupHist = replication_histogram(df %>% filter(subtype == "MediumDup"))

ggplot(df) + 
  geom_histogram(aes(replication, y=0.01 * ..density..), binwidth = 0.01, position = "stack", color = singleBlue) + 
  facet_wrap( ~subtype, nrow = 3) + xlab("Replication") + ylab("Density") + xlim(0.3, 0.9)



load(file = "/Users/jon/hmf/analysis/svPaper/averageReplication.RData")
df = bind_rows(averageReplication_100k %>% mutate(bucket = "100k"), averageReplication_1M %>% mutate(bucket = "1M"))
df = bind_rows(df, averageReplication_10M %>% mutate(bucket = "10M"))


ggplot(df) + 
  geom_histogram(aes(averageReplication, y=0.01 * ..density..), binwidth = 0.01, position = "stack", color = singleBlue) + 
  facet_wrap( ~bucket, nrow = 3) + xlab("Replication") + ylab("Density")





ggplot(replication %>% mutate(replication = replication / 100) %>%  filter(replication > 0.01) )  + 
  geom_histogram(aes(replication, y=0.01 * ..density..), binwidth = 0.01, position = "stack", color = singleBlue) + 
  xlab("Replication") + ylab("Density")



domain = seq(0, 1, 0.01)


head(jon)


replication_histogram <- function(df, binwidth = 0.01) {
  result = 
}




domain = seq(0 + binwidth / 2, 10 + binwidth / 2, binwidth)

myInitialSomatics = somatics %>% filter(somaticPloidy > (binwidth / 2), somaticPloidy < (10 + binwidth / 2))
myEmptyHist = data.frame(bucket = domain, n = 0)

myInitialHist = myInitialSomatics %>%
  mutate(
    bucketFactor = cut(somaticPloidy, breaks = domain),
    bucket = (as.numeric(bucketFactor) - 1) * binwidth + binwidth / 2
  ) %>%
  group_by(bucket) %>%
  count() %>%
  bind_rows(myEmptyHist) %>%
  group_by(bucket) %>%
  summarise(n = sum(n)) %>%
  arrange(bucket)





load(file = "~/hmf/analysis/cohort/hpcCopyNumbers.RData")
head(hpcCopyNumbers)

jon = hpcCopyNumbers %>% arrange(-copyNumber)
head(jon)
