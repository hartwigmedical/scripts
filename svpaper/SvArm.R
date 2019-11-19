library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

localPath = '~/hmf/analyses/SVAnalysis/'
sharedPath = '~/Dropbox/HMF Australia team folder/RData/'
hmftoolsPath = '~/hmf/repos/hmftools/'

##### ARM locations
centromeres = read.csv(paste0(hmftoolsPath,"hmf-common/src/main/resources/refgenome/hg19_centromere.tsv"), sep = "\t", header= F)
names(centromeres) <- c("chromosome", "centromere")

lengths = read.csv(paste0(hmftoolsPath,"hmf-common/src/main/resources/refgenome/hg19_len.tsv"), sep = "\t", header= F)
names(lengths) <- c("chromosome", "length")

chromosomes = left_join(centromeres, lengths, by = "chromosome") %>%
  mutate(pStart = 1, pEnd = centromere, qStart =centromere +1, qEnd = length, pLength = pEnd - pStart + 1, qLength = qEnd - qStart + 1) %>%
  select(chromosome, pStart, pEnd, qStart, qEnd, pLength, qLength)

arms = bind_rows(
  chromosomes %>% select(chromosome, start = pStart, end = pEnd, armLength = pLength) %>% mutate(arm = paste0("p", chromosome)),
  chromosomes %>% select(chromosome, start = qStart, end = qEnd, armLength = qLength) %>% mutate(arm = paste0("q", chromosome)))

rm(centromeres, lengths, chromosomes)

##### AVERAGE COPY NUMBER
load(file = paste0(sharedPath,"hpcCopyNumbers.RData"))
armRegions = GRanges(arms$chromosome, ranges = IRanges(start = arms$start, end = arms$end)) 
copyNumberRegions = GRanges(hpcCopyNumbers$chromosome, ranges = IRanges(start = hpcCopyNumbers$start, end = hpcCopyNumbers$end)) 
ol = as.matrix(findOverlaps(armRegions, copyNumberRegions, type = "any"))
hpcCopyNumbers[ol[, 2], "arm"] = arms[ol[, 1], "arm"]

armAveragedCopyNumber = hpcCopyNumbers %>% 
  group_by(sampleId, arm) %>%
  mutate(length = end - start + 1) %>%
  summarise(averageCopyNumber = sum(length * copyNumber) / sum(length)) %>%
  group_by(arm) %>%
  summarise(averageCopyNumber = sum(averageCopyNumber) / n())

##### FUNCTIONS
count_summmary <- function(x,byCluster=T) {
  if(byCluster) {
      x=x %>% group_by(SampleId, feature, ClusterId, Chr, arm, armLength) %>% count() %>% ungroup() 
  }
  x %>%
    filter(!Chr %in% c('X','Y'), !arm %in% c("p13","p14","p15","p20","p21")) %>%
    left_join(armAveragedCopyNumber, by = "arm") %>%
    group_by(feature, arm, armLength) %>%
    summarise(n = n(), nAdj = sum(1 / averageCopyNumber))
}
plot_arm_by_feature <- function(x) {
  ggplot(x, aes(x = armLength, y = nAdj)) + 
    geom_point() + 
    geom_smooth(method=lm,se=TRUE, fullrange=F) +
    geom_text(aes(label = arm), size = 2, nudge_x = 7000000, color = "red") +
    theme(panel.grid.major.y=element_line(linetype = 8,size=0.1),panel.grid.major.x=element_line(linetype = 8,size=0.1)) +
    facet_wrap(~feature, scales = "free_y") + 
    expand_limits(x = 0, y = 0)
}


##### Annotate FEATURES
load(file = paste0(localPath,"featuredBreakends.RData"))  

armData = featuredBreakends %>% 
  mutate(
    feature = ifelse(grepl("Foldback", feature), "Foldback", as.character(feature)),
    feature = ifelse(grepl("Dup", feature), "Dup", as.character(feature)),
    feature = ifelse(grepl("Del", feature), "Del", as.character(feature))) %>%
  filter(IsStart | !feature %in% c("Dup","Del")) 
armRegions = GRanges(arms$chromosome, ranges = IRanges(start = arms$start, end = arms$end)) 
armDataRegions = GRanges(armData$Chr, ranges = IRanges(start = armData$Pos, end = armData$Pos)) 
ol = as.matrix(findOverlaps(armRegions, armDataRegions, type = "any"))
armData[ol[, 2], "armLength"] = arms[ol[, 1], "armLength"]
armData[ol[, 2], "arm"] = arms[ol[, 1], "arm"]

armDataSummary = count_summmary(armData,byCluster=T)
plot_arm_by_feature(armDataSummary)

# FOLDBACKS and centromeric SGLs tend to occur on the same arms...
ggplot(armDataSummary %>% select(feature,n,arm) %>% filter(feature=='CentromericSGL'|feature=='Foldback') %>% spread(feature,n), aes(x = Foldback, y = CentromericSGL)) + 
  geom_point() + 
  geom_smooth(method=lm,se=TRUE, fullrange=F) +
  geom_text(aes(label = arm), size = 2, nudge_x = 5, color = "red") +
  theme(panel.grid.major.y=element_line(linetype = 8,size=0.1),panel.grid.major.x=element_line(linetype = 8,size=0.1)) +
  expand_limits(x = 0, y = 0)

#### FEATURES PER CLUSTERID / ARM
load(paste0(localPath,"resolveTypeBreakends.RData"))

armData = resolveTypeBreakends %>% mutate(feature = ResolvedType) %>% filter(IsStart | !feature %in% c("DUP","DEL","INV"))
armRegions = GRanges(arms$chromosome, ranges = IRanges(start = arms$start, end = arms$end)) 
armDataRegions = GRanges(armData$Chr, ranges = IRanges(start = armData$Pos, end = armData$Pos)) 
ol = as.matrix(findOverlaps(armRegions, armDataRegions, type = "any"))
armData[ol[, 2], "armLength"] = arms[ol[, 1], "armLength"]
armData[ol[, 2], "arm"] = arms[ol[, 1], "arm"]

armDataSummary = count_summmary(armData,byCluster = T)
plot_arm_by_feature(armDataSummary)

#### Complex
load(file = paste0(sharedPath,"highestPurityCohort.RData"))

armData = resolveTypeBreakends  %>% filter(ResolvedType == 'COMPLEX', !ClusterContainsDriver) %>% 
  left_join(highestPurityCohort %>% select(sampleId, cancerType), by = c("SampleId" = "sampleId")) %>% mutate(feature = cancerType) 
  #mutate(ClusterVariantCount = cut(ClusterVariantCount, c(0, 2, 4, 8, 16, 32, 64, 128, 8192), include.lowest = F, ordered = T)) %>%
  #mutate(feature = ifelse(ClusterContainsBND > 0, "NonDriverWithBND", "NonDriverNoBND"))
  #mutate(feature = paste0("ClusterVariantCount", as.character(ClusterVariantCount)))
armRegions = GRanges(arms$chromosome, ranges = IRanges(start = arms$start, end = arms$end)) 
armDataRegions = GRanges(armData$Chr, ranges = IRanges(start = armData$Pos, end = armData$Pos)) 
ol = as.matrix(findOverlaps(armRegions, armDataRegions, type = "any"))
armData[ol[, 2], "armLength"] = arms[ol[, 1], "armLength"]
armData[ol[, 2], "arm"] = arms[ol[, 1], "arm"]

armDataSummary = count_summmary(armData,byCluster = T)
plot_arm_by_feature(armDataSummary)

View(merge(featuredBreakends %>% filter(feature=='TelomericSGL') %>% group_by(SampleId) %>% count() %>% arrange(-n),highestPurityCohort %>% select(SampleId=sampleId,cancerType),by='SampleId' ))
View(svDriver)

View(merge(featuredBreakends %>% filter(feature=='TelomericSGL') %>% group_by(SampleId) %>% count() %>% filter(n>20) %>%
        arrange(-n),svDrivers %>% filter(DriverLikelihood>0.5) %>% select(SampleId,Gene),by='SampleId') %>% group_by(Gene) %>% summarise(sum(n),n()))


View(merge(svDrivers,highestPurityCohort %>% select(SampleId=sampleId,cancerType),by='SampleId') %>% filter(cancerType=='Bone/Soft tissue'))# %>% group_by(Gene) %>% count)
View(svDrivers)
