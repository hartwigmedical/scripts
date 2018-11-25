rm(list=setdiff(ls(), "allSomatics_p1"))

library(tidyr)
library(dplyr)
library(GenomicRanges)

sampleIdMap = read.csv(file = "~/hmf/secure/SampleIdMap.csv", stringsAsFactors = F)

bedFile = read.table(file = "~/hmf/analysis/polak/hg19.1Mb.sorted.bed", header = F) %>% 
  select(-V5) %>%
  mutate(V1 = substring(V1, 4))
colnames(bedFile) <- c("chromosome", "start", "end", "bin")

bedOverlap <- function(bedFile, somatics) {
  require(GenomicRanges)
  
  somatics = somatics %>% filter(type == 'SNP')
  somaticsRange = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position))
  
  bedFileRange = GRanges(bedFile$chromosome, IRanges(bedFile$start,bedFile$end))
  
  ol = data.frame(findOverlaps(bedFileRange, somaticsRange, type="any", select="all"))
  ol$bin = bedFile[ol$queryHits, "bin"]
  ol$sampleId = somatics[ol$subjectHits, "sampleId"]
  ol = ol %>% 
    select(sampleId, bin) %>%
    group_by(sampleId, bin) %>% 
    count() %>%
    ungroup()
  
  return (ol)
} 

#load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
bedOverlap_p1 = bedOverlap(bedFile, allSomatics_p1) %>%
  left_join(sampleIdMap, by = "sampleId") %>%
  select(hmfSampleId, bin, n)
save(bedOverlap_p1, file = "~/hmf/analysis/polak/bedOverlap_p1.RData")
rm(allSomatics_p1)


#load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
bedOverlap_p2 = bedOverlap(bedFile, allSomatics_p2) %>%
  left_join(sampleIdMap, by = "sampleId") %>%
  select(hmfSampleId, bin, n)
save(bedOverlap_p2, file = "~/hmf/analysis/polak/bedOverlap_p2.RData")
rm(allSomatics_p2)

load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary = highestPurityCohortSummary %>% left_join(sampleIdMap, by = "sampleId") %>%
  select(hmfSampleId, cancerType)

bedOverlapComplete = bind_rows(bedOverlap_p1, bedOverlap_p2) %>% inner_join(highestPurityCohortSummary, by = "hmfSampleId") 
write.csv(bedOverlapComplete, file = "~/hmf/analysis/polak/BedOverlap.csv", row.names = F) 


head(bedOverlapComplete)
