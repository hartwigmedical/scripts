library(tidyr)
library(dplyr)
library(RMySQL)
library(purple)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


####### PREPARE DATA
load(file = "~/data/sv/cohort.RData")
View(multipleBiopsyCohort)

allSvs = read.csv("~/data/sv/SVA_SVS.csv")
mbcSvs = allSvs %>% filter(SampleId %in% multipleBiopsyCohort$sampleId)
nrow(mbcSvs)
save(mbcSvs, file = "~/data/sv/mbcSvs.RData")

doubleBiopsy = multipleBiopsyCohort %>% group_by(patientId) %>% filter(n() == 2) %>% mutate(n = row_number())
dbSvs = mbcSvs %>% filter(SampleId %in% doubleBiopsy$sampleId)
nrow(dbSvs)

reformattedSVS = dbSvs %>% select(sampleId = SampleId, everything()) %>% left_join(doubleBiopsy %>% select(sampleId, n), by = "sampleId")

query = reformattedSVS %>% filter(n == 1)
subject = reformattedSVS %>% filter(n == 2)
nrow(query)
nrow(subject)
View(query)

# tidy up preparation objects
rm(reformattedSVS)
rm(dbSvs)
rm(allSvs)
rm(mbcSvs)


save(query, subject, file = "~/data/sv/multipleBiopsy.RData")

####### ANALYSIS
sv_overlaps <- function(query, subject, maxgap = -1) {
  require(tidyr)
  require(dplyr)
  require(GenomicRanges)
  
  queryStartRange <- GRanges(paste0(query$patientId, query$ChrStart, query$OrientStart), IRanges(query$PosStart, query$PosStart))
  subjectStartRange <- GRanges(paste0(subject$patientId, subject$ChrStart, subject$OrientStart), IRanges(subject$PosStart, subject$PosStart))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = maxgap))
  
  queryEndRange <- GRanges(paste0(query$patientId, query$ChrEnd, query$OrientEnd), IRanges(query$PosEnd, query$PosEnd))
  subjectEndRange <- GRanges(paste0(subject$patientId, subject$ChrEnd, subject$OrientEnd), IRanges(subject$PosEnd, subject$PosEnd))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = maxgap))
  
  overlaps = inner_join(startOverlaps, endOverlaps, by = c("queryHits", "subjectHits"))
  
  overlapQueryData = query[overlaps$queryHits, ] %>%
    mutate(queryHits = overlaps$queryHits) %>%
    select(queryHits, sampleId, ChrStart, ChrEnd, PosStart, PosEnd, OrientStart, OrientEnd)
  
  overlapSubjectData = subject[overlaps$subjectHits, ] %>%
    mutate(subjectHits = overlaps$subjectHits) %>%
    select(subjectHits, subjectPosStart = PosStart, subjectPosEnd = PosEnd, subjectOrientStart = OrientStart, subjectOrientEnd = OrientEnd)
  
  overlapsData = bind_cols(overlapQueryData, overlapSubjectData) %>%
    filter(OrientStart == subjectOrientStart, OrientEnd == subjectOrientEnd) %>%
    select(-subjectOrientStart, -subjectOrientEnd) %>%
    mutate(PosStartDiff = abs(PosStart - subjectPosStart), PosEndDiff = abs(PosEnd - subjectPosEnd), positionDiff = PosStartDiff + PosEndDiff) %>%
    group_by(ChrStart, ChrEnd, PosStart, PosEnd,OrientStart,OrientEnd) %>%
    top_n(1, -positionDiff) %>%
    group_by(queryHits) %>%
    top_n(1, -subjectHits)
  
  return (overlapsData %>% select(queryHits, subjectHits))
}

load(file = "~/data/sv/multipleBiopsy.RData")
overlaps = sv_overlaps(query, subject, 10)
nrow(overlaps)
View(overlaps)

subject$scope <- "Private"
query$scope <- "Private"
query$ClusterCountMatch <- NA
query$ClusterCountDesc <- NA

subject[overlaps$subjectHits, "scope"] <- "Shared"
query[overlaps$queryHits, "scope"] <- "Shared"

query[overlaps$queryHits, "ClusterCountMatch"] <- query[overlaps$queryHits, "ClusterCount"] == subject[overlaps$subjectHits, "ClusterCount"]
query[overlaps$queryHits, "ClusterDescMatch"] <- query[overlaps$queryHits, "ClusterDesc"] == subject[overlaps$subjectHits, "ClusterDesc"]
query[overlaps$queryHits, "ResolvedTypeMatch"] <- query[overlaps$queryHits, "ResolvedType"] == subject[overlaps$subjectHits, "ResolvedType"]
query[overlaps$queryHits, "ResolvedTypeOther"] <- subject[overlaps$subjectHits, "ResolvedType"]

query[overlaps$queryHits, "SubjectClusterCount"] <- subject[overlaps$subjectHits, "ClusterCount"]
query[overlaps$queryHits, "SubjectClusterId"] <- subject[overlaps$subjectHits, "ClusterId"]
subject[overlaps$subjectHits, "QueryClusterCount"] <- query[overlaps$queryHits, "ClusterCount"]
subject[overlaps$subjectHits, "QueryClusterId"] <- query[overlaps$queryHits, "ClusterId"]

subjectClusterCount = subject %>% select(patientId, SubjectClusterId = ClusterId, SubjectClusterCount = ClusterCount) %>% distinct()
subjectClusterCount %>% group_by(patientId, SubjectClusterId) %>% count() %>% filter(n > 1)

View(subject)

queryClusterMap = query %>% 
  filter(scope == 'Shared') %>% 
  select(sampleId, patientId, ClusterId, ClusterCount, SubjectClusterId) %>%
  distinct() %>%
  left_join(subjectClusterCount, by = c("patientId","SubjectClusterId")) %>%
  group_by(sampleId, ClusterId, ClusterCount) %>% 
  summarise(SubjectClusterIdCount = n(), SubjectClusterTotalCount = sum(SubjectClusterCount))

queryClusterStatus = query %>% 
  group_by(sampleId, ClusterId, ClusterCount, ResolvedType, scope) %>% 
  summarise(n = n()) %>% ungroup() %>%
  spread(scope, n, fill = 0) %>%
  left_join(queryClusterMap, by = c("sampleId", "ClusterId", "ClusterCount")) %>%
  mutate(
    SubjectClusterIdCount = ifelse(is.na(SubjectClusterIdCount), 0, SubjectClusterIdCount),
    SubjectClusterTotalCount = ifelse(is.na(SubjectClusterTotalCount), 0, SubjectClusterTotalCount)
  ) %>%
  mutate(
    status = "Mixed",
    status = ifelse(Shared == 0, "Private", status),
    status = ifelse(Private == 0 & SubjectClusterIdCount == 1 & ClusterCount == SubjectClusterTotalCount, "Exact", status),
    status = ifelse(Private > 0 & SubjectClusterIdCount == 1 & Shared == SubjectClusterTotalCount, "SimpleSuperset", status),
    status = ifelse(Private >= 0 & SubjectClusterIdCount > 1 & Shared == SubjectClusterTotalCount, "ComplexSuperset", status),
    status = ifelse(Private == 0 & SubjectClusterIdCount == 1 & ClusterCount < SubjectClusterTotalCount, "Subset", status)
  ) 

queryClusterStatusSummary = queryClusterStatus %>% group_by(ResolvedType,status) %>% count() %>% spread(status,n)
View(queryClusterStatusSummary)

subjectClusterMap = subject %>% 
  filter(scope == 'Shared') %>% 
  select(sampleId, ClusterId, ClusterCount, QueryClusterId) %>% 
  distinct() %>% 
  arrange(sampleId, ClusterId, ClusterCount, QueryClusterId) %>%
  group_by(sampleId, ClusterId, ClusterCount) %>% summarise(OtherClusterIdCount = n(), QueryClusterId = paste0(QueryClusterId, collapse = ",")) %>% 
  select(sampleId,  ClusterId,  OtherClusterIdCount)

queryClusterStatus = query %>% 
  group_by(sampleId,ClusterId,ClusterCount,SubjectClusterCount, ResolvedType,scope) %>% 
  count() %>% ungroup() %>%
  spread(scope, n, fill = 0) %>% mutate(SubjectClusterCount = ifelse(is.na(SubjectClusterCount), 0, SubjectClusterCount)) %>%
  left_join(queryClusterMap, by = c("sampleId", "ClusterId")) %>%
  mutate(OtherClusterIdCount = ifelse(is.na(OtherClusterIdCount), 0, OtherClusterIdCount)) %>%
  mutate(
    status = "Mixed",
    status = ifelse(Shared == 0, "Private", status),
    status = ifelse(Private == 0 & OtherClusterIdCount == 1 & ClusterCount == SubjectClusterCount, "Exact", status),
    status = ifelse(Private > 0 & OtherClusterIdCount == 1 & Shared == SubjectClusterCount, "Superset", status),
    status = ifelse(Private == 0 & OtherClusterIdCount == 1 & ClusterCount < SubjectClusterCount, "Subset", status)
  )




View(query %>% filter(scope=='Shared',ResolvedTypeMatch==F) %>% group_by(ResolvedType,scope, ResolvedTypeOther) %>% count() %>% spread(ResolvedTypeOther,n))
View(query %>% filter(scope=='Shared') %>% group_by(ResolvedType,scope, ResolvedTypeOther) %>% count() %>% spread(ResolvedTypeOther,n))
View(query %>% filter(scope=='Shared') %>% group_by(ClusterCount,scope, ClusterCountOther) %>% count() %>% spread(ClusterCountOther,n))
View(query %>% group_by(sampleId,ClusterId,ClusterCount,scope) %>% tally %>% spread(scope,n) %>% filter(!is.na(Private),!is.na(Shared)))

#private, Exact,Subset,Superset





####### TEMP

query_patient_id_lookup<-function(dbConnect) {
  query = paste(
    "SELECT CPCTCNT.patientId AS sampleId, concat('CPCT02', CPCTCNT.itemValue, LPAD(RIGHT(CPCTPN.itemValue,4), 4, '0')) as patientId",
    "  FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN ",
    " WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId ",
    "   AND CPCTCNT.item = 'FLD.CPCTCNT' AND CPCTCNT.itemValue != ''",
    "   AND CPCTPN.item = 'FLD.CPCTPN' AND CPCTPN.itemValue != ''",
    "   AND CTCT2YN.item = 'FLD.CTCT2YN' AND CTCT2YN.itemValue = 'Yes'",
    sep = "")
  
  result = dbGetQuery(dbConnect, query)
  return (result)
}


sample_to_patient_id<-function(sampleId, lookup) {
  colnames(lookup) <- c("truncatedSampleIds", "patientIds")
  
  lookup = rbind(manual_patient_id(), lookup)
  
  index = match(substr(sampleId, 1, 12) , substr(lookup[[1]], 1, 12))
  if (is.na(index)) {
    substr(sampleId, 1, 12)
  } else {
    lookup[index, c(2)]
  }
}

manual_patient_id<-function() {
  truncatedSampleIds  = c("CPCT02020192", "CPCT02030224", "DRUP01010007", "DRUP01070024", "DRUP01050008",
                          "DRUP01010065", "DRUP01330002", "DRUP01340004", "DRUP01340003", "DRUP01340002", "DRUP01070008")
  patientIds = c("CPCT02020438", "CPCT02030292", "DRUP01010044", "CPCT02070110", "CPCT02050116",
                 "CPCT02010639", "CPCT02330049", "CPCT02340029", "CPCT02340014", "CPCT02340026", "CPCT02070023")
  return (data.frame(truncatedSampleIds, patientIds, stringsAsFactors = FALSE))
}

