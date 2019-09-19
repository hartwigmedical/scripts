library(tidyr)
library(dplyr)

svDrivers = read.csv(file = "/Users/jon/hmf/analysis/fusions/SVA_DRIVERS.csv")
svDriverSummary = svDrivers %>% group_by(SampleId, ClusterId, ResolvedType) %>% count() %>% filter(ResolvedType != "") %>% select(-n) %>% mutate(IsDriver = T)

svData = read.csv(file = "/Users/jon/hmf/analysis/fusions/SVA_SVS.csv")
svData = left_join(svData, svDriverSummary, by = c("SampleId", "ClusterId", "ResolvedType")) 
svData = svData %>% 
  mutate(
    ClusterContainsDriver = !is.na(IsDriver), 
    ClusterContainsBND = grepl("BND", ClusterDesc)) %>%
  group_by(SampleId, ClusterId, ResolvedType) %>%
  mutate(ClusterVariantCount = n()) %>%
  ungroup() %>%
  select(-IsDriver)

beData = rbind(
  svData %>% mutate(LnkLen = LnkLenStart, IsFoldback = Type == 'INV' & FoldbackLnkStart == FoldbackLnkEnd & FoldbackLnkStart > 0, RepeatClass, RepeatType, FoldbackLnk=FoldbackLnkStart,LocTopType=LocTopTypeStart,Chr=ChrStart,Pos=PosStart,Orient=OrientStart,IsStart=T,LE=LEStart, replication = RepOriginStart, Len=PosEnd - PosStart + 1, PosMid = round(PosStart  + (PosEnd - PosStart)/2), CN = CNStart),
  svData %>% filter(ChrEnd != 0) %>% mutate(LnkLen = LnkLenEnd, IsFoldback = Type == 'INV' & FoldbackLnkStart == FoldbackLnkEnd & FoldbackLnkStart > 0, RepeatClass, RepeatType, FoldbackLnk=FoldbackLnkEnd,LocTopType=LocTopTypeEnd,Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,IsStart=F, LE=LEEnd, replication = RepOriginEnd, Len=PosEnd - PosStart + 1, PosMid = round(PosStart  + (PosEnd - PosStart)/2), CN = CNEnd)) %>%
  select(SampleId,Id,IsStart,ClusterId,Type,ResolvedType,RepeatType,RepeatClass,FoldbackLnk,LocTopType,Chr,Pos,Orient,LE,ClusterCount,Len,replication, IsFoldback, LnkLen, CN, ClusterContainsDriver, ClusterContainsBND, ClusterVariantCount) 
  
lengthLabels = c("Short", "Medium", "Long", "VeryLong") 
lengthBreaks = c(0, 500, 10000, 500000, 1e100)
dupLabels = paste0(lengthLabels, "Dup")
delLabels = paste0(lengthLabels, "Del")
featureLevels = c(dupLabels, delLabels, "TiSource", "CentromericSGL", "PeriCentromericSGL", "TelomericSGL", "LeftFoldback", "RightFoldback", "LineInsertion", "LineSource")

dupBreakends = beData %>% filter(ResolvedType == 'DUP', ClusterCount == 1) %>% 
  mutate(feature = cut(Len, lengthBreaks, labels = paste0(lengthLabels, "Dup"), ordered_result = F))
delBreakends = beData %>% filter(ResolvedType == 'DEL', ClusterCount == 1) %>% 
  mutate(feature = cut(Len, lengthBreaks, labels = paste0(lengthLabels, "Del"), ordered_result = F))
sglBreakends = beData %>% 
  mutate(
    feature = ifelse(RepeatClass == 'Satellite/centr', "CentromericSGL", NA), 
    feature = ifelse(RepeatType == 'HSATII', "PeriCentromericSGL", feature),
    feature = ifelse(RepeatType=='(CCCTAA)n'|RepeatType=='(TTAGGG)n', "TelomericSGL", feature)) %>%
  filter(!is.na(feature))
lineBreakends = beData %>% filter(ResolvedType == 'LINE') %>%
  mutate(
    feature = ifelse(LE == "None", "LineInsertion", NA), 
    feature = ifelse(LE != "None", "LineSource", feature)) %>%
  filter(!is.na(feature))
tiBreakends = beData %>% filter(LocTopType=='TI_ONLY', LnkLen < 1000) %>% 
  mutate(feature = "TiSource") 
foldbackBreakends = beData %>% filter(IsFoldback) %>%
  mutate(feature = ifelse(Orient == 1, "LeftFoldback", "RightFoldback"))

featuredBreakends = bind_rows(dupBreakends, delBreakends) %>% bind_rows(sglBreakends) %>% bind_rows(lineBreakends) %>% bind_rows(tiBreakends) %>% bind_rows(foldbackBreakends)
featuredBreakends = featuredBreakends %>% mutate(feature = factor(feature, levels = featureLevels, ordered = T))

averageGC_1k = read.table("~/hmf/resources/GC_profile.1000bp.cnp", sep = '\t', header = F, stringsAsFactors = F) %>%
  mutate(chromosome = V1, start = V2 + 1, end = start + 1000, gc = V3) %>%
  select(chromosome, start, end, gc) %>%
  filter(gc > -1)

gcRegions = GRanges(averageGC_1k$chromosome, ranges = IRanges(start = averageGC_1k$start, end = averageGC_1k$end)) 
featuredBreakendRegions = GRanges(featuredBreakends$Chr, ranges = IRanges(start = featuredBreakends$Pos, end = featuredBreakends$Pos)) 
ol = as.matrix(findOverlaps(gcRegions, featuredBreakendRegions, type = "any"))
featuredBreakends[ol[, 2], "gc"] = averageGC_1k[ol[, 1], "gc"]

save(featuredBreakends, file = "/Users/jon/hmf/analysis/svPaper/featuredBreakends.RData")
