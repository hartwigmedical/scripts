library(RMySQL)
library(purple)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
patientIdLookups = query_patient_id_lookup(dbPilot)
dbDisconnect(dbPilot)
rm(dbPilot)

sv_set_common_fields<-function(cluster){ 
  cluster %>% mutate( 
  IsLINE = ifelse(LEStart!='None'|LEEnd!='None',T,F),
  IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F),
  IsGenicStart = ifelse(GeneStart!='',T,F),
  IsGenicEnd = ifelse(GeneEnd!='',T,F),
  Length = ifelse(as.character(ChrStart)!=as.character(ChrEnd)|Type=='INS'|ArmEnd!=ArmStart, -1, PosEnd-PosStart),
  #DoubleDupBE = ifelse(DupBEStart=='true'&DupBEEnd=='true',T,F),
  #SingleDupBE = ifelse(DoubleDupBE==0&(DupBEStart=='true'|DupBEEnd=='true'),T,F),
  #TICount = ifelse(LnkTypeStart=='TI',0.5,0)+ifelse(LnkTypeEnd=='TI',0.5,0),
  DBCount = ifelse(DBLenStart>=0,0.5,0)+ifelse(DBLenEnd>=0,0.5,0),
  #IsSglTI = ifelse(LnkTypeStart=='SGL',0.5,0),
  AsmbTICount = ifelse(AsmbMatchStart=='MATCH',0.5,0)+ifelse(AsmbMatchEnd=='MATCH',0.5,0),
  #InferTICount = TICount - AsmbTICount,
  #ShortTICount=ifelse(LnkTypeStart=='TI'&LnkLenStart<=1000,0.5,0)+ifelse(LnkTypeEnd=='TI'&LnkLenEnd<=1000,0.5,0),
  ClusterSize = ifelse(ClusterCount==1,'Single',ifelse(ClusterCount<=4,'Small','Large')),
  IsConsistent = ifelse(Consistency==0,T,F),
  IsChained = (ChainCount>=1),
  IsFoldBack = FoldbackLenStart>=0|FoldbackLenEnd>=0,
  RepeatedChainLink = (cluster$ChainCount>0 & grepl(';',cluster$ChainIndex)),
  IsPolyA = grepl('TTTTTTTTTT',InsertSeq)|grepl('AAAAAAAAAA',InsertSeq),
  IsLowQual = ResolvedType=='LowQual',
  CustomType=case_when(grepl('Chain',ResolvedType) ~ "Chain",grepl('Sgl',ResolvedType) ~ 'PairedSGL', grepl('TI',ResolvedType) ~ 'SyntheticDelDup',TRUE ~ ResolvedType))
}


### DATABASE
query_entire_cohort <- function(dbConnect) {
    query = paste(
    "SELECT p.*",
    " FROM purity p WHERE qcStatus = 'PASS' and purity > 0.2",
    sep = " ")
    return (dbGetQuery(dbConnect, query))
}


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

allPurity = query_entire_cohort(dbProd)
allPurity$patientId <- sapply(allPurity$sampleId, function(x) {purple::sample_to_patient_id(x, patientIdLookups)})
save(allPurity, file = "/Users/jon/hmf/analysis/svEnrichment/allPurity.RData")

hpc = allPurity %>% group_by(patientId) %>% arrange(patientId, -purity) %>% filter(row_number() == 1)
save(hpc, file = "/Users/jon/hmf/analysis/svEnrichment/hpc.RData")

hpcCopyNumbers = purple::query_copy_number(dbProd, hpc)
save(hpcCopyNumbers, file = "/Users/jon/hmf/analysis/svEnrichment/hpcCopyNumbers.RData")

dbDisconnect(dbProd)
rm(dbProd)



########################## BINS
library(BSgenome.Hsapiens.UCSC.hg19)

bins_100k <- tileGenome(seqinfo(Hsapiens), tilewidth=100000, cut.last.tile.in.chrom=TRUE)
bins_1M <- tileGenome(seqinfo(Hsapiens), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
bins_10M <- tileGenome(seqinfo(Hsapiens), tilewidth=10000000, cut.last.tile.in.chrom=TRUE)


########################## Average Copy Number
load(file = "/Users/jon/hmf/analysis/svEnrichment/hpcCopyNumbers.RData")
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

########################## N Count
genome <- BSgenome.Hsapiens.UCSC.hg19
averageCopyNumber_100k = averageCopyNumber(bins_100k, hpcCopyNumbers)
averageCopyNumber_1M = averageCopyNumber(bins_1M, hpcCopyNumbers)
averageCopyNumber_10M = averageCopyNumber(bins_10M, hpcCopyNumbers)

save(averageCopyNumber_100k, averageCopyNumber_1M, averageCopyNumber_10M, file = "/Users/jon/hmf/analysis/svEnrichment/averageCopyNumbers.RData")

mappability <- function(bins) {
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

mappability_100k = mappability(bins_100k)
mappability_1M = mappability(bins_1M)
mappability_10M = mappability(bins_10M)

save(mappability_100k, mappability_1M, mappability_10M, file = "/Users/jon/hmf/analysis/svEnrichment/mappability.RData")


########################## Uniquely Mappable
mappability_150 = read.table(file = "/Users/jon/hmf/resources/unique.mappability.bed", sep = "\t")
colnames(mappability_150) <- c("chromosome", "start0", "end", "id", "mappability") 
mappability_150 = mappability_150 %>% 
  filter(mappability >= 1.00, chromosome != 'MT') %>%
  mutate(start = start0 + 1) %>%
  select(chromosome, start, end, mappability)

#jon = mappability_150 %>% filter(chromosome == 6, end >= 61000000, start <= 62000000)

bins = bins_1M
mappability = jon

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
save(uniquelyMappable_100k, uniquelyMappable_1M, uniquelyMappable_10M, file = "/Users/jon/hmf/analysis/svEnrichment/uniquelyMappable.RData")


jon2 = uniquelyMappable_1M %>% filter(chromosome == 6, binEnd >= 61000000, binStart <= 62000000)



#shortSimpleDels = svs %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4)
#shortSimpleDels = svs %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4)
#rm(svs)
#load(file = "/Users/jon/hmf/analysis/svEnrichment/shortSimpleDels.RData")
#save(shortSimpleDels, file = "/Users/jon/hmf/analysis/svEnrichment/shortSimpleDels.RData")


############################################ PREP ############################################ 
#mappability =  mappability_10M
#averageCopyNumber = averageCopyNumber_10M
#uniquelyMappable = uniquelyMappable_10M
#criteria = end
#type = "end"
load(file = "/Users/jon/hmf/analysis/svEnrichment/hpc.RData")
allSvs = read.csv(file = "/Users/jon/hmf/analysis/svEnrichment/SVA_SVS.csv")
hpcSvs = svs %>% filter(SampleId %in% hpc$sampleId)

rm(list=setdiff(ls(), c("allSvs", "hpcSvs")))


load(file = "/Users/jon/hmf/analysis/svEnrichment/mappability.RData")
load(file = "/Users/jon/hmf/analysis/svEnrichment/averageCopyNumbers.RData")
load(file = "/Users/jon/hmf/analysis/svEnrichment/uniquelyMappable.RData")

enrich <- function(criteria, mappability, averageCopyNumber, uniquelyMappable) {
  result = data.frame()
  
  middle = criteria %>% filter(as.character(ChrStart) == as.character(ChrEnd)) 
  start = criteria %>% filter(as.character(ChrStart) != as.character(ChrEnd)) 
  end = criteria %>% filter(as.character(ChrStart) != as.character(ChrEnd), !Type %in% c('NONE','SGL'), !ChrEnd %in% c('0','MT'))
  
  result = bind_rows(result, enrichInner(middle, mappability, averageCopyNumber, uniquelyMappable, "middle"))
  result = bind_rows(result, enrichInner(start, mappability, averageCopyNumber, uniquelyMappable, "start"))
  result = bind_rows(result, enrichInner(end, mappability, averageCopyNumber, uniquelyMappable, "end"))
  
  return (result)
}


enrichInner <- function(criteria, mappability, averageCopyNumber, uniquelyMappable, type) {
  criteria = criteria %>% mutate(PosMid = (PosEnd + PosStart) / 2)
  
  if (type == "middle") {
    svsRegions = GRanges(criteria$ChrStart, ranges = IRanges(start = criteria$PosMid, end = criteria$PosMid))  
  } else if (type == "start") {
    svsRegions = GRanges(criteria$ChrStart, ranges = IRanges(start = criteria$PosStart, end = criteria$PosStart))  
  } else if (type == "end") {
    svsRegions = GRanges(criteria$ChrEnd, ranges = IRanges(start = criteria$PosEnd, end = criteria$PosEnd))  
  }
  
  binRegions = GRanges(mappability$chromosome, ranges = IRanges(start = mappability$start, end = mappability$end)) 
  uniquelyMappableRegions = GRanges(as.character(uniquelyMappable$chromosome), ranges = IRanges(start = uniquelyMappable$binStart, end = uniquelyMappable$binEnd)) 
  
  ol = as.matrix(findOverlaps(svsRegions, binRegions, type = "within"))
  ol2 = as.matrix(findOverlaps(svsRegions, uniquelyMappableRegions, type = "within"))
  
  criteria[ol[, 1], "N"] <- mappability[ol[, 2], ]$N
  criteria[ol[, 1], "averageCopyNumber"] <- averageCopyNumber[ol[, 2], ]$averageCopyNumber
  criteria$binChromosome <- mappability[ol[, 2], ]$chromosome
  criteria$binStart <- mappability[ol[, 2], ]$start
  criteria$binEnd <- mappability[ol[, 2], ]$end
  
  criteria[ol2[, 1], "uniquelyMappablePercentage"] <- uniquelyMappable[ol2[, 2], ]$uniquelyMappablePercentage
  criteria$uniquelyMappablePercentage <- ifelse(is.na(criteria$uniquelyMappablePercentage), 0, criteria$uniquelyMappablePercentage)
  
  return (criteria)
}

normalise <- function(svs, buckets, meanCopyNumber) {
  normalised_svs = svs %>% 
    mutate(nrows=n(),meanMappable=mean(uniquelyMappablePercentage)) %>%
    group_by(binChromosome, binStart, binEnd, N,nrows, averageCopyNumber, uniquelyMappablePercentage,meanMappable) %>% 
    summarise(
      unnormalisedBucketCount=n(),   
      percentFS=sum(IsFS)/unnormalisedBucketCount,
      avgRep=mean(RepOriginStart)) %>%
    mutate(
      proportionOfN = N / (binEnd - binStart), 
      expectedBucketCount = nrows / buckets * averageCopyNumber / meanCopyNumber * uniquelyMappablePercentage / meanMappable,  
      p=ppois(unnormalisedBucketCount, expectedBucketCount, FALSE),
      q=p.adjust(p,"BH",buckets),
      buckets = buckets)
  
  return (normalised_svs)
}


############################################ Executions ############################################ 
result = data.frame()

#filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart<1e3) %>% mutate(criteria = 'UltraShortDup') %>%
#filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>5e3,PosEnd-PosStart<8e4) %>% mutate(criteria = 'BRCADup') %>%
#filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e5,PosEnd-PosStart<1e6) %>% mutate(criteria = 'LongDup') %>%
#filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e6) %>% mutate(criteria = 'VeryLongDup') %>%
#filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<5e2) %>% mutate(criteria = 'UltraShortDel') %>%
#filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart>5e2,PosEnd-PosStart<1e4) %>% mutate(criteria = 'ShortDel') %>%
#filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart>2e4,PosEnd-PosStart<5e5) %>% mutate(criteria = 'FSLikeDel') %>%
#filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart>5e5) %>% mutate(criteria = 'LongDel') %>%
#filter(ResolvedType=='Line', as.character(ChrStart) == as.character(ChrEnd)) %>% mutate(criteria = 'Line') %>%

#filter(ResolvedType=='Line', as.character(ChrStart) != as.character(ChrEnd)) %>% mutate(criteria = 'Line') %>%
#filter(ResolvedType=='RecipTrans', as.character(ChrStart) != as.character(ChrEnd))%>% mutate(criteria = 'Recip Translocation') %>%

criteria = hpcSvs %>% 
  filter(ResolvedType=='RecipTrans', as.character(ChrStart) != as.character(ChrEnd))%>% mutate(criteria = 'Recip Translocation') %>%
  mutate(IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F))

enriched_100k = enrich(criteria, mappability_100k, averageCopyNumber_100k, uniquelyMappable_100k)
enriched_1M = enrich(criteria, mappability_1M, averageCopyNumber_1M, uniquelyMappable_1M)
enriched_10M = enrich(criteria, mappability_10M, averageCopyNumber_10M, uniquelyMappable_10M)

normalised_100k = normalise(enriched_100k, 28000, mean(averageCopyNumber_100k$averageCopyNumber)) %>% filter(q < 0.01)
normalised_1M = normalise(enriched_1M, 2800, mean(averageCopyNumber_1M$averageCopyNumber)) %>% filter(q < 0.01)
normalised_10M = normalise(enriched_10M, 280, mean(averageCopyNumber_10M$averageCopyNumber)) %>% filter(q < 0.01)

criteriaResult = data.frame()
criteriaResult = bind_rows(criteriaResult, normalised_100k)
criteriaResult = bind_rows(criteriaResult, normalised_1M)
criteriaResult = bind_rows(criteriaResult, normalised_10M)
criteriaResult$criteria <- first(criteria$criteria)

result = bind_rows(result, criteriaResult)
unique(result$criteria)


multiChromosomeResult = result
save(multiChromosomeResult, file = "/Users/jon/hmf/analysis/svEnrichment/multiChromosomeResult.RData")


singleChromosomeResult = result
save(singleChromosomeResult, file = "/Users/jon/hmf/analysis/svEnrichment/singleChromosomeResult.RData")


singleChromosomeResult$method <- "midPoint"
multiChromosomeResult$method <- "breakPoint"


combinedResult = bind_rows(singleChromosomeResult, multiChromosomeResult)
save(combinedResult, file = "/Users/jon/hmf/analysis/svEnrichment/combinedResult.RData")


#save(result, file = "/Users/jon/hmf/analysis/svEnrichment/result.RData")

#translocationResults = result
#save(translocationResults, file = "/Users/jon/hmf/analysis/svEnrichment/translocationResults.RData")

#resultWithoutLineOrRecipTrans = result
#save(resultWithoutLineOrRecipTrans, file = "/Users/jon/hmf/analysis/svEnrichment/resultWithoutLineOrRecipTrans.RData")
#load(file = "/Users/jon/hmf/analysis/svEnrichment/resultWithoutLineOrRecipTrans.RData")
#result = bind_rows(result, resultWithoutLineOrRecipTrans)
#resultWithoutLineOrRecipTrans
############################################ Plot ############################################ 
#ggplot(data=normalised_ssd_1M,aes(x=expectedBucketCount,y=unnormalisedBucketCount)) + geom_point() + xlim(0,50) + ylim(0,50)


#ggplot(data=normalised_ssd_1M) + 
#  geom_point(aes(binStart,expectedBucketCount,colour='expected'),shape="x") +
#  geom_point(aes(binStart,unnormalisedBucketCount,colour='raw'),shape='x') + 
#  facet_wrap(~binChromosome)

#load(file = "/Users/jon/hmf/analysis/svEnrichment/result.RData")
#View(result %>% mutate(bucketWidth=binEnd-binStart+1) %>% filter(unnormalisedBucketCount>5,round(bucketWidth,-3)==bucketWidth,bucketWidth==1e5) %>% group_by(binChromosome,binStart,criteria,isFS=percentFS>0) %>% summarise(q=sum(q)) %>% spread(criteria,q))

#View(translocationResults %>% mutate(bucketWidth=binEnd-binStart+1) %>% filter(unnormalisedBucketCount>5,round(bucketWidth,-3)==bucketWidth,bucketWidth==1e5) %>% 
#       group_by(binChromosome,binStart,criteria,isFS=percentFS>0) %>% summarise(n=sum(unnormalisedBucketCount)) %>% spread(criteria,n))

#View(result %>% mutate(bucketWidth=binEnd-binStart+1) %>% filter(unnormalisedBucketCount>5,round(bucketWidth,-3)==bucketWidth,bucketWidth<=1e6,criteria=='UltraShortDel') %>% 
#       group_by(binChromosome,round(binStart,-6),bucket=ifelse(bucketWidth==1e6,1e6,binStart-round(binStart,-6)),isFS=percentFS>0) %>% summarise(n=sum(unnormalisedBucketCount)) %>% spread(bucket,n))
