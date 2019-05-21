library(RMySQL)
library(purple)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)

########################## BINS
library(BSgenome.Hsapiens.UCSC.hg19)
bins_100k <- tileGenome(seqinfo(Hsapiens), tilewidth=100000, cut.last.tile.in.chrom=TRUE)
bins_1M <- tileGenome(seqinfo(Hsapiens), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
bins_10M <- tileGenome(seqinfo(Hsapiens), tilewidth=10000000, cut.last.tile.in.chrom=TRUE)

########################## Average Mappability
mappability <- function(bins, genome) {
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

mappability_100k = mappability(bins_100k, BSgenome.Hsapiens.UCSC.hg19)
mappability_1M = mappability(bins_1M, BSgenome.Hsapiens.UCSC.hg19)
mappability_10M = mappability(bins_10M, BSgenome.Hsapiens.UCSC.hg19)

save(mappability_100k, mappability_1M, mappability_10M, file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")


########################## Average Copy Number 
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

load(file = "/Users/jon/hmf/analysis/cohort/hpcCopyNumbers.RData")
averageCopyNumber_100k = averageCopyNumber(bins_100k, hpcCopyNumbers)
averageCopyNumber_1M = averageCopyNumber(bins_1M, hpcCopyNumbers)
averageCopyNumber_10M = averageCopyNumber(bins_10M, hpcCopyNumbers)
save(averageCopyNumber_100k, averageCopyNumber_1M, averageCopyNumber_10M, file = "/Users/jon/hmf/analysis/svPaper/averageCopyNumbers.RData")

########################## Average Replication
replication = read.table(file = "~/hmf/analysis/svPaper/heli_rep_origins.bed", sep = '\t', header = F, stringsAsFactors = F) %>% 
  mutate(chromosome = substring(V1, 4), start = V2 + 1, end = V3, replication = V4) %>%
  select(chromosome, start, end, replication)

averageReplication <- function(bins, replication) {
  result = averageCopyNumber(bins, replication %>% select(chromosome, start, end, copyNumber = replication))
  result = result %>% 
    mutate(averageReplication = averageCopyNumber / 100) %>%
    select(chromosome, start, end, averageReplication)
  return (result)
}


averageReplication_100k = averageReplication(bins_100k, replication)
averageReplication_1M = averageReplication(bins_1M, replication)
averageReplication_10M = averageReplication(bins_10M, replication)
save(averageReplication_100k, averageReplication_1M, averageReplication_10M, file = "/Users/jon/hmf/analysis/svPaper/averageReplication.RData")


########################## Uniquely Mappable
mappability_150 = read.table(file = "/Users/jon/hmf/resources/unique.mappability.bed", sep = "\t")
colnames(mappability_150) <- c("chromosome", "start0", "end", "id", "mappability") 
mappability_150 = mappability_150 %>% 
  filter(mappability >= 1.00, chromosome != 'MT') %>%
  mutate(start = start0 + 1) %>%
  select(chromosome, start, end, mappability)

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
save(uniquelyMappable_100k, uniquelyMappable_1M, uniquelyMappable_10M, file = "/Users/jon/hmf/analysis/svPaper/uniquelyMappable.RData")


############################################ HPC ############################################ 
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
allowedBins = mappability_100k %>% mutate(mappable = (100000 - N) / 100000) %>% filter(mappable >= 0.8)
binRegions = GRanges(allowedBins$chromosome, ranges = IRanges(start = allowedBins$start, end = allowedBins$end)) 
regions_1M = GRanges(mappability_1M$chromosome, ranges = IRanges(start = mappability_1M$start, end = mappability_1M$end)) 
regions_10M = GRanges(mappability_10M$chromosome, ranges = IRanges(start = mappability_10M$start, end = mappability_10M$end)) 

bins_100k = nrow(allowedBins)
bins_1M = length(unique(as.matrix(findOverlaps(regions_1M, binRegions, type = "any"))[, 1]))
bins_10M = length(unique(as.matrix(findOverlaps(regions_10M, binRegions, type = "any"))[, 1]))
save(bins_100k, bins_1M, bins_10M, file = "/Users/jon/hmf/analysis/svPaper/binSizes.RData")

load(file = "/Users/jon/hmf/analysis/cohort/cohort.RData")
allSvs = read.csv(file = "/Users/jon/hmf/analysis/svPaper/SVA_SVS.csv")

hpcDelsDups =  allSvs %>% 
  filter(SampleId %in% highestPurityCohort$sampleId, ResolvedType=='SimpleSV', Type %in% c('DEL', 'DUP')) %>%
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
  ) 


svsRegions = GRanges(hpcDelsDups$ChrStart, ranges = IRanges(start = (hpcDelsDups$PosStart + hpcDelsDups$PosEnd) / 2, end = (hpcDelsDups$PosStart + hpcDelsDups$PosEnd ) / 2))  
ol = as.matrix(findOverlaps(svsRegions, binRegions, type = "any"))
hpcDelsDups = hpcDelsDups[ol[, 1], ]

save(hpcDelsDups, file = "/Users/jon/hmf/analysis/svPaper/hpcDelsDups.RData")
rm(list=setdiff(ls(), c("hpcDelsDups")))



############################################ PREP ############################################ 
load(file = "/Users/jon/hmf/analysis/svPaper/binSizes.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/hpcDelsDups.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/averageCopyNumbers.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/uniquelyMappable.RData")
load(file = "~/hmf/analysis/svPaper/replicationFactor.RData")
load(file = "/Users/jon/hmf/analysis/svPaper/averageReplication.RData")

replication_bucket <- function(replication) {
  cut(replication, breaks = c(0, seq(0.062, 0.78, 0.02), 1))
}

enrich <- function(criteria, mappability, averageCopyNumber, uniquelyMappable, averageReplication) {
  result = data.frame()
  
  middle = criteria %>% filter(as.character(ChrStart) == as.character(ChrEnd), subtype != 'Line') 
  start = criteria %>% filter(as.character(ChrStart) != as.character(ChrEnd) || subtype == 'Line') 
  end = criteria %>% filter(as.character(ChrStart) != as.character(ChrEnd) || subtype == 'Line', !Type %in% c('NONE','SGL'), !ChrEnd %in% c('0','MT'))
  
  result = bind_rows(result, enrichInner(middle, mappability, averageCopyNumber, uniquelyMappable, averageReplication, "middle"))
  result = bind_rows(result, enrichInner(start, mappability, averageCopyNumber, uniquelyMappable, averageReplication, "start"))
  result = bind_rows(result, enrichInner(end, mappability, averageCopyNumber, uniquelyMappable, averageReplication, "end"))
  
  return (result)
}

enrichInner <- function(criteria, mappability, averageCopyNumber, uniquelyMappable, averageReplication, type) {
  if (nrow(criteria) == 0) {
    return (criteria)
  }
  
  criteria = criteria %>% mutate(PosMid = (PosEnd + PosStart) / 2, breakend = type)
  
  if (type == "middle") {
    svsRegions = GRanges(criteria$ChrStart, ranges = IRanges(start = criteria$PosMid, end = criteria$PosMid))  
  } else if (type == "start") {
    svsRegions = GRanges(criteria$ChrStart, ranges = IRanges(start = criteria$PosStart, end = criteria$PosStart))  
  } else if (type == "end") {
    svsRegions = GRanges(criteria$ChrEnd, ranges = IRanges(start = criteria$PosEnd, end = criteria$PosEnd))  
  }
  
  binRegions = GRanges(mappability$chromosome, ranges = IRanges(start = mappability$start, end = mappability$end)) 
  uniquelyMappableRegions = GRanges(as.character(uniquelyMappable$chromosome), ranges = IRanges(start = uniquelyMappable$binStart, end = uniquelyMappable$binEnd)) 
  replicationRegions = GRanges(as.character(averageReplication$chromosome), ranges = IRanges(start = averageReplication$start, end = averageReplication$end)) 
  
  ol = as.matrix(findOverlaps(svsRegions, binRegions, type = "within"))
  ol2 = as.matrix(findOverlaps(svsRegions, uniquelyMappableRegions, type = "within"))
  ol3 = as.matrix(findOverlaps(svsRegions, replicationRegions, type = "within"))
  
  criteria[ol[, 1], "N"] <- mappability[ol[, 2], ]$N
  criteria[ol[, 1], "averageCopyNumber"] <- averageCopyNumber[ol[, 2], ]$averageCopyNumber
  criteria$binChromosome <- mappability[ol[, 2], ]$chromosome
  criteria$binStart <- mappability[ol[, 2], ]$start
  criteria$binEnd <- mappability[ol[, 2], ]$end
  
  criteria[ol2[, 1], "uniquelyMappablePercentage"] <- uniquelyMappable[ol2[, 2], ]$uniquelyMappablePercentage
  criteria$uniquelyMappablePercentage <- ifelse(is.na(criteria$uniquelyMappablePercentage), 0, criteria$uniquelyMappablePercentage)
  
  criteria[ol3[, 1], "averageReplication"] <- averageReplication[ol3[, 2], ]$averageReplication
  criteria$averageReplication <- ifelse(is.na(criteria$averageReplication), 0.00001, criteria$averageReplication)
  criteria$averageReplication <- ifelse(criteria$averageReplication < 0.00001, 0.00001, criteria$averageReplication)
  criteria$replicationBucket = replication_bucket(criteria$averageReplication)
  
  criteria = criteria %>% left_join(replicationFactor, by = c("subtype", "replicationBucket"))
  
  return (criteria)
}


normalise <- function(svs, buckets, meanCopyNumber, meanMappable) {

  normalised_svs = svs %>% 
    ungroup() %>%
    mutate(nrows=n()) %>%
    group_by(binChromosome, binStart, binEnd, N, nrows, averageCopyNumber, uniquelyMappablePercentage, averageReplication, replicationFactor) %>% 
    summarise(
      unnormalisedBucketCount=n(),   
      percentFS=sum(IsFS)/unnormalisedBucketCount) %>%
    mutate(
      #proportionOfN = N / (binEnd - binStart), 
      expectedBucketCount = nrows / buckets * replicationFactor * averageCopyNumber / meanCopyNumber,  
      #expectedBucketCount = nrows / buckets * replicationFactor,  
      p=ppois(unnormalisedBucketCount, expectedBucketCount, FALSE),
      q=p.adjust(p,"BH",buckets),
      buckets = buckets)
  
  sum(normalised_svs$expectedBucketCount)
  sum(normalised_svs$unnormalisedBucketCount)
  
  
  return (normalised_svs)
}

enrich_and_normalise <- function(criteria, qFilter = 0.01) {
  
  enriched_100k = enrich(criteria, mappability_100k, averageCopyNumber_100k, uniquelyMappable_100k, averageReplication_100k)
  enriched_1M = enrich(criteria, mappability_1M, averageCopyNumber_1M, uniquelyMappable_1M, averageReplication_1M)
  enriched_10M = enrich(criteria, mappability_10M, averageCopyNumber_10M, uniquelyMappable_10M, averageReplication_10M)
  
  normalised_100k = normalise(enriched_100k, bins_100k, mean(averageCopyNumber_100k$averageCopyNumber), mean(uniquelyMappable_100k$uniquelyMappablePercentage)) %>% filter(q < qFilter)
  normalised_1M = normalise(enriched_1M, bins_1M, mean(averageCopyNumber_1M$averageCopyNumber), mean(uniquelyMappable_1M$uniquelyMappablePercentage)) %>% filter(q < qFilter)
  normalised_10M = normalise(enriched_10M, bins_10M, mean(averageCopyNumber_10M$averageCopyNumber), mean(uniquelyMappable_10M$uniquelyMappablePercentage)) %>% filter(q < qFilter)
  
  criteriaResult = data.frame()
  criteriaResult = bind_rows(criteriaResult, normalised_100k)
  criteriaResult = bind_rows(criteriaResult, normalised_1M)
  criteriaResult = bind_rows(criteriaResult, normalised_10M)
  criteriaResult$subtype <- dplyr::first(criteria$subtype)
  
  return (criteriaResult)
}


############################################ Executions ############################################ 
#filter(ResolvedType=='Line') %>% mutate(criteria = 'Line') %>%
#filter(ResolvedType=='RecipTrans', as.character(ChrStart) != as.character(ChrEnd))%>% mutate(criteria = 'Recip Translocation') %>%


df = hpcDelsDups
result = data.frame()
for (selectedSubtype in unique(df$subtype)) {
  cat ("Processing", selectedSubtype, "\n")
  criteria = df %>% 
    filter(subtype == selectedSubtype) %>%
    mutate(IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F))
  result = bind_rows(result, enrich_and_normalise(criteria, 10000000))
}

svEnrichment = result
save(svEnrichment, file = "/Users/jon/hmf/analysis/svPaper/svEnrichment.RData")


############################################ Add Genes ############################################ 
load( file = "/Users/jon/hmf/analysis/svPaper/svEnrichment.RData")
svEnrichment$gene <- NA

load(file = "~/hmf/analysis/svPaper/canonicalTranscripts.RData")

transcriptRegions = GRanges(canonicalTranscripts$chromosome, ranges = IRanges(start = canonicalTranscripts$geneStart, end = canonicalTranscripts$geneEnd))  
enrichmentRegions = GRanges(svEnrichment$binChromosome, ranges = IRanges(start = svEnrichment$binStart, end = svEnrichment$binEnd))  
ol = as.matrix(findOverlaps(enrichmentRegions, transcriptRegions, type = "any"))
genes = canonicalTranscripts[ol[, 2], ]
genes$index <- ol[, 1]
genes = genes %>% mutate(length = geneEnd - geneStart) %>% group_by(index) %>% arrange(-length) %>% filter(row_number() == 1)

svEnrichment[genes$index, "gene"] <- genes$gene
svEnrichment[genes$index, "geneLength"] <- genes$length
save(svEnrichment, file = "/Users/jon/hmf/analysis/svPaper/svEnrichment.RData")

View(svEnrichment %>% filter(binEnd-binStart+1==1000000) %>% group_by(gene,binChromosome,binStart,subtype,round(averageReplication,1),isFS=percentFS>0) %>% summarise(q=mean(q)) %>% spread(subtype,q) )
View(svEnrichment %>% filter(binEnd-binStart+1==100000) %>% group_by(gene,binChromosome,binStart,subtype,round(averageReplication,1),isFS=percentFS>0) %>% summarise(q=mean(q)) %>% spread(subtype,q) )

############################################ 1M Genes ############################################ 
load(file = "/Users/jon/hmf/analysis/svPaper/mappability.RData")
transcriptRegions = GRanges(canonicalTranscripts$chromosome, ranges = IRanges(start = canonicalTranscripts$geneStart, end = canonicalTranscripts$geneEnd))  
enrichmentRegions = GRanges(mappability_1M$chromosome, ranges = IRanges(start = mappability_1M$start, end = mappability_1M$end))  
ol = as.matrix(findOverlaps(enrichmentRegions, transcriptRegions, type = "any"))
genes = canonicalTranscripts[ol[, 2], ]
genes$index <- ol[, 1]
genes = genes %>% mutate(length = geneEnd - geneStart) %>% group_by(index) %>% arrange(-length) %>% filter(row_number() == 1)
mappability_1M[genes$index, "gene"] <- genes$gene
mappability_1M[genes$index, "geneLength"] <- genes$length

mappability_1M_with_genes = mappability_1M
save(mappability_1M_with_genes, file = "~/hmf/analysis/svPaper/mappability_1M_with_genes.RData")

#selectedSubtype = "MediumDup"
#jon = enrich(criteria, mappability_100k, averageCopyNumber_100k, uniquelyMappable_100k, averageReplication_100k)
#jon = enrichInner(criteria, mappability_100k, averageCopyNumber_100k, uniquelyMappable_100k, averageReplication_100k, "end")
#jon$replicationBucket = replication_bucket(jon $>$ )
#mappability = mappability_100k;
#averageCopyNumber = averageCopyNumber_100k
#uniquelyMappable = uniquelyMappable_100k 
#averageReplication = averageReplication_100k;
#type ="middle"
#unique(result$criteria)


View(svEnrichment %>% group_by(binChromosome,binStart,binWidth=binEnd-binStart+1,isFS=percentFS>0,subtype) %>% summarise(n=sum(q)) %>% spread(subtype,n))

load(file = "/Users/jon/hmf/analysis/svPaper/combinedResult.RData")
View(combinedResult %>% group_by(binChromosome,binStart,buckets,isFS=percentFS>0,criteria) %>% summarise(n=sum(q)) %>% spread(criteria,n))

View(svEnrichment %>% group_by(subtype) %>% summarise(sum(unnormalisedBucketCount),sum(expectedBucketCount)))
View(svEnrichment %>% filter(binEnd-binStart+1==1000000) %>% group_by(binChromosome,binStart,subtype,round(averageReplication,1),isFS=percentFS>0) %>% summarise(q=mean(q)) %>% spread(subtype,q) )


load(file = "/Users/jon/hmf/analysis/svPaper/svEnrichment.RData")
new = svEnrichment

load(file = "/Users/jon/hmf/analysis/svPaper/old/svEnrichment.RData")
old = svEnrichment

##### CIRCOS

jon = combinedResult %>% group_by(criteria, buckets) %>% count()

copyNumberTrack = combinedResult %>% filter(buckets == 28000, criteria == "LongDup") %>% mutate(chromosome = paste0("hs", binChromosome)) %>% select(chromosome, binStart, binEnd, averageCopyNumber)
write.table(copyNumberTrack, file = "/Users/jon/hmf/analysis/svPaper/svEnrichment/copyNumberTrack.circos", quote = F, row.names = F, col.names = F, sep = "\t")

unnormalisedBucketCountTrack = combinedResult %>% filter(buckets == 28000, criteria == "LongDup") %>% mutate(chromosome = paste0("hs", binChromosome)) %>% select(chromosome, binStart, binEnd, unnormalisedBucketCount)
write.table(unnormalisedBucketCountTrack, file = "/Users/jon/hmf/analysis/svPaper/svEnrichment/unnormalisedBucketCountTrack.circos", quote = F, row.names = F, col.names = F, sep = "\t")

expectedBucketCountTrack = combinedResult %>% filter(buckets == 28000, criteria == "LongDup") %>% mutate(chromosome = paste0("hs", binChromosome)) %>% select(chromosome, binStart, binEnd, expectedBucketCount)
write.table(expectedBucketCountTrack, file = "/Users/jon/hmf/analysis/svPaper/svEnrichment/expectedBucketCountTrack.circos", quote = F, row.names = F, col.names = F, sep = "\t")

View(svEnrichment %>% filter(buckets==28000) %>% group_by(binChromosome,binStart,criteria,round(avgRep,1),isFS=percentFS>0) %>% summarise(q=mean(q)) %>% spread(criteria,q) )

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
