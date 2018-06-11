library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)


getSampleIdsStr<-function(samples)
{
  sampleIdsStr = ""

  for(i in 1:nrow(samples))
  {
    sampleId <- samples[i,1]

    if(i > 1)
      sampleIdStr = paste(",'", sampleId, "'", sep="")
    else
      sampleIdStr = paste("'", sampleId, "'", sep="")

    sampleIdsStr = paste(sampleIdsStr, sampleIdStr)
  }

  return (sampleIdsStr)
}


getGeneCopyData<-function(dbConnect,sampleId,gene)
{
  sql = paste("select SampleId, Chromosome, Gene, Start, End, MinCopyNumber, MaxCopyNumber, SomaticRegions, MinRegionStart, MinRegionEnd, ")
  sql = paste(sql, "MinRegionStartSupport, MinRegionEndSupport from geneCopyNumber")
  sql = paste(sql, " where SampleId = '", sampleId, "' and gene = '", gene, "'", sep='')
  return ((dbGetQuery(dbConnect, sql)))
}

getGeneCopyDataMultiple<-function(dbConnect,sampleIds,gene)
{
  sql = paste("select SampleId, Chromosome, Gene, Start, End, MinCopyNumber, MaxCopyNumber, SomaticRegions, MinRegionStart, MinRegionEnd, ")
  sql = paste(sql, "MinRegionStartSupport, MinRegionEndSupport, MinRegionMethod from geneCopyNumber")
  sql = paste(sql, " where gene = '", gene, "' and SampleId in (", sampleIds, ")", sep='')
  return ((dbGetQuery(dbConnect, sql)))
}

####################
### CODE STARTS HERE
####################


# SV data file
#svData = read.csv('~/logs/CLUSTER_V23.csv')
#svData = svData %>% filter(PONCount<2)
# rely on SV data in another script

load("~/data/driversByGene.RData")
View(driversByGene)

load("~/data/hpcDriversByGene.RData")
View(hpcDriversByGene)
driversByGene = hpcDriversByGene

View(driversByGene %>% group_by(type) %>% summarise(Count=n()))
View(driversByGene %>% group_by(gene) %>% summarise(Count=n()))

cancerTypes = driversByGene %>% group_by(cancerType) %>% summarise(Count=n())
View(cancerTypes)

svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))

# Correlation between genes using a probability of co-occurence

getSampleGeneList<-function(dgData)
{
  # create a dataframe with a row for each sampleId and a string list of the genes involved
  sgList = data.frame(matrix(ncol = 3, nrow = 0))
  sgList = setNames(sgList, c("SampleId", "GeneList", "GeneCount"))

  curGeneList = ""
  curSample = ""
  curGeneCount = 0

  for(i in 1:nrow(dgData))
  {
    sgRow = dgData[i,]

    if(curSample != sgRow$sampleId)
    {
      if(nrow(sgList) > 0)
      {
        rowIndex = nrow(sgList)
        sgList[rowIndex,2] = paste(curGeneList, ",", sep='') # append a comma for the grep search below
        sgList[rowIndex,3] = curGeneCount
      }

      curSample = sgRow$sampleId
      curGeneList = paste(",", sgRow$gene, sep='')
      rowIndex = nrow(sgList)+1
      sgList[rowIndex,1] = curSample
      sgList[rowIndex,2] = curGeneList
      sgList[rowIndex,3] = 1
      curGeneCount = 1
    }
    else
    {
      curGeneList = paste(curGeneList, ",", sgRow$gene, sep='')
      curGeneCount = curGeneCount + 1
    }
  }

  return (sgList)
}

calcGeneProbabilities<-function(gsCounts, sgList, sampleCount)
{
  ggPPResults = data.frame(matrix(ncol = 10, nrow = 0))
  ggPPResults = setNames(ggPPResults, c("Gene1", "Gene2", "Gene1SC", "Gene2SC", "BothGenesSC", "BothGenesExpected", "Binomial", "Poisson", "GeneChr1", "GeneChr2"))

  geneCount = n_distinct(gsCounts$gene)

  for(i in 1:nrow(gsCounts))
  {
    geneRow = gsCounts[i,]
    gene1 = geneRow$gene
    scWithGene1 = geneRow$SampleCount
    gene1SamplesPerc = round(scWithGene1/sampleCount,4)

    # if(i > 4) # temp..
    #   break

    # now search each gene and filter on each signature about a significant threshold
    for(j in i+1:nrow(gsCounts))
    {
      if(j > nrow(gsCounts))
        break

      gene2Row = gsCounts[j,]
      gene2 = gene2Row$gene
      scWithGene2 = gene2Row$SampleCount
      gene2SamplesPerc = round(scWithGene2/sampleCount,4)

      # how many samples have both genes
      grep1 = paste(',', gene1, ',', sep='')
      grep2 = paste(',', gene2, ',', sep='')
      scWithGene1And2 = nrow(sgList %>% filter(grepl(grep1, GeneList)&grepl(grep2, GeneList)))

      # expected count of gene 2 in samples with gene 1
      ggExpectedCount = round(sampleCount *gene1SamplesPerc * gene2SamplesPerc,4)

      if(scWithGene1And2 > 0 & scWithGene1And2 > ggExpectedCount)
      {
        # ie prob of not getting 1 less than the actual count or power
        ggBinProb = 1-pbinom(scWithGene1And2-1, scWithGene2, gene1SamplesPerc)
        ggPoisProb = 1 - ppois(scWithGene1And2-1, ggExpectedCount)
      }
      else
      {
        # ie prob ofgetting the actual count or lower
        ggBinProb = pbinom(scWithGene1And2, scWithGene2, gene1SamplesPerc)
        ggPoisProb = ppois(scWithGene1And2, ggExpectedCount)
      }

      # print(paste("Gene1=", gene1, ", Gene2=", gene2, ", gene1SC=", scWithGene1, ", gene2SC=", scWithGene2, ", actualBothSC=", scWithGene1And2, ", expectedBothSC=", ggExpectedCount, ", binProb=", ggBinProb, ", poisProb=", ggPoisProb, sep=''))

      if(ggPoisProb <= 0.2 || ggBinProb <= 0.2)
      {
        if(ggPoisProb < 0.001 || ggBinProb < 0.001)
          print(paste("Gene1=", gene1, ", Gene2=", gene2, ", gene1SC=", scWithGene1, ", gene2SC=", scWithGene2, ", actualBothSC=", scWithGene1And2, ", expectedBothSC=", ggExpectedCount, ", binProb=", round(ggBinProb,6), ", poisProb=", round(ggPoisProb,6), sep=''))

        rowIndex = nrow(ggPPResults)+1
        ggPPResults[rowIndex,1] = gene1
        ggPPResults[rowIndex,2] = gene2
        ggPPResults[rowIndex,3] = scWithGene1
        ggPPResults[rowIndex,4] = scWithGene2
        ggPPResults[rowIndex,5] = scWithGene1And2
        ggPPResults[rowIndex,6] = ggExpectedCount
        ggPPResults[rowIndex,7] = ggBinProb
        ggPPResults[rowIndex,8] = ggPoisProb
        ggPPResults[rowIndex,9] = geneRow$Chromosome
        ggPPResults[rowIndex,10] = gene2Row$Chromosome
      }
    }
  }

  ggPPResults$Count_GT_Exp = ggPPResults$BothGenesSC > ggPPResults$BothGenesExpected

  ggPPResults = ggPPResults %>% arrange(Binomial)

  ggPPResults$Rank = 0
  colIndex = length(colnames(ggPPResults))

  # set ranking values (for use in BH tests)
  for(i in 1:nrow(ggPPResults))
  {
    ggPPResults[i,colIndex] = i
  }

  # ggPPResults$RankPerc = ggPPResults$Rank/geneCountSq
  ggPPResults$GeneCount = geneCount

  return (ggPPResults)
}

driversNonFusions = driversByGene %>% filter(type!='FUSION'&driverLikelihood>=0.5) %>% arrange(sampleId)
nrow(driversNonFusions)

View(driversNonFusions %>% group_by(gene) %>% count())
n_distinct(driversNonFusions$gene)

ggAllProbs = data.frame(matrix(ncol = 12, nrow = 0))

for(i in 1:nrow(cancerTypes)+1)
{
  if(i > nrow(cancerTypes))
  {
    cancerTypeStr = "All"
    dgData = driversNonFusions
  }
  else
  {
    cancerRow = cancerTypes[i,]
    cancerTypeStr = cancerRow$cancerType
    dgData = driversNonFusions %>% filter(cancerType==cancerTypeStr)
  }

  if(TRUE) # cancerTypeStr == 'All'
  #if(cancerTypeStr == 'CNS')
    {
    sampleGeneList = getSampleGeneList(dgData)

    sgList = dgData %>% group_by(gene) %>% summarise(SampleCount=n(), Chromosome=first(chromosome), Type=first())
    sampleCount = n_distinct(dgData$sampleId)

    print(paste(i, ": cancer=", cancerTypeStr, ", sampleCount=", sampleCount, ", geneRecords=", nrow(dgData), sep=''))

    ggProbs = calcGeneProbabilities(sgList, sampleGeneList, sampleCount)

    ctStr = stringi::stri_replace_all_fixed(cancerTypeStr, '/', '')
    ggFilename = paste("~/logs/r_output/ggProb_", ctStr, ".csv", sep='')
    write.csv(ggProbs, ggFilename)

    ggProbs$CancerType = cancerTypeStr
    ggProbs$GeneCount = n_distinct(sgList)

    ggAllProbs = rbind(ggAllProbs, ggProbs)
  }

  # if(i > 2)
  #   break
}

pValueThreshold = 0.05

# geneCountSq = geneCount*geneCount/2
# ggAllProbs$RankPerc = ggAllProbs$BHValue*2
# ggAllProbs$BHValue = ggAllProbs$RankPerc*pValueThreshold

# View(ggProbs)
nrow(ggAllProbs)
View(ggAllProbs)
write.csv(ggAllProbs, "~/logs/r_output/ggProb_all_by_type2.csv", row.names=F, quote=F)
# ggAllProbs = read.csv("~/logs/r_output/ggProb_all_by_type.csv")


# filter for significant results
nrow(ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<pValueThreshold&Binomial<BHValue))
nrow(ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<pValueThreshold&Binomial<BHValue))

nrow(ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<pValueThreshold&Binomial<BHValue))

# merge to get gene type
geneTypes = driversNonFusions %>% group_by(gene) %>% summarise(Type=first(type))
View(geneTypes)

ggAllProbs2 = merge(ggAllProbs, geneTypes, by.x='Gene1', by.y='gene')
ggAllProbs2$Gene1Type = ggAllProbs2$Type
ggAllProbs2 = within(ggAllProbs2, rm(Type))
View(ggAllProbs2)
ggAllProbs2 = merge(ggAllProbs2, geneTypes, by.x='Gene2', by.y='gene')
ggAllProbs2$Gene2Type = ggAllProbs2$Type
ggAllProbs2 = within(ggAllProbs2, rm(Type))
View(ggAllProbs2)

# are there any positively correctly oncogenes and TSGs?
nrow(ggAllProbs2 %>% filter(Count_GT_Exp==T&Binomial<0.01&Gene1Type!=Gene2Type))
View(ggAllProbs2 %>% filter(Count_GT_Exp==T&Binomial<0.001&Gene1Type!=Gene2Type))

# calculate Benjamini-Hochberg level for pan-cancer driver genes
ggPValsLTBH = ggAllProbs %>% filter(CancerType=='All'&Binomial<BHValue) 
View(ggPValsLTBH)
ggMaxPValue_LT_BH = max(ggPValsLTBH$Binomial)
print(ggMaxPValue_LT_BH)
nrow(ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<ggMaxPValue_LT_BH))
nrow(ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<ggMaxPValue_LT_BH))


# Bonferroni correction
bfThreshold = pValueThreshold / (401*401*0.5)
print(bfThreshold)
nrow(ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<bfThreshold))
nrow(ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<bfThreshold))
View(ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<bfThreshold))
View(ggAllProbs %>% filter(GeneChr1!=GeneChr2&Count_GT_Exp==T&Binomial<bfThreshold))
View(ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<bfThreshold))
View(ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<0.01))

# positive relationship
positiveResults = ggAllProbs %>% filter(Count_GT_Exp==T&Binomial<pValueThreshold&Binomial<BHValue)
View(positiveResults)
View(positiveResults %>% filter(GeneChr1!=GeneChr2) %>% arrange(Gene1,Gene2))

# same pairing appearing more than once
View(positiveResults %>% filter(GeneChr1!=GeneChr2) %>% group_by(Gene1,Gene2) %>% summarise(Count=n()) %>% filter(Count > 1))

# negative relationship
negativeResults = ggAllProbs %>% filter(Count_GT_Exp==F&Binomial<pValueThreshold&Binomial<BHValue)
View(negativeResults)
View(negativeResults %>% arrange(Gene1,Gene2))
View(negativeResults %>% arrange(Gene1,Binomial))
View(negativeResults %>% group_by(Gene1,Gene2) %>% summarise(Count=n()) %>% filter(Count > 1))

# backing out distinct gene count by cancer group (but now recorded)
tmpAll = ggAllProbs %>% filter(CancerType=='All')
nrow(tmpAll)
nrow(driversNonFusions %>% group_by(gene) %>% count())

tmpAll$DistinctGenes = sqrt((1/(tmpAll$RankPerc/tmpAll$Rank))*2)
View(tmpAll)
# rankPerc = rank/(geneCount*geneCount*0.5)

geneCancerCounts = driversNonFusions %>% group_by(cancerType,gene) %>% count()
View(geneCancerCounts)
cancerDistinctGeneCounts = geneCancerCounts %>% group_by(cancerType) %>% count()
allRowCount = nrow(cancerDistinctGeneCounts)+1
cancerDistinctGeneCounts[allRowCount,1] = 'All'
cancerDistinctGeneCounts[allRowCount,2] = nrow(driversNonFusions %>% group_by(gene) %>% count())
View(cancerDistinctGeneCounts)

# merge with distinct gene counts and clean-up columns
ggFinalProbs = merge(ggAllProbs %>% filter(Binomial<=0.01), cancerDistinctGeneCounts, by.x='CancerType', by.y='cancerType')
ggFinalProbs = within(ggFinalProbs, rm(Poisson))
ggFinalProbs = within(ggFinalProbs, rm(BHValue))
ggFinalProbs = within(ggFinalProbs, rm(RankPerc))
names(ggFinalProbs)[names(ggFinalProbs) == 'nn'] <- 'GeneCount'
names(ggFinalProbs)[names(ggFinalProbs) == 'Binomial'] <- 'PValue'
View(ggFinalProbs)
ggFinalProbs$QValue = ggFinalProbs$PValue * (ggFinalProbs$GeneCount*ggFinalProbs$GeneCount)/2 / ggFinalProbs$Rank
write.csv(ggFinalProbs, "~/logs/r_output/genePairPValues.csv")

rm(ggAllProbs2)

# calculate a FDR

# BH was rank/hypoTestCount*0.05
tmpAll$HypoTestCount = (tmpAll$DistinctGenes*tmpAll$DistinctGenes)/2
tmpAll$FDR = tmpAll$Binomial*tmpAll$HypoTestCount/tmpAll$Rank
View(tmpAll)

tmpPos = tmpAll %>% filter(Count_GT_Exp==T&GeneChr1!=GeneChr2&FDR<=0.2) %>% arrange(FDR)
View(tmpPos)


# plotting with VennDiagram
# install.packages('VennDiagram')
library(VennDiagram)

topPosResults = head(positiveResults %>% filter(GeneChr1!=GeneChr2) %>% arrange(Binomial),5)
# topPosResults = head(positiveResults %>% filter(CancerType!='All'&GeneChr1!=GeneChr2) %>% arrange(Binomial),5)
View(topPosResults)
genePairInputs = topPosResults
genePairInputs = negAPCResults

for(i in 1:nrow(genePairInputs))
{
  genePairRow = genePairInputs[i,]
  
  grid.newpage()
  draw.pairwise.venn(genePairRow$Gene1SC, genePairRow$Gene2SC, genePairRow$BothGenesSC, category = c(as.character(genePairRow$Gene1), as.character(genePairRow$Gene2)), 
                     lty = rep("blank", 2), 
                     fill = c("light blue", "pink"), 
                     alpha = rep(0.5, 2), 
                     cat.pos = c(0, 0), 
                     cat.dist = rep(0.025, 2))
}



# over-lapping circles

grid.newpage()
draw.pairwise.venn(area1 = 22, area2 = 20, cross.area = 11, category = c("Dog People","Cat People"))

# grid.newpage()
draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), 
                   lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), 
                   cat.pos = c(0, 0), 
                   cat.dist = rep(0.025, 2))

# two non-overlapping circles

grid.newpage()
draw.pairwise.venn(area1 = 22, area2 = 6, cross.area = 0, category = c("Dog People", 
                                                                       "Snake People"), lty = rep("blank", 2), fill = c("light blue", "green"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 180), euler.d = TRUE, sep.dist = 0.03, 
                   rotation.degree = 45)


# bucket all results into p-value buckets
ggAllProbs$PValBucket =  10**round(log(ggAllProbs$Binomial,10))
ggPValStats = ggAllProbs %>% group_by(PValBucket) %>% summarise(Count=n()) %>% arrange(PValBucket)
ggPValStats$BucketPerc = round(ggPValStats$Count / sum(ggPValStats$Count),3)
View(ggPValStats)

# Gene to cancer type relationship
gsCounts = driversNonFusions %>% group_by(gene) %>% summarise(SampleCount=n(), Chromosome=first(chromosome))
View(gsCounts)
geneCount = n_distinct(gsCounts$gene)
print(geneCount) # 401

totalSampleCount = n_distinct(driversNonFusions$sampleId)
print(totalSampleCount) # 2352

gcPPResults = data.frame(matrix(ncol = 8, nrow = 0))
gcPPResults = setNames(gcPPResults, c("CancerType", "Gene", "CancerSC", "GeneSC", "CancerAndGeneSC", "CancerAndGeneExpected", "Binomial", "Poisson"))

for(i in 1:nrow(cancerTypes))
{
  # if(i > 4)
  #   break

  cancerRow = cancerTypes[i,]
  cancerTypeStr = cancerRow$cancerType
  dgData = driversNonFusions %>% filter(cancerType==cancerTypeStr)

  cancerSC = n_distinct(dgData$sampleId)

  print(paste(i, ": cancer=", cancerTypeStr, ", sampleCount=", cancerSC, ", geneRecords=", nrow(dgData), sep=''))

  for(j in 1:nrow(gsCounts))
  {
    # if(j > 4)
    #   break

    geneRow = gsCounts[j,]
    geneName = geneRow$gene
    scWithGene = geneRow$SampleCount
    geneSamplesPerc = round(scWithGene/totalSampleCount,4)

    # how many samples of this cancer type also have this gene
    scWithCancerAndGene = nrow(dgData %>% filter(gene==geneName))

    # expected count of gene 2 in samples with gene 1
    gcExpectedCount = round(cancerSC * geneSamplesPerc,4)

    if(scWithCancerAndGene > 0 & scWithCancerAndGene > gcExpectedCount)
    {
      # ie prob of not getting 1 less than the actual count or lower
      gcBinProb = 1 - pbinom(scWithCancerAndGene-1, cancerSC, geneSamplesPerc)
      gcPoisProb = 1 - ppois(scWithCancerAndGene-1, gcExpectedCount)
    }
    else
    {
      # ie prob ofgetting the actual count or lower
      gcBinProb = pbinom(scWithCancerAndGene, cancerSC, geneSamplesPerc)
      gcPoisProb = ppois(scWithCancerAndGene, gcExpectedCount)
    }

    # print(paste("cancer=", cancerTypeStr, ", gene=", geneName, ", geneSC=", scWithGene, ", actualBothSC=", scWithCancerAndGene, ", expectedBothSC=", gcExpectedCount, ", binProb=", round(gcBinProb,6), ", poisProb=", round(gcPoisProb,6), sep=''))

    if(TRUE)
    # if(gcPoisProb <= 0.2 || gcBinProb <= 0.2)
    {
      if(gcPoisProb < 0.001 || gcBinProb < 0.001)
        print(paste("cancer=", cancerTypeStr, ", gene=", geneName, ", geneSC=", scWithGene, ", actualBothSC=", scWithCancerAndGene, ", expectedBothSC=", gcExpectedCount, ", binProb=", round(gcBinProb,6), ", poisProb=", round(gcPoisProb,6), sep=''))

      rowIndex = nrow(gcPPResults)+1
      gcPPResults[rowIndex,1] = cancerTypeStr
      gcPPResults[rowIndex,2] = geneName
      gcPPResults[rowIndex,3] = cancerSC
      gcPPResults[rowIndex,4] = scWithGene
      gcPPResults[rowIndex,5] = scWithCancerAndGene
      gcPPResults[rowIndex,6] = gcExpectedCount
      gcPPResults[rowIndex,7] = gcBinProb
      gcPPResults[rowIndex,8] = gcPoisProb
    }
  }
}

View(gcPPResults)

gcPPResults$Count_GT_Exp = gcPPResults$CancerAndGeneSC > gcPPResults$CancerAndGeneExpected
gcPPResults = gcPPResults %>% arrange(Binomial)
gcPPResults$Rank = 0
colIndex = length(colnames(gcPPResults))

# set ranking BH values
for(i in 1:nrow(gcPPResults))
{
  gcPPResults[i,colIndex] = i
}

hyposCount = geneCount * nrow(cancerTypes)
print(hyposCount)
gcPPResults$RankPerc = gcPPResults$Rank/hyposCount

pValueThreshold = 0.05
gcPPResults$BHValue = gcPPResults$RankPerc*pValueThreshold

pValuesLTBH = gcPPResults %>% filter(Binomial<BHValue) 
View(pValuesLTBH)
maxPValue_LT_BH = max(pValuesLTBH$Binomial)
print(maxPValue_LT_BH)

# View(ggProbs)
nrow(gcPPResults)
View(gcPPResults)
write.csv(gcPPResults, "~/logs/r_output/geneCancerProbs.csv", row.names=F, quote=F)


# filter for significant results
gcPosResults = gcPPResults %>% filter(Count_GT_Exp==T&Binomial<maxPValue_LT_BH)
View(gcPosResults)
gcNegResults = gcPPResults %>% filter(Count_GT_Exp==F&Binomial<maxPValue_LT_BH)
View(gcNegResults)

# temp: link back up to gene and position
genePositions = driversByGene %>% group_by(gene) %>% summarise(Count=n(), Chr=first(chromosome), Start=first(start), End=first(end))
View(genePositions)

gcPPResultsTmp = merge(gcPPResults, genePositions, by.x='Gene', by.y='gene')
View(gcPPResultsTmp %>% filter(Count_GT_Exp==T&CancerType=='Colon/Rectum') %>% filter(Chr==8) %>% arrange(Binomial))

# Bonferroni correction
bfThreshold = pValueThreshold / hyposCount
print(bfThreshold)
nrow(gcPPResults %>% filter(Binomial<bfThreshold))
View(gcPPResults %>% filter(Binomial<bfThreshold))

# bucket all results into p-value buckets
gcPPResults$PValBucket =  10**round(log(gcPPResults$Binomial,10))
View(gcPPResults)
gcPValStats = gcPPResults %>% group_by(PValBucket) %>% summarise(Count=n()) %>% arrange(PValBucket)
gcPValStats$BucketPerc = round(gcPValStats$Count / sum(gcPValStats$Count),3)
View(gcPValStats)

pvStats = merge(ggPValStats, gcPValStats, by='PValBucket')
View(pvStats)

pValFreqPlot = (ggplot(data = pvStats %>% filter(PValBucket<0.001), aes(x = PValBucket))
               + geom_line(aes(y=Count.x, colour='GenePairs'))
               + geom_line(aes(y=Count.y, colour='CancerGene'))
               + scale_x_log10()
               # + facet_wrap(as.formula(paste("~", facetWrap)))
               + ylab("Frequency") + labs(title = "P-Value Frequency")
)

print(pValFreqPlot)
    


# WGD and gene correlations
load("~/data/highestPurityCohortSummary.RData")
View(highestPurityCohortSummary)

nrow(highestPurityCohortSummary)
wgdSamples = highestPurityCohortSummary %>% filter(WGD==T) %>% select(sampleId)
nrow(wgdSamples)
View(wgdSamples)


gwPPResults = data.frame(matrix(ncol = 8, nrow = 0))
gwPPResults = setNames(gwPPResults, c("CancerType", "Gene", "WGDSC", "GeneSC", "WGDAndGeneSC", "WGDAndGeneExp", "Binomial", "Poisson"))

for(i in 1:nrow(cancerTypes)+1)
{
  if(i > nrow(cancerTypes))
  {
    cancerTypeStr = 'All'
    dgData = driversNonFusions
  }
  else
  {
    cancerRow = cancerTypes[i,]
    cancerTypeStr = cancerRow$cancerType
    dgData = driversNonFusions %>% filter(cancerType==cancerTypeStr)
  }
  
  cancerSC = n_distinct(dgData$sampleId)
  
  gsCounts = dgData %>% group_by(gene) %>% summarise(SampleCount=n())
  geneCount = nrow(gsCounts)

  dgWGDData = dgData %>% filter(sampleId %in% wgdSamples$sampleId)
  
  # WGD rate for this cancer type
  wgdCount = n_distinct(dgWGDData$sampleId)

  print(paste(i, ": cancer=", cancerTypeStr, ", sampleCount=", cancerSC, ", wgdCount=", wgdCount, ", geneCount=", geneCount, sep=''))
  
  for(j in 1:nrow(gsCounts))
  {
    # if(j > 10)
    #   break
    
    geneRow = gsCounts[j,]
    geneName = geneRow$gene
    scWithGene = geneRow$SampleCount
    geneSamplesPerc = round(scWithGene/cancerSC,4) # % of samples with this gene driver
    
    # how many samples of this cancer type also have this gene
    scWithWGDAndGene = nrow(dgWGDData %>% filter(gene==geneName))
    
    # expected count of this gene within samples with WGD
    gwExpectedCount = round(wgdCount * geneSamplesPerc,2)
    
    if(scWithWGDAndGene > 0 & scWithWGDAndGene > gwExpectedCount)
    {
      # ie prob of not getting 1 less than the actual count or lower
      gwBinProb = 1 - pbinom(scWithWGDAndGene-1, wgdCount, geneSamplesPerc)
      gwPoisProb = 1 - ppois(scWithWGDAndGene-1, gwExpectedCount)
    }
    else
    {
      # ie prob of getting the actual count or lower
      gwBinProb = pbinom(scWithWGDAndGene, wgdCount, geneSamplesPerc)
      gwPoisProb = ppois(scWithWGDAndGene, gwExpectedCount)
    }

    # if(TRUE)
    if(gwPoisProb <= 0.2 || gwBinProb <= 0.2)
    {
      if(gwPoisProb < 0.01 || gwBinProb < 0.01)
        print(paste("cancer=", cancerTypeStr, ", gene=", geneName, ", geneSC=", scWithGene, ", actualBothSC=", scWithWGDAndGene, ", expectedBothSC=", gwExpectedCount, ", binProb=", round(gwBinProb,6), ", poisProb=", round(gwPoisProb,6), sep=''))
      
      rowIndex = nrow(gwPPResults)+1
      gwPPResults[rowIndex,1] = cancerTypeStr
      gwPPResults[rowIndex,2] = geneName
      gwPPResults[rowIndex,3] = wgdCount
      gwPPResults[rowIndex,4] = scWithGene
      gwPPResults[rowIndex,5] = scWithWGDAndGene
      gwPPResults[rowIndex,6] = gwExpectedCount
      gwPPResults[rowIndex,7] = gwBinProb
      gwPPResults[rowIndex,8] = gwPoisProb
    }
  }
}

View(gwPPResults)

# validation
nrow(driversNonFusions %>% filter(sampleId %in% wgdSamples$sampleId))

# 481 breast cancer samples
nrow(driversNonFusions %>% filter(cancerType=='Breast') %>% group_by(sampleId) %>% count())

# 275 of which have WGD, so about 60%
nrow(driversNonFusions %>% filter(sampleId %in% wgdSamples$sampleId & cancerType=='Breast') %>% group_by(sampleId) %>% count())

# gene count and percent within breast cancer = 222
nrow(driversNonFusions %>% filter(cancerType=='Breast' & gene=='TP53'))

# so would expect to see 148* 275/481 = 
nrow(driversNonFusions %>% filter(sampleId %in% wgdSamples$sampleId & cancerType=='Breast' & gene=='TP53'))


gwPPResults$Count_GT_Exp = gwPPResults$WGDAndGeneSC > gwPPResults$WGDAndGeneExp
gwPPResults = gwPPResults %>% arrange(Binomial)
gwPPResults$Rank = 0
colIndex = length(colnames(gwPPResults))

# set ranking BH values
for(i in 1:nrow(gwPPResults))
{
  gwPPResults[i,colIndex] = i
}

View(gwPPResults)

write.csv(gwPPResults %>% filter(CancerType!='All'), "~/logs/r_output/geneWGDCorrelations.csv", row.names=F, quote=F)


## Gene and NMF Signatures Correlation
# NMF SV Signatures linked with Driver Genes

# bring in NMF data to combine at sample level- append SV and Signature counts to driver genes
nmfSampleSigData = read.csv("~/logs/r_output/nmf_all_11_sampleSigData.csv")
View(nmfSampleSigData)

sampleSvCounts = nmfSampleSigData %>% group_by(SampleId)%>% summarise(SampleCount=first(SampleSvCount))
driverGenesWithSigs = (merge(driversByGene, sampleSvCounts, by.x="sampleId", by.y="SampleId", all.x=TRUE))
samSigSpread = nmfSampleSigData %>% select(SampleId,SigName,SvCount) %>% spread(SigName,SvCount)
driverGenesWithSigs = merge(driverGenesWithSigs, samSigSpread, by.x="sampleId", by.y="SampleId", all.x=TRUE)
View(driverGenesWithSigs)

# now calculate a Poisson probability of each gene occuring within each signature
load("~/logs/r_output/snvSampleSigData.RData")
nrow(sampleSigData) # 30*2405
View(sampleSigData)
sampleCount = n_distinct(sampleSigData$SampleId)
print(sampleCount)

# first determine a level of significance for each signature (ie SNV count above which a sample would be enriched for this signature)
sigNames = unique(sampleSigData$SigName)
print(sigNames)
nrow(sigNames)

topNPercent = 0.1
# topNIndex = round(topNPercent * sampleCount) 
# print(topNIndex)

# this data gives SNV Count by signature (and a total SNV count per sample)
# from which we can calculate a signature rate for each sample

allGeneSnvSigProbs = data.frame(matrix(ncol = 10, nrow = 0))

for(i in 1:nrow(cancerTypes))
{
  cancerRow = cancerTypes[i,]
  cancerTypeStr = cancerRow$cancerType
  dgData = driversNonFusions %>% filter(cancerType==cancerTypeStr)
  cancerSampleSigData = sampleSigData %>% filter(CancerType==cancerTypeStr)
  cancerSC = n_distinct(dgData$sampleId)
  
  gsCounts = dgData %>% group_by(gene) %>% summarise(SampleCount=n())

  print(paste(i, ": cancerType=", cancerTypeStr, ", sampleCount=", cancerSC, sep=''))
  
  geneSnvSigProbs = data.frame(matrix(ncol = 7, nrow = 0))
  geneSnvSigProbs = setNames(geneSnvSigProbs, c("SigName", "Gene", "SigSC", "GeneSC", "SigAndGeneSC", "SigAndGeneExp", "Binomial"))

  for(sigName in sigNames)
  {
    # ssData = sampleSigData %>% filter(SigName==sigName) %>% arrange(-SnvCount) # samples with highest SNV count for this sig first
    # samplesAboveSigThreshold = head(ssData, topNIndex)
    # ssDataNth = slice(ssData, topNIndex:topNIndex+1)
    # thresSnvCount = ssDataNth$SnvCount
  
    # take samples above 10% by signature contribution
    samplesAboveSigThreshold = cancerSampleSigData %>% filter(SigName==sigName & SigPercent>=topNPercent) # samples with significant contribution from this sig
    scAboveSigThreshold = nrow(samplesAboveSigThreshold)
  
    print(paste("sigName=", sigName, ", SC Above Threshold=", scAboveSigThreshold, sep=''))
    
    dgSamplesWithSig = dgData %>% filter(sampleId %in% samplesAboveSigThreshold$SampleId)
  
    for(j in 1:nrow(gsCounts))
    {
      geneData = gsCounts[j,]
      geneName = geneData$gene
      scWithGene = geneData$SampleCount
      
      geneSamplesPerc = round(scWithGene/cancerSC,4) # % of samples with this gene driver
  
      # how many samples have enriched sig and this gene
      scWithSigAndGene = nrow(dgSamplesWithSig %>% filter(gene==geneName))
      
      # expected count of this gene within samples with WGD
      geneSigExpected = round(scAboveSigThreshold * geneSamplesPerc,2)
        
      if(scWithSigAndGene > 0 & scWithSigAndGene > geneSigExpected)
      {
        # ie prob of not getting 1 less than the actual count or lower
        binomialProb = 1 - pbinom(scWithSigAndGene-1, scAboveSigThreshold, geneSamplesPerc)
      }
      else
      {
        # ie prob of getting the actual count or lower
        binomialProb = pbinom(scWithSigAndGene, scAboveSigThreshold, geneSamplesPerc)
      }
      
      if(binomialProb <= 0.2)
      {
        if(binomialProb < 0.01)
          print(paste("sigName=", sigName, ", gene=", geneName, ", geneSC=", scWithGene, ", actualBothSC=", scWithSigAndGene, ", expectedBothSC=", geneSigExpected, ", binProb=", round(binomialProb,6), sep=''))
        
        rowIndex = nrow(geneSnvSigProbs)+1
        geneSnvSigProbs[rowIndex,1] = sigName
        geneSnvSigProbs[rowIndex,2] = geneName
        geneSnvSigProbs[rowIndex,3] = scAboveSigThreshold
        geneSnvSigProbs[rowIndex,4] = scWithGene
        geneSnvSigProbs[rowIndex,5] = scWithSigAndGene
        geneSnvSigProbs[rowIndex,6] = geneSigExpected
        geneSnvSigProbs[rowIndex,7] = binomialProb
      }
    }
  }
  
  geneSnvSigProbs$Count_GT_Exp = geneSnvSigProbs$SigAndGeneSC > geneSnvSigProbs$SigAndGeneExp
  geneSnvSigProbs = geneSnvSigProbs %>% arrange(Binomial)

  # set ranking values
  geneSnvSigProbs$Rank = 0
  colIndex = length(colnames(geneSnvSigProbs))
  
  for(j in 1:nrow(geneSnvSigProbs))
  {
    geneSnvSigProbs[j,colIndex] = j
  }
  
  numTests = cancerSC

  geneSnvSigProbs$CancerType = cancerTypeStr
  allGeneSnvSigProbs = rbind(allGeneSnvSigProbs, geneSnvSigProbs)
}

length(sigNames)

View(allGeneSnvSigProbs)

View(samplesAboveSigThreshold)
View(geneSnvSigProbs)
nrow(geneSnvSigProbs)

View(driversNonFusions %>% filter(gene=='TP53'&sampleId %in% samplesAboveSigThreshold$SampleId))

View(driversNonFusions)

# validate
View(sampleSigData %>% filter(SigName=="01"&SigPercent>=0.1))

padjust

geneSnvSigProbs$Count_GT_Exp = geneSnvSigProbs$SigAndGeneSC > geneSnvSigProbs$SigAndGeneExp
geneSnvSigProbs = geneSnvSigProbs %>% arrange(Binomial)
geneSnvSigProbs$Rank = 0
colIndex = length(colnames(geneSnvSigProbs))

# set ranking BH values
for(i in 1:nrow(geneSnvSigProbs))
{
  geneSnvSigProbs[i,colIndex] = i
}

View(geneSnvSigProbs)

write.csv(geneSnvSigProbs, "~/logs/r_output/geneSigCor_All.csv", row.names=F, quote=F)




# Gene Info and Prevalence info for reference purposes
driverGeneInfo = (driversByGene %>% group_by(gene,cancerType) 
                 %>% summarise(SampleCount=n(),
                               Type=first(type),
                               Chromosome=first(chromosome),
                               PosStart=first(start),
                               PosEnd=first(end)) 
                 %>% arrange(gene,-SampleCount))

View(driverGeneInfo)
View(driverGeneInfo %>% group_by(gene) %>% count())

driverGeneRefInfo = data.frame(matrix(ncol = ncol(driverGeneInfo), nrow = 0))
driverGeneRefInfo = setNames(driverGeneRefInfo, colnames(driverGeneInfo))
# View(driverGeneRefInfo)
curGene = ""
geneCancerList = ""
curCancerCount = 0
for(i in 1:nrow(driverGeneInfo))
{
  driverGeneRow = driverGeneInfo[i,]
  
  if(is.na(driverGeneRow$gene))
  {
    print(paste("NA gene found at index=", i, ", while curGene=", curGene,  sep=''))
  }
  else
  {
    if(curGene != driverGeneRow$gene)
    {
      curGene = driverGeneRow$gene
      # driverGeneRefInfo = rbind(driverGeneRefInfo, driverGeneRow)
      rowIndex = nrow(driverGeneRefInfo)+1
      driverGeneRefInfo[rowIndex,] = driverGeneRow
      geneCancerList = paste(driverGeneRow$cancerType, "=", driverGeneRow$SampleCount, sep='')
      driverGeneRefInfo[rowIndex,2] = geneCancerList
      curCancerCount = 1
    }
    else
    {
      rowIndex = nrow(driverGeneRefInfo)
      driverGeneRefInfo[rowIndex,3] = driverGeneRefInfo[rowIndex,3] + driverGeneRow$SampleCount # cumulative
  
      # take top 3 cancer contributors
      curCancerCount = curCancerCount + 1
      if(curCancerCount <= 3)
      {
        geneCancerList = paste(geneCancerList, " ", driverGeneRow$cancerType, "=", driverGeneRow$SampleCount, sep='')
        driverGeneRefInfo[rowIndex,2] = geneCancerList
      }
    }
  }
}

driverGeneRefInfo = driverGeneRefInfo %>% arrange(-SampleCount)
View(driverGeneRefInfo)
write.csv(driverGeneRefInfo, "~/logs/r_output/driverGeneInfo.csv", row.names=F, quote=F)


# Retrieve Gene Copy Number data for each driver gene
dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")

# first TSGs
tsgSamples = driversByGene %>% filter(type=='TSG') %>% select(sampleId,gene)

tsgGenes = driversByGene %>% filter(type=='TSG') %>% group_by(gene) %>% summarise(Count=n())
nrow(tsgGenes)
View(tsgGenes)
tsgGeneSamples = driversByGene %>% filter(type=='TSG') %>% group_by(gene,sampleId) %>% summarise(Count=n())
nrow(tsgGeneSamples)

# then oncogenes
oncoGenes = driversByGene %>% filter(type=='ONCO') %>% group_by(gene) %>% summarise(Count=n())
nrow(oncoGenes)
oncoGeneSamples = driversByGene %>% filter(type=='ONCO') %>% group_by(gene,sampleId) %>% summarise(Count=n())
nrow(oncoGeneSamples)

View(oncoGenes)

geneList = oncoGenes
geneSampleList = oncoGeneSamples

View(geneSampleList)

geneRow = geneList[i,1]
View(geneRow)
geneDataSamples = geneSampleList %>% filter(gene==geneRow$gene)
View(geneDataSamples)
geneSamples = data.frame(SampleId = geneDataSamples$sampleId)
View(geneSamples)


gcnDataAll = data.frame()

for(i in 1:nrow(geneList))
{
  geneRow = geneList[i,]
  geneDataSamples = geneSampleList %>% filter(gene==geneRow$gene)
  geneSamples = data.frame(SampleId = geneDataSamples$sampleId)

  print(paste("retreived gene(", geneRow$gene, ") recordCount=", nrow(geneDataSamples), ", sampleCount=", nrow(geneSamples), sep=''))

  if(nrow(geneSamples) > 0)
  {
    samplesStr = getSampleIdsStr(geneSamples)

    gcnData = getGeneCopyDataMultiple(dbProd, samplesStr, geneRow$gene)
    gcnDataAll = rbind(gcnDataAll, gcnData)

    print(paste("retreived gene(", geneRow$gene, ") for sampleCount=", nrow(gcnData), sep=''))
  }
  else
  {
    print(paste("gene(", geneRow$gene, ") has no samples", sep=''))
  }

  # if(i > 5)
  #   break
}

# nrow(gcnDataAll %>% group_by(Gene) %>% summarise(Count=n()))
# nrow(gcnDataAll %>% group_by(SampleId) %>% summarise(Count=n()))

nrow(gcnDataAll)
View(gcnDataAll)
write.csv(gcnDataAll, "~/logs/r_output/onocGenesWithCN.csv", row.names = F, quote = F)

# now combine back with drive info
gcnDataAll$Id = paste(gcnDataAll$SampleId, "_", gcnDataAll$Gene, sep='')

tsgDrivers = driversByGene %>% filter(type=='TSG')
tsgDrivers$Id = paste(tsgDrivers$sampleId, "_", tsgDrivers$gene, sep='')
tsgAllData = (merge(tsgDrivers, gcnDataAll, by.x="Id", by.y="Id", all.x=TRUE))
nrow(tsgAllData)
View(tsgAllData)
tsgAllData = within(tsgAllData, rm(Id))
tsgAllData = within(tsgAllData, rm(SampleId))
tsgAllData = within(tsgAllData, rm(Chromosome))
tsgAllData = within(tsgAllData, rm(Gene))
tsgAllData = within(tsgAllData, rm(Start))
tsgAllData = within(tsgAllData, rm(End))
tsgAllData = within(tsgAllData, rm(StartPosId))
tsgAllData = within(tsgAllData, rm(EndPosId))
write.csv(tsgAllData, "~/logs/r_output/tsg_cn_data.csv", row.names = F, quote = F)

oncoDrivers = driversByGene %>% filter(type=='ONCO')
oncoDrivers$Id = paste(oncoDrivers$sampleId, "_", oncoDrivers$gene, sep='')
oncoAllData = (merge(oncoDrivers, gcnDataAll, by.x="Id", by.y="Id", all.x=TRUE))

oncoAllData = within(oncoAllData, rm(Id))
oncoAllData = within(oncoAllData, rm(SampleId))
oncoAllData = within(oncoAllData, rm(Chromosome))
oncoAllData = within(oncoAllData, rm(Gene))
oncoAllData = within(oncoAllData, rm(Start))
oncoAllData = within(oncoAllData, rm(End))
write.csv(oncoAllData, "~/logs/r_output/onco_cn_data.csv", row.names = F, quote = F)


View(tsgAllData %>% filter(type=='TSG'&driver=='Del') %>% group_by(gene,sampleId) %>% summarise(Count=n()))
View(tsgAllData %>% filter(type=='TSG'&driver=='Del'&MinRegionStartSupport!='NONE'&MinRegionStartSupport!='UNKNOWN'&MinRegionEndSupport!='NONE'&MinRegionEndSupport!='UNKNOWN'))


# DriverGene and SV Data linked
svDgData = read.csv('~/logs/CLUSTER_V24.csv')
View(svDgData)
nrow(svDgData %>% filter(GeneDriverStart=='Del'|GeneDriverEnd=='Del'))
tsgDels = svDgData %>% filter(GeneDriverStart=='Del'|GeneDriverEnd=='Del')
View(tsgDels)

# a single SV covers the gene
tsgDels$GeneBothMatched = tsgDels$GeneDriverStart=='Del' & tsgDels$GeneDriverEnd=='Del'
nrow(tsgDels %>% filter(GeneBothMatched))
View(tsgDels %>% filter(GeneBothMatched&!GeneDoubleup) %>% group_by(Type) %>% summarise(Count=n()))

# a single SV covers 2 genes
tsgDels$GeneDoubleup = tsgDels$GeneBothMatched & as.character(tsgDels$GeneStart)!=as.character(tsgDels$GeneEnd)
nrow(tsgDels %>% filter(GeneDoubleup))

# all others should be either only 1 side of the gene matched or most commonly pairs of SVs
multiSVTsgDel = tsgDels %>% filter(!GeneBothMatched)
multiSVTsgDel$GeneName = ifelse(multiSVTsgDel$GeneDriverStart=='Del',as.character(multiSVTsgDel$GeneStart),as.character(multiSVTsgDel$GeneEnd))

multiSVTsgGroups = (multiSVTsgDel %>% group_by(SampleId,GeneName)
                    %>% summarise(SvCount=n(),
                                  IsDB=sum((GeneDriverStart=='Del'&LnkTypeStart=='DB')|(GeneDriverEnd=='Del'&LnkTypeEnd=='DB')),
                                  DistinctClusters=n_distinct(ClusterCount),
                                  Cluster1=first(ClusterId),
                                  Cluster2=last(ClusterId),
                                  ClusterCount1=first(ClusterCount),
                                  ClusterCount2=last(ClusterCount),
                                  Pos1=first(ifelse(GeneDriverStart=='Del',PosStart,PosEnd)),
                                  Pos2=last(ifelse(GeneDriverStart=='Del',PosStart,PosEnd)),
                                  DBLength=abs(Pos2-Pos1))
                    %>% arrange(SampleId,GeneName))

View(multiSVTsgGroups)

# just focused on 2 SVs in different clusters
View(multiSVTsgGroups %>% filter(SvCount==2&DistinctClusters==2))

multiSVTsgStats = (multiSVTsgGroups %>% group_by(SvCount,DistinctClusters)
                    %>% summarise(Count=n())
                    %>% arrange(SvCount,DistinctClusters))

View(multiSVTsgStats)


tsgDelStats




# svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
# svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',1,0)
# svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
# svData$IsStressed = ifelse((svData$ArmCountStart>=svData$ArmExpStart*2.5&svData$ArmCountStart>=10)|(svData$ArmCountEnd>=svData$ArmExpEnd*2.5&svData$ArmCountEnd>=10),1,0)
# svData$IsDB=ifelse(svData$NearestDBLen>-1&(svData$NearestDBLen<svData$NearestTILen|svData$NearestTILen<30),1,0)
# svData$IsTI=ifelse(svData$NearestTILen>=30&svData$IsDB==0,1,0)
# svData$ClusterNone=ifelse(svData$ClusterCount==1,1,0)
# svData$ClusterSmall=ifelse(svData$ClusterCount>1&svData$ClusterCount<=3,1,0)
# svData$ClusterLarge=ifelse(svData$ClusterCount>3,1,0)
# View(svData)


# Analysis of SV by recurrent positions (to nearest 10M bases)

# associate with cancer type
patientCancerTypes = read.csv('~/data/patient_cancertypes.csv')
svData = (merge(svData, patientCancerTypes, by.x="SampleId", by.y="SampleId", all.x=TRUE))
svData$CancerType = ifelse(is.na(svData$CancerType), 'N/A', paste(svData$CancerType, sep=""))

cancerTypes = svData %>% group_by(CancerType,SampleId) %>% summarise(Count=n())
cancerTypes = cancerTypes %>% group_by(CancerType) %>% summarise(SampleCount=n(), SvCount=sum(Count))
View(cancerTypes)


# first group SVs by distinct chromosomal arm
svStartData = svData
svStartData$Chr = svStartData$ChrStart
svStartData$Arm = svStartData$ArmStart
svStartData$Position = svStartData$PosStart

svEndData = svData
svEndData$Chr = svData$ChrEnd
svEndData$Arm = svData$ArmEnd
svEndData$Position = svData$PosEnd

# merge rows resulting in a set of breakends
beData = rbind(svStartData, svEndData)

beData$PositionBucket = paste(round(beData$Position/1e7)*10, "M", sep='')

beStats = (beData %>% group_by(Chr, PositionBucket)
            %>% summarise(SvCount=n())
            %>% arrange(-SvCount))

View(beStats)

# plot the results
beStats = unite(beStats, "Chr_Pos", Chr, PositionBucket, sep="_")
nrow(beStats)

bePositionPlot = (ggplot(data = beStats[1:100,], aes(x = reorder(Chr_Pos, -SvCount), y = SvCount), fill = Chr_Pos)
                 + geom_bar(stat = "identity", colour = "black", size = 0.2)
                 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                 + ylab("SV Count") + xlab("Chr_Pos")
                 + theme(legend.position="none"))

print(bePositionPlot)


# beCancerStats = (beData %>% group_by(CancerType, Chr, PositionBucket)
#            %>% summarise(SvCount=n())
#            %>% arrange(CancerType, Chr, PositionBucket))
#
# View(beCancerStats)

# plot the results
# DATA OUTPUT TO PDF
pdf(file="~/logs/r_output/sv_by_chr_position.pdf", height = 14, width = 20)
par(mar=c(1,1,1,1))

topNPosCount = 100

beStats = (beData %>% group_by(Chr, PositionBucket)
           %>% summarise(SvCount=n())
           %>% arrange(-SvCount))

beStats = unite(beStats, "Chr_Pos", Chr, PositionBucket, sep="_")

bePositionPlot = (ggplot(data = beStats[1:100,], aes(x = reorder(Chr_Pos, -SvCount), y = SvCount), fill = Chr_Pos)
                  + geom_bar(stat = "identity", colour = "black", size = 0.2)
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                  + ylab("SV Count") + xlab("Chr_Pos") + ggtitle("Position by 10M Buckets")
                  + theme(legend.position="none"))

print(bePositionPlot)

# and plot all by cancer type
plotIndex = 1
chrPosPlotList = list()
for(cancerType in cancerTypes$CancerType)
{
  bePerCancerStats = (beData %>% filter(CancerType==cancerType) %>% group_by(Chr, PositionBucket)
             %>% summarise(SvCount=n())
             %>% arrange(-SvCount))

  bePerCancerStats = unite(bePerCancerStats, "Chr_Pos", Chr, PositionBucket, sep="_")

  title = paste("Position by 10M Buckets: ", cancerType, sep="")

  bePositionPlot = (ggplot(data = bePerCancerStats[1:topNPosCount,], aes(x = reorder(Chr_Pos, -SvCount), y = SvCount), fill = Chr_Pos)
                    + geom_bar(stat = "identity", colour = "black", size = 0.2)
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                    + ylab("SV Count") + xlab("Chr_Pos") + ggtitle(title)
                    + theme(legend.position="none"))

  chrPosPlotList[[plotIndex]] <- bePositionPlot

  if(plotIndex >= 4)
  {
    multiplot(plotlist = chrPosPlotList, cols = 2)
    chrPosPlotList = list()
    plotIndex = 1
  }
  else
  {
    plotIndex = plotIndex + 1
  }
}

if(plotIndex > 1)
{
  # now print all plots for this cancer type
  multiplot(plotlist = chrPosPlotList, cols = 2)
}


# additionally plot positions which are dominated by specific cancer types
topNPosCount = 20
topNPositions = head(beStats,topNPosCount)
View(topNPositions$Chr_Pos)

beCancerStats = (beData %>% group_by(CancerType, Chr, PositionBucket)
           %>% summarise(SvCount=n())
           %>% arrange(CancerType, Chr, PositionBucket))

beCancerStats = unite(beCancerStats, "Chr_Pos", Chr, PositionBucket, sep="_")
beCancerStats = beCancerStats %>% filter(Chr_Pos %in% topNPositions$Chr_Pos)
View(beCancerStats)

title = paste("SV Count by cancer for top ", topNPosCount, " positions", sep='')

beCancerTopNPlot <- (ggplot(beCancerStats, aes(x = reorder(Chr_Pos, -SvCount), y = SvCount, fill = CancerType))
                  + geom_bar(stat = "identity", colour = "black")
                  + labs(x = "", y = "SV count by Chromosomal Position")
                  + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                  + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7))
                  + ggtitle(title))

print(beCancerTopNPlot)


dev.off()








# determine IsStressed using poisson distribution using expected vs actual SV counts per arm
# combinedArmData$ArmExpected = ifelse(combinedArmData$ArmExpected>0,combinedArmData$ArmExpected,0.1)
# combinedArmData$StressedPoissonProb = round(1 - ppois(combinedArmData$ArmCount - 1, combinedArmData$ArmExpected),4)
# combinedArmData$IsStressedOld = ifelse(combinedArmData$ArmCount>=10 & combinedArmData$ArmCount >= 2.5 * combinedArmData$ArmExpected,1,0)
# combinedArmData$IsStressed = ifelse(combinedArmData$StressedPoissonProb <= 0.001 & combinedArmData$ArmCount >= 10,1,0)


# prepare arm stats
armStats = (combinedArmData %>% group_by(SampleId, Chr, Arm)
            %>% summarise(SvCount=n(),
                          MaxCN=round(max((AdjCNStart+AdjCNEnd)*0.5),2),
                          AvgCN=round(sum((AdjCNStart+AdjCNEnd)*0.5)/n(),2),
                          AvgPloidy=round(sum(Ploidy)/n(),2),
                          MaxInvCN=round(max(ifelse(Type=='INV',(AdjCNStart+AdjCNEnd)*0.5,0)),2),
                          InvMinPosStart=min(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosStart,3e8)),
                          InvMaxPosEnd=max(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosEnd,-1)),
                          ClusteredPerc=round(sum(ClusterCount>1)/n(),2),
                          MaxClusterCount=max(ClusterCount),
                          AvgClusterCount=round(sum(ClusterCount)/n(),0),
                          LEPerc=round(sum(IsLINE=='true')/n(),3),
                          FSPerc=round(sum(IsFS=='true')/n(),3),
                          BndCount=sum(Type=='BND'),
                          BndPerc=round(sum(Type=='BND')/n(),2),
                          IsStressed=max(IsStressed))
            %>% arrange(SampleId, Chr, Arm))

View(armStats)







##########################
## bucket and arm analysis

library("svnmf")

chr = "17"
arm = 'P'
chrLen = get_chromosome_length(chr)
centroPos = get_centromere_position(chr, arm)
armLen = get_arm_length(chr,  arm)
View(chrLen)
View(centroPos)
View(get_centromere_position(chr, 'Q'))
View(armLen)
View(get_arm_length(chr,  'Q'))


combinedArmData$ClusterNone=ifelse(combinedArmData$ClusterCount==1,1,0)
combinedArmData$ClusterSmall=ifelse(combinedArmData$ClusterCount>1&combinedArmData$ClusterCount<=3,1,0)
combinedArmData$ClusterLarge=ifelse(combinedArmData$ClusterCount>3,1,0)

armBucketSummary = (combinedArmData %>% group_by(SampleId, Chr, Arm)
                    %>% summarise(LINE=sum(IsLINE==1),
                                  Stressed=sum(IsLINE==0&IsStressed==1),
                                  Del_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterNone==1),
                                  Del_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Del_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Del_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterLarge==1),
                                  Del_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Del_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Del_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Del_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Del_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Del_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Del_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Del_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Del_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Del_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Del_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Del_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Del_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterNone==1),
                                  Del_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Del_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Del_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterLarge==1),
                                  Dup_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterNone==1),
                                  Dup_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Dup_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Dup_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterLarge==1),
                                  Dup_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Dup_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Dup_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Dup_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Dup_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Dup_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Dup_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Dup_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Dup_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Dup_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Dup_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Dup_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Dup_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterNone==1),
                                  Dup_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Dup_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Dup_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterLarge==1),
                                  Inv_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterNone==1),
                                  Inv_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Inv_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Inv_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterLarge==1),
                                  Inv_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Inv_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Inv_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Inv_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Inv_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Inv_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Inv_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Inv_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Inv_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Inv_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Inv_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Inv_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Inv_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterNone==1),
                                  Inv_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Inv_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Inv_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterLarge==1),
                                  Bnd_CN=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterNone==1),
                                  Bnd_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsDB==1),
                                  Bnd_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsTI==1),
                                  Bnd_CL=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterLarge==1)))

View(armBucketSummary)
armBucketSummary = armSummary

# convert to stats for each arm regardless of sample

armBucketStats = (armBucketSummary %>% group_by(Chr,Arm)
                  %>% summarise(LINE=sum(LINE),
                                Stressed=sum(Stressed),
                                Del_LT10K_CN=sum(Del_LT10K_CN),
                                Del_LT10K_CS_DB=sum(Del_LT10K_CS_DB),
                                Del_LT10K_CS_TI=sum(Del_LT10K_CS_TI),
                                Del_LT10K_CL=sum(Del_LT10K_CL),
                                Del_10Kto100K_CN=sum(Del_10Kto100K_CN),
                                Del_10Kto100K_CS_DB=sum(Del_10Kto100K_CS_DB),
                                Del_10Kto100K_CS_TI=sum(Del_10Kto100K_CS_TI),
                                Del_10Kto100K_CL=sum(Del_10Kto100K_CL),
                                Del_100Kto500K_CN=sum(Del_100Kto500K_CN),
                                Del_100Kto500K_CS_DB=sum(Del_100Kto500K_CS_DB),
                                Del_100Kto500K_CS_TI=sum(Del_100Kto500K_CS_TI),
                                Del_100Kto500K_CL=sum(Del_100Kto500K_CL),
                                Del_500Kto5M_CN=sum(Del_500Kto5M_CN),
                                Del_500Kto5M_CS_DB=sum(Del_500Kto5M_CS_DB),
                                Del_500Kto5M_CS_TI=sum(Del_500Kto5M_CS_TI),
                                Del_500Kto5M_CL=sum(Del_500Kto5M_CL),
                                Del_GT5M_CN=sum(Del_GT5M_CN),
                                Del_GT5M_CS_DB=sum(Del_GT5M_CS_DB),
                                Del_GT5M_CS_TI=sum(Del_GT5M_CS_TI),
                                Del_GT5M_CL=sum(Del_GT5M_CL),
                                Dup_LT10K_CN=sum(Dup_LT10K_CN),
                                Dup_LT10K_CS_DB=sum(Dup_LT10K_CS_DB),
                                Dup_LT10K_CS_TI=sum(Dup_LT10K_CS_TI),
                                Dup_LT10K_CL=sum(Dup_LT10K_CL),
                                Dup_10Kto100K_CN=sum(Dup_10Kto100K_CN),
                                Dup_10Kto100K_CS_DB=sum(Dup_10Kto100K_CS_DB),
                                Dup_10Kto100K_CS_TI=sum(Dup_10Kto100K_CS_TI),
                                Dup_10Kto100K_CL=sum(Dup_10Kto100K_CL),
                                Dup_100Kto500K_CN=sum(Dup_100Kto500K_CN),
                                Dup_100Kto500K_CS_DB=sum(Dup_100Kto500K_CS_DB),
                                Dup_100Kto500K_CS_TI=sum(Dup_100Kto500K_CS_TI),
                                Dup_100Kto500K_CL=sum(Dup_100Kto500K_CL),
                                Dup_500Kto5M_CN=sum(Dup_500Kto5M_CN),
                                Dup_500Kto5M_CS_DB=sum(Dup_500Kto5M_CS_DB),
                                Dup_500Kto5M_CS_TI=sum(Dup_500Kto5M_CS_TI),
                                Dup_500Kto5M_CL=sum(Dup_500Kto5M_CL),
                                Dup_GT5M_CN=sum(Dup_GT5M_CN),
                                Dup_GT5M_CS_DB=sum(Dup_GT5M_CS_DB),
                                Dup_GT5M_CS_TI=sum(Dup_GT5M_CS_TI),
                                Dup_GT5M_CL=sum(Dup_GT5M_CL),
                                Inv_LT10K_CN=sum(Inv_LT10K_CN),
                                Inv_LT10K_CS_DB=sum(Inv_LT10K_CS_DB),
                                Inv_LT10K_CS_TI=sum(Inv_LT10K_CS_TI),
                                Inv_LT10K_CL=sum(Inv_LT10K_CL),
                                Inv_10Kto100K_CN=sum(Inv_10Kto100K_CN),
                                Inv_10Kto100K_CS_DB=sum(Inv_10Kto100K_CS_DB),
                                Inv_10Kto100K_CS_TI=sum(Inv_10Kto100K_CS_TI),
                                Inv_10Kto100K_CL=sum(Inv_10Kto100K_CL),
                                Inv_100Kto500K_CN=sum(Inv_100Kto500K_CN),
                                Inv_100Kto500K_CS_DB=sum(Inv_100Kto500K_CS_DB),
                                Inv_100Kto500K_CS_TI=sum(Inv_100Kto500K_CS_TI),
                                Inv_100Kto500K_CL=sum(Inv_100Kto500K_CL),
                                Inv_500Kto5M_CN=sum(Inv_500Kto5M_CN),
                                Inv_500Kto5M_CS_DB=sum(Inv_500Kto5M_CS_DB),
                                Inv_500Kto5M_CS_TI=sum(Inv_500Kto5M_CS_TI),
                                Inv_500Kto5M_CL=sum(Inv_500Kto5M_CL),
                                Inv_GT5M_CN=sum(Inv_GT5M_CN),
                                Inv_GT5M_CS_DB=sum(Inv_GT5M_CS_DB),
                                Inv_GT5M_CS_TI=sum(Inv_GT5M_CS_TI),
                                Inv_GT5M_CL=sum(Inv_GT5M_CL),
                                Bnd_CN=sum(Bnd_CN),
                                Bnd_CS_DB=sum(Bnd_CS_DB),
                                Bnd_CS_TI=sum(Bnd_CS_TI),
                                Bnd_CL=sum(Bnd_CL)))

View(armBucketStats)

bucketNames = tail(colnames(armBucketStats), length(colnames(armBucketStats))-2)
View(bucketNames)

armBucketRelStats = armBucketStats

for(i in 1:nrow(armBucketStats))
{
  chr = armBucketStats[i,1]
  arm = armBucketStats[i,2]

  armLen = get_arm_length(chr, arm)

  for(j in 1:length(bucketNames))
  {
    armBucketRelStats[i,j+2] = armBucketRelStats[i,j+2] / armLen * 1000000
  }
}

armBucketRelStats = unite(armBucketRelStats, "Chr_Arm", Chr, Arm, sep="_")
View(armBucketRelStats)

armBucketPlot = (ggplot(data = armBucketRelStats, aes(x = reorder(Chr_Arm, -Del_LT10K_CN), y = Del_LT10K_CN), fill = Chr_Arm)
                 + geom_bar(stat = "identity", colour = "black", size = 0.2)
                 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                 + ylab("Norm. SV Count") + xlab("Chr_Arm")
                 + theme(legend.position="none"))

print(armBucketPlot)


# DATA OUTPUT TO PDF
pdf(file="~/logs/r_output/bucket_by_arm.pdf", height = 14, width = 20)
par(mar=c(1,1,1,1))

plotIndex = 1
armBucketPlotList = list()
for(bucketName in bucketNames)
{
  selCols = armBucketRelStats %>% select(Chr_Arm, bucketName)
  colnames(selCols) <- c("Chr_Arm", "Count")

  armBucketPlot = (ggplot(data = selCols, aes(x = reorder(Chr_Arm, -Count), y = Count), fill = Chr_Arm)
                   + geom_bar(stat = "identity", colour = "black", size = 0.2)
                   + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                   + ylab(paste("Norm. SV Count: ", bucketName, sep="")) + xlab("Chr_Arm")
                   + theme(legend.position="none"))

  armBucketPlotList[[plotIndex]] <- armBucketPlot
  # print(armBucketPlot)

  if(plotIndex >= 4)
  {
    multiplot(plotlist = armBucketPlotList, cols = 2)
    armBucketPlotList = list()
    plotIndex = 1
  }
  else
  {
    plotIndex = plotIndex + 1
  }
}

if(plotIndex > 1)
{
  # now print all plots for this cancer type
  multiplot(plotlist = armBucketPlotList, cols = 2)
}

dev.off()


svData$HomologyLen = nchar(as.character(svData$Homology))
svData$HomologyLen = ifelse(svData$HomologyLen <= 8, svData$HomologyLen, min(round(svData$HomologyLen/10)*10,100))
View(svData)


# now for homology
homBucketSummary = (svData %>% group_by(SampleId, HomologyLen)
                    %>% summarise(LINE=sum(IsLINE==1),
                                  Stressed=sum(IsLINE==0&IsStressed==1),
                                  Del_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterNone==1),
                                  Del_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Del_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Del_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterLarge==1),
                                  Del_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Del_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Del_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Del_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Del_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Del_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Del_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Del_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Del_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Del_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Del_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Del_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Del_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterNone==1),
                                  Del_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Del_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Del_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterLarge==1),
                                  Dup_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterNone==1),
                                  Dup_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Dup_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Dup_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterLarge==1),
                                  Dup_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Dup_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Dup_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Dup_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Dup_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Dup_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Dup_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Dup_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Dup_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Dup_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Dup_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Dup_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Dup_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterNone==1),
                                  Dup_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Dup_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Dup_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterLarge==1),
                                  Inv_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterNone==1),
                                  Inv_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsDB==1),
                                  Inv_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsTI==1),
                                  Inv_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterLarge==1),
                                  Inv_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterNone==1),
                                  Inv_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                                  Inv_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                                  Inv_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterLarge==1),
                                  Inv_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterNone==1),
                                  Inv_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                                  Inv_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                                  Inv_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterLarge==1),
                                  Inv_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterNone==1),
                                  Inv_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                                  Inv_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                                  Inv_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterLarge==1),
                                  Inv_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterNone==1),
                                  Inv_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsDB==1),
                                  Inv_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsTI==1),
                                  Inv_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterLarge==1),
                                  Bnd_CN=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterNone==1),
                                  Bnd_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsDB==1),
                                  Bnd_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsTI==1),
                                  Bnd_CL=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterLarge==1)))

View(homBucketSummary)

# convert to stats for each arm regardless of sample

homBucketStats = (homBucketSummary %>% group_by(HomologyLen)
                  %>% summarise(LINE=sum(LINE),
                                Stressed=sum(Stressed),
                                Del_LT10K_CN=sum(Del_LT10K_CN),
                                Del_LT10K_CS_DB=sum(Del_LT10K_CS_DB),
                                Del_LT10K_CS_TI=sum(Del_LT10K_CS_TI),
                                Del_LT10K_CL=sum(Del_LT10K_CL),
                                Del_10Kto100K_CN=sum(Del_10Kto100K_CN),
                                Del_10Kto100K_CS_DB=sum(Del_10Kto100K_CS_DB),
                                Del_10Kto100K_CS_TI=sum(Del_10Kto100K_CS_TI),
                                Del_10Kto100K_CL=sum(Del_10Kto100K_CL),
                                Del_100Kto500K_CN=sum(Del_100Kto500K_CN),
                                Del_100Kto500K_CS_DB=sum(Del_100Kto500K_CS_DB),
                                Del_100Kto500K_CS_TI=sum(Del_100Kto500K_CS_TI),
                                Del_100Kto500K_CL=sum(Del_100Kto500K_CL),
                                Del_500Kto5M_CN=sum(Del_500Kto5M_CN),
                                Del_500Kto5M_CS_DB=sum(Del_500Kto5M_CS_DB),
                                Del_500Kto5M_CS_TI=sum(Del_500Kto5M_CS_TI),
                                Del_500Kto5M_CL=sum(Del_500Kto5M_CL),
                                Del_GT5M_CN=sum(Del_GT5M_CN),
                                Del_GT5M_CS_DB=sum(Del_GT5M_CS_DB),
                                Del_GT5M_CS_TI=sum(Del_GT5M_CS_TI),
                                Del_GT5M_CL=sum(Del_GT5M_CL),
                                Dup_LT10K_CN=sum(Dup_LT10K_CN),
                                Dup_LT10K_CS_DB=sum(Dup_LT10K_CS_DB),
                                Dup_LT10K_CS_TI=sum(Dup_LT10K_CS_TI),
                                Dup_LT10K_CL=sum(Dup_LT10K_CL),
                                Dup_10Kto100K_CN=sum(Dup_10Kto100K_CN),
                                Dup_10Kto100K_CS_DB=sum(Dup_10Kto100K_CS_DB),
                                Dup_10Kto100K_CS_TI=sum(Dup_10Kto100K_CS_TI),
                                Dup_10Kto100K_CL=sum(Dup_10Kto100K_CL),
                                Dup_100Kto500K_CN=sum(Dup_100Kto500K_CN),
                                Dup_100Kto500K_CS_DB=sum(Dup_100Kto500K_CS_DB),
                                Dup_100Kto500K_CS_TI=sum(Dup_100Kto500K_CS_TI),
                                Dup_100Kto500K_CL=sum(Dup_100Kto500K_CL),
                                Dup_500Kto5M_CN=sum(Dup_500Kto5M_CN),
                                Dup_500Kto5M_CS_DB=sum(Dup_500Kto5M_CS_DB),
                                Dup_500Kto5M_CS_TI=sum(Dup_500Kto5M_CS_TI),
                                Dup_500Kto5M_CL=sum(Dup_500Kto5M_CL),
                                Dup_GT5M_CN=sum(Dup_GT5M_CN),
                                Dup_GT5M_CS_DB=sum(Dup_GT5M_CS_DB),
                                Dup_GT5M_CS_TI=sum(Dup_GT5M_CS_TI),
                                Dup_GT5M_CL=sum(Dup_GT5M_CL),
                                Inv_LT10K_CN=sum(Inv_LT10K_CN),
                                Inv_LT10K_CS_DB=sum(Inv_LT10K_CS_DB),
                                Inv_LT10K_CS_TI=sum(Inv_LT10K_CS_TI),
                                Inv_LT10K_CL=sum(Inv_LT10K_CL),
                                Inv_10Kto100K_CN=sum(Inv_10Kto100K_CN),
                                Inv_10Kto100K_CS_DB=sum(Inv_10Kto100K_CS_DB),
                                Inv_10Kto100K_CS_TI=sum(Inv_10Kto100K_CS_TI),
                                Inv_10Kto100K_CL=sum(Inv_10Kto100K_CL),
                                Inv_100Kto500K_CN=sum(Inv_100Kto500K_CN),
                                Inv_100Kto500K_CS_DB=sum(Inv_100Kto500K_CS_DB),
                                Inv_100Kto500K_CS_TI=sum(Inv_100Kto500K_CS_TI),
                                Inv_100Kto500K_CL=sum(Inv_100Kto500K_CL),
                                Inv_500Kto5M_CN=sum(Inv_500Kto5M_CN),
                                Inv_500Kto5M_CS_DB=sum(Inv_500Kto5M_CS_DB),
                                Inv_500Kto5M_CS_TI=sum(Inv_500Kto5M_CS_TI),
                                Inv_500Kto5M_CL=sum(Inv_500Kto5M_CL),
                                Inv_GT5M_CN=sum(Inv_GT5M_CN),
                                Inv_GT5M_CS_DB=sum(Inv_GT5M_CS_DB),
                                Inv_GT5M_CS_TI=sum(Inv_GT5M_CS_TI),
                                Inv_GT5M_CL=sum(Inv_GT5M_CL),
                                Bnd_CN=sum(Bnd_CN),
                                Bnd_CS_DB=sum(Bnd_CS_DB),
                                Bnd_CS_TI=sum(Bnd_CS_TI),
                                Bnd_CL=sum(Bnd_CL)))

View(homBucketStats)

write.csv(homBucketStats, "~/logs/r_output/hom_bucket_stats.csv")




# position bucketing to nearest 1M bases by cancer type
