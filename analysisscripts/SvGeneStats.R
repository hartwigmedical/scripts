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

View(driversByGene %>% group_by(type) %>% summarise(Count=n()))
View(driversByGene %>% group_by(gene) %>% summarise(Count=n()))

# experiment with CDKN2A first
View(driversByGene %>% filter(grepl('CDKN2A',gene)))
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))

View(driversByGene %>% filter(type=='TSG') %>% group_by(gene) %>% summarise(Count=n()))
View(driversByGene %>% filter(gene=='CDKN2A'))

View(driversByGene %>% filter(sampleId=='CPCT02010035T'))
View(tsgAllData %>% filter(sampleId=='CPCT02010035T'))


# NMF SV Signatures linked with Driver Genes

# bring in NMF data to combine at sample level- append SV and Signature counts to driver genes
nmfSampleSigData = read.csv("~/logs/r_output/nmf_all_11_sampleSigData.csv")
View(nmfSampleSigData)

sampleSvCounts = nmfSampleSigData %>% group_by(SampleId)%>% summarise(SampleCount=first(SampleSvCount))
driverGenesWithSigs = (merge(driversByGene, sampleSvCounts, by.x="sampleId", by.y="SampleId", all.x=TRUE))
samSigSpread = nmfSampleSigData %>% select(SampleId,SigName,SvCount) %>% spread(SigName,SvCount)
driverGenesWithSigs = (merge(driverGenesWithSigs, samSigSpread, by.x="sampleId", by.y="SampleId", all.x=TRUE))
View(driverGenesWithSigs)

# now calculate a Poisson probability of each gene occuring within each signature
sampleCount = n_distinct(driverGenesWithSigs$sampleId)

# for each driver gene, work out prevalence in cohort
distinctGenes = driverGenesWithSigs %>% group_by(gene) %>% summarise(Count=n())
nrow(distinctGenes)
View(distinctGenes)

sigNames = nmfSampleSigData %>% group_by(SigName) %>% summarise(Count=n())
View(sigNames)


topNValue = 0.02
sgPPResults = data.frame(matrix(ncol = 9, nrow = 0))
sgPPResults = setNames(sgPPResults, c("Gene", "SigName", "GeneSamplePerc", "ThresSvCount", "SampleCountAboveThres", "GeneCountAboveThres", "ExpSvCount", "Binomial", "Poisson"))

for(i in 1:nrow(distinctGenes))
{
  geneRow = distinctGenes[i,]

  # if(geneRow$gene == 'CDK12')
  # {
    geneSamplePerc = round(geneRow$Count/sampleCount, 4)
    # print(paste("Gene=", geneRow$gene, ", count=", geneRow$Count, ", perc=", geneSamplePerc, sep=''))

    # now search each gene and filter on each signature about a significant threshold
    for(j in 1:nrow(sigNames))
    {
      sigRow = sigNames[j,]

      sigColumn = paste("`", sigRow$SigName, "`", sep='')
      sigArrange = paste("-`", sigRow$SigName, "`", sep='')
      sigCond1 = paste("`", sigRow$SigName, "`", ">0", sep='')
      thresSvCount = tail(head(driverGenesWithSigs %>% arrange_(sigArrange), nrow(driverGenesWithSigs %>% filter_(sigCond1))*topNValue),1) %>% select_(sigColumn)
      sigCond2 = paste("`", sigRow$SigName, "`", ">=", thresSvCount, sep='')

      # how many samples have elevated levels of this signature
      samplesAboveThres = nrow(driverGenesWithSigs %>% filter_(sigCond2))

      # how many of those relate to this specific gene
      geneCountAboveThres = nrow(driverGenesWithSigs %>% filter(gene==geneRow$gene) %>% filter_(sigCond2))

      sigPercent = samplesAboveThres / sampleCount
      geneExpPercent = round(geneSamplePerc*sigPercent,6)

      print(paste("Gene=", geneRow$gene, ", Sig=", sigRow$SigName, ", samplesAboveThres=", samplesAboveThres, ", geneCountAboveThres=", geneCountAboveThres, ", median=", thresSvCount, sep=''))

      binomialProb = round(1-pbinom(geneCountAboveThres, sampleCount, geneExpPercent), 6)
      print(paste("Gene=", geneRow$gene, ", Sig=", sigRow$SigName, ", geneExpPercent=", geneExpPercent, ", binomialProb=", binomialProb, sep=''))

      gsPoissonProb = 0
      if(samplesAboveThres > 0 & geneCountAboveThres > 0)
      {
        # now work out significance
        geneExpCount = round(geneSamplePerc * samplesAboveThres,0)

        logData = F
        if(binomialProb < 0.01)
        {
          logData = T
        }

        if(geneCountAboveThres > geneExpCount & geneExpCount > 0)
        {
          gsPoissonProb = round(1 - ppois(geneCountAboveThres, geneExpCount), 8)

          logData = T
          if(gsPoissonProb < 0.001)
          {
            print(paste("Gene=", geneRow$gene, ", Sig=", sigRow$SigName, ", samplesAboveThres=", samplesAboveThres, ", sgGeneCount=", geneCountAboveThres, ", expCount=", geneExpCount, ", PP=", gsPoissonProb, sep=''))
          }
        }
        else if(geneCountAboveThres > 0 & geneExpCount == 0)
        {
          logData = T
          print(paste("Gene=", geneRow$gene, ", Sig=", sigRow$SigName, ", samplesAboveThres=", samplesAboveThres, ", sgGeneCount=", geneCountAboveThres, " vs zero expCount", sep=''))
        }
        else
        {
          print(paste("Gene=", geneRow$gene, ", Sig=", sigRow$SigName, ", samplesAboveThres=", samplesAboveThres, ", sgGeneCount=", geneCountAboveThres, ", expCount=", geneExpCount, sep=''))
        }

        if(logData)
        {
          rowIndex = nrow(sgPPResults)+1
          sgPPResults[rowIndex,1] = geneRow$gene
          sgPPResults[rowIndex,2] = stringi::stri_replace_all_fixed(sigColumn, "`", "")
          sgPPResults[rowIndex,3] = geneSamplePerc
          sgPPResults[rowIndex,4] = thresSvCount
          sgPPResults[rowIndex,5] = samplesAboveThres
          sgPPResults[rowIndex,6] = geneCountAboveThres
          sgPPResults[rowIndex,7] = geneExpCount
          sgPPResults[rowIndex,8] = binomialProb
          sgPPResults[rowIndex,9] = gsPoissonProb
        }
      }
    }

  #   if(geneRow$gene == 'CDK12')
  #      break
  # }

}

View(sgPPResults)
write.csv(sgPPResults, "~/logs/r_output/sgPPResults_0.02.csv")


View(driverGenesWithSigs %>% filter(gene=='CDK12'))
View(driverGenesWithSigs %>% filter(`02_LongDUP` > 40,gene %in% c('CDK12','CCNE1','ERBB2','MYC')) %>% group_by(gene,sampleId) %>% count() %>% spread(gene,n))


# Retrieve Gene Copy Number data for each driver gene
dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")

# first TSGs
tsgSamples = driversByGene %>% filter(type=='TSG') %>% select(sampleId,gene)

tsgGenes = driversByGene %>% filter(type=='TSG') %>% group_by(gene) %>% summarise(Count=n())
nrow(tsgGenes)
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

# try linking to SV data using sample, chr and position
nrow(svData)

svData$StartPosId = paste(svData$SampleId, "_", svData$ChrStart, "_", svData$PosStart, sep='')
svData$EndPosId = paste(svData$SampleId, "_", svData$ChrEnd, "_", svData$PosEnd, sep='')
View(svData)
tsgAllData$StartPosId = paste(tsgAllData$SampleId, "_", tsgAllData$Chromosome, "_", tsgAllData$MinRegionStart, sep='')
tsgAllData$EndPosId = paste(tsgAllData$SampleId, "_", tsgAllData$Chromosome, "_", (tsgAllData$MinRegionEnd+1), sep='')

# join svData to TSG
tsgSvData = (merge(tsgAllData, svData %>% filter(IsSpan==0&OrientStart==1), by.x="StartPosId", by.y="StartPosId", all.x=TRUE))
View(tsgSvData)
tsgSvData = within(tsgSvData, rm(SampleId.x))
nrow(tsgSvData %>% filter(!is.na(ClusterId)))
tsgSvData = (merge(tsgSvData, svData %>% filter(IsSpan==0&OrientEnd==-1), by.x="EndPosId.x", by.y="EndPosId", all.x=TRUE))
nrow(tsgSvData %>% filter(!is.na(ClusterId.x)|!is.na(ClusterId.y)))

# remove redundant columns
tsgSvData = within(tsgSvData, rm(StartPosId.x))
tsgSvData = within(tsgSvData, rm(StartPosId.y))
tsgSvData = within(tsgSvData, rm(EndPosId.x))
tsgSvData = within(tsgSvData, rm(EndPosId.y))
# tsgSvData = within(tsgSvData, rm(Id.x))
tsgSvData = within(tsgSvData, rm(SampleId.y))
View(tsgSvData)

tsgSvData$MissedSv = (tsgSvData$MinRegionStartSupport!="NONE"&is.na(tsgSvData$Id.x)) | (tsgSvData$MinRegionEndSupport!="NONE"&is.na(tsgSvData$Id.y))

tsgSvData$HasSV = !is.na(tsgSvData$Id.x) | !is.na(tsgSvData$Id.y)
tsgSvData$HasBothSVs = !is.na(tsgSvData$Id.x) & !is.na(tsgSvData$Id.y)
tsgSvData$HasSingleSV = tsgSvData$HasBothSVs & (tsgSvData$Id.x==tsgSvData$Id.y)
tsgSvData$HasDiffSVs = tsgSvData$HasBothSVs & (tsgSvData$Id.x!=tsgSvData$Id.y)
tsgSvData$HasLinkedSVs = (tsgSvData$HasBothSVs & tsgSvData$HasDiffSVs
                          & (tsgSvData$LnkSvStart.x==tsgSvData$Id.y|tsgSvData$LnkSvStart.y==tsgSvData$Id.x|tsgSvData$LnkSvEnd.x==tsgSvData$Id.y|tsgSvData$LnkSvEnd.y==tsgSvData$Id.x))

tsgSvData$SvLink = (ifelse(tsgSvData$HasLinkedSVs,'Linked_SVs',ifelse(tsgSvData$HasDiffSVs,'Diff_SVs',
                      ifelse(tsgSvData$HasSingleSV,'Same_SV',ifelse(tsgSvData$HasSV,'SingleMatch','NoMatch')))))


View(tsgSvData %>% filter(MinRegionStartSupport!="NONE"&is.na(Id.x)))

nrow(tsgSvData %>% filter(!HasSV))
nrow(tsgSvData %>% filter(HasSV))
nrow(tsgSvData %>% filter(HasBothSVs))
nrow(tsgSvData %>% filter(HasSingleSV))
nrow(tsgSvData %>% filter(HasDiffSVs))
nrow(tsgSvData %>% filter(HasLinkedSVs))

# summary stats for TSGs linked to SVs
tsgSvStats = (tsgSvData %>% group_by(SvLink)
              %>% summarise(Count=n(),
                            MissedSVCount=sum(MissedSv),
                            DelCount=sum(Type=='DEL'),
                            DupCount=sum(Type=='DUP'),
                            InvCount=sum(Type=='INV'),
                            BndCount=sum(Type=='BND'),
                            BndCount=sum(Type=='INS'))
              %>% arrange(SvLink))

View(tsgSvStats)



tsgSvData$RegionLength = tsgSvData$MinRegionEnd-tsgSvData$MinRegionStart+1

# names(tsgSvData)[names(tsgSvData) == 'Id.y'] <- 'SvId'

nrow(tsgSvData %>% filter(!is.na(Id.x)))
nrow(tsgSvData %>% filter(!is.na(Id.y)))
nrow(tsgSvData %>% filter(!is.na(Id.x)|!is.na(Id.y)))
nrow(tsgSvData %>% filter(!is.na(Id.x)&!is.na(Id.y)))
nrow(tsgSvData %>% filter(!is.na(Id.x)&!is.na(Id.y)&Id.x!=Id.y))

svPairTsgData = tsgSvData %>% filter(!is.na(Id.x)&!is.na(Id.y)&Id.x!=Id.y)

# linked pairs matching the deleted region
View(svPairTsgData %>% filter(LnkSvStart.x==Id.y|LnkSvStart.y==Id.x|LnkSvEnd.x==Id.y|LnkSvEnd.y==Id.x))

# unlinked pairs matching the deleted region
View(svPairTsgData %>% filter(LnkSvStart.x!=Id.y&LnkSvStart.y!=Id.x&LnkSvEnd.x!=Id.y&LnkSvEnd.y!=Id.x))

View(svPairTsgData %>% filter(LnkLenStart.x==RegionLength|LnkLenEnd.x==RegionLength|LnkLenStart.y==RegionLength|LnkLenEnd.y==RegionLength))

write.csv(tsgSvData, "~/logs/r_output/tsgSvData.csv")


svGcnData = (merge(svData, tsgAllData, by.x="StartPosId", by.y="StartPosId", all.x=TRUE))
View(svGcnData %>% filter(SampleId.x=='CPCT02010022T'&ChrStart==9))
View(svGcnData)
nrow(svGcnData %>% filter(!is.na(Id.x)))
nrow(svGcnData %>% filter(!is.na(Id.y)))
nrow(svGcnData %>% filter(is.na(Id.y)))
nrow(svGcnData %>% filter(type=='TSG'))

svGcnData2 = (merge(svData, tsgAllData, by.x="EndPosId", by.y="EndPosId", all.x=TRUE))
nrow(svGcnData2 %>% filter(!is.na(Id.y)))
nrow(svGcnData2 %>% filter(is.na(Id.y)))
nrow(svGcnData2 %>% filter(type=='TSG'))


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
