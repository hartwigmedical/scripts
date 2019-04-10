library(dplyr)
library(tidyr)
library(stringi)
library("pracma") # for ceil function

# plotting
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(scater)

# Load sample counts in matrix form
snvMatrixData = read.csv('~/data/sigs/snv_prod_matrix_data_20190409.csv')
nrow(snvMatrixData)
ncol(snvMatrixData)
View(snvMatrixData[,1:10])
View(snvMatrixData[,1146:1147])

prodSamples = colnames(snvMatrixData)
View(prodSamples)

# re-create SNV buckets
snvBuckets = create_empty_signature()
snvBuckets = snvBuckets %>% mutate(Bucket=paste(type,context,sep='_'))
snvBuckets = snvBuckets %>% select(Bucket)
View(snvBuckets)

# create sample counts from matrix data
snvSampleCounts = matrix_to_sample_counts(snvMatrixData,snvBuckets)
View(snvSampleCounts)
View(snvSampleCounts %>% group_by(SampleId) %>% count()) # 4311 samples (at 09/04/2019), without QC

# limit to relevant samples in HPC
load('~/data/r_data/highestPurityCohortSummary.RData')
nrow(highestPurityCohortSummary)

sampleCancerTypes = read.csv('~/data/sigs/sample_cancer_types.csv')
sampleCancerTypes[is.na(sampleCancerTypes)] <- 'Unknown'

snvDpSampleCounts = snvSampleCounts %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
View(snvDpSampleCounts %>% group_by(SampleId) %>% count()) # 2370, not 2406 for some reason?
snvDpSampleTotals = snvDpSampleCounts %>% group_by(SampleId) %>% count()
nrow(snvDpSampleTotals)

# no missing DP samples 
missingDpSamples = highestPurityCohortSummary %>% filter(!(sampleId %in% snvDpSampleTotals$SampleId))
View(missingDpSamples)
View(missingDpSamples$sampleId)

# write out DP sample matrix data
snvDpMatrixData = sample_counts_to_matrix(snvDpSampleCounts)
View(snvDpMatrixData[,1:10])
ncol(snvDpMatrixData)

write.csv(snvDpMatrixData, '~/data/sigs/snv_dp_matrix_data.csv', row.names = F, quote = F)


# load in results from SigAnalyser NMF fit with PCAWG sigs
snvDpPcawgContribs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_fit_nmf_contribs.csv', stringsAsFactors = F))
snvDpPcawgSigs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_fit_nmf_sigs.csv', stringsAsFactors = F))

snvPcawgSigsNamesStr = get_pcawg_named_sigs()
length(snvPcawgSigsNamesStr)

View(snvDpPcawgSigs)
View(snvDpPcawgContribs[,1:20])
#evaluate_nmf_data<-function(runType, runId, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
#                            sigNamesNamed, plotByCancerType = T, viewResults = F, bgSigCount = 0, printAllPlots = T)
  
evaluate_nmf_data("SNV", "dp_pcawg_fit", snvDpPcawgSigs, snvDpPcawgContribs, snvDpMatrixData, snvDpSampleCounts,
                  sampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, F, F, 0, F)





# only skin and lung
View(sampleCancerTypes %>% group_by(CancerType) %>% count())

skinAndLungSamples = sampleCancerTypes %>% filter(CancerType=='Skin'|CancerType=='Lung')
nrow(skinAndLungSamples)
View(skinAndLungSamples %>% group_by(CancerType) %>% count())
snvSkinLungSampleCounts = snvDpSampleCounts %>% filter(SampleId %in% skinAndLungSamples$SampleId)
nrow(snvSkinLungSampleCounts %>% group_by(SampleId) %>% count())

# write out DP sample matrix data
snvSkinLungMatrixData = sample_counts_to_matrix(snvSkinLungSampleCounts)
View(snvSkinLungMatrixData[,1:10])
ncol(snvSkinLungMatrixData)

write.csv(snvSkinLungMatrixData, '~/data/sigs/snv_dp_skin_lung_matrix_data.csv', row.names = F, quote = F)


snvSkinLungPcawgContribs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_skin_lung_fit_nmf_contribs.csv', stringsAsFactors = F))
snvSkinLungPcawgSigs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_skin_lung_fit_nmf_sigs.csv', stringsAsFactors = F))

View(snvSkinLungPcawgContribs[,1:20])
ncol(snvSkinLungPcawgContribs)

# snvSkinLungSampleCounts = matrix_to_sample_counts(snvSkinLungMatrixData, snvBuckets)

evaluate_nmf_data("SNV", "skin_lung_pcawg_fit", snvSkinLungPcawgSigs, snvSkinLungPcawgContribs, snvSkinLungMatrixData, snvSkinLungSampleCounts,
                  sampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, T, F)


specificSample = snvDpSampleCounts %>% filter(SampleId=='CPCT02010772T')
View(specificSample)
specificSampleMatrixData = sample_counts_to_matrix(specificSample)
View(specificSampleMatrixData)
write.csv(specificSampleMatrixData, '~/data/sigs/snv_dp_syd895_matrix_data.csv', row.names = F, quote = F)



## plotting

# install_github("taiyun/corrplot")
library(corrplot)

snvDpSampleSigData[is.na(snvDpSampleSigData)] <- "Unknown"


cancerSampleCounts = sampleCancerTypes %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
rowIndex = data.frame(as.numeric(as.character(rownames(cancerSampleCounts))))
colnames(rowIndex) <- c("CancerIndex")
cancerSampleCounts = cbind(rowIndex, cancerSampleCounts)
View(cancerSampleCounts)

cancerSampleSigSummary = (snvDpSampleSigData %>% filter(SigName!='Excess'&SigName!='Unalloc') %>% group_by(CancerType,SigName) 
                          %>% summarise(CancerSigSampleCount=sum(Count>0),
                                        CancerSigSampleMutLoad=round(sum(Count),0)))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSampleCounts,by='CancerType',all.x=T)
View(cancerSampleSigSummary)

cancerTypesBySig = cancerSampleSigSummary %>% group_by(SigName) %>% summarise(SigCancerTypeCount=sum(CancerSigSampleCount>0))
rowIndex = data.frame(as.numeric(as.character(rownames(cancerTypesBySig))))
colnames(rowIndex) <- c("SigIndex")
cancerTypesBySig = cbind(rowIndex, cancerTypesBySig)
View(cancerTypesBySig)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerTypesBySig,by='SigName',all.x=T)
View(cancerSampleSigSummary)

# work out median mutational load per sig and cancer type for samples with the signature
cancerSigMedianCounts = data.frame(ncol(3))
rowIndex = 1
for(sig in cancerTypesBySig$SigName)
{
  for(cancerType in cancerSampleCounts$CancerType)
  {
    sampleData = snvDpSampleSigData %>% filter(SigName==sig&CancerType==cancerType&Count>0)
    cancerSigMedianCounts[rowIndex,1] = sig
    cancerSigMedianCounts[rowIndex,2] = cancerType
    cancerSigMedianCounts[rowIndex,3] = round(median(sampleData$Count),0)
    rowIndex = rowIndex + 1
  }
}

colnames(cancerSigMedianCounts) = c('SigName','CancerType','MedianSampleCount')
View(cancerSigMedianCounts)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSigMedianCounts,by=c('SigName','CancerType'),all.x=T)
View(cancerSampleSigSummary)

cancerSampleSigSummary = cancerSampleSigSummary %>% mutate(SamplesWithSig=round(CancerSigSampleCount/CancerSampleCount,3),
                                                           MedianLoadPerMb=ifelse(MedianSampleCount>0,round(MedianSampleCount/3e3,3),0))
                                                           

cancerSampleSigSummary = (cancerSampleSigSummary %>% mutate(MedianLoadBucket=round(0.5**round(log(MedianLoadPerMb,0.5)),2),
                                                            MedianLoadGroup=ifelse(MedianLoadBucket<0.075,0.05,
                                                                                   ifelse(MedianLoadBucket<0.017,'0.1',
                                                                                   ifelse(MedianLoadBucket<0.37,'0.25',
                                                                                   ifelse(MedianLoadBucket<0.75,'0.5',
                                                                                   ifelse(MedianLoadBucket<1.75,'1',
                                                                                   ifelse(MedianLoadBucket<3.75,'2.5',
                                                                                   ifelse(MedianLoadBucket<7.5,'5',
                                                                                   ifelse(MedianLoadBucket<17.5,'10','25'))))))))))
                                                            

View(cancerSampleSigSummary)
# View(cancerSampleSigSummary %>% group_by(SigName,CancerType,MedianLoadBucket) %>% summarise(SamplesWithSig=first(SamplesWithSig)) %>% spread(MedianLoadBucket,SamplesWithSig))

plot = (ggplot(cancerSampleSigSummary %>% filter(SamplesWithSig>0), aes(x=CancerType, y=reorder(SigName, -SigIndex))) 
        + geom_point(aes(col=MedianLoadGroup, size=SamplesWithSig)) 
        # + geom_smooth(method="loess", se=F) 
        + theme_gray()
        # + theme(axis.text.x = element_text(data=cancerSampleCounts$CancerType))
        + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
        + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7))
        +  labs(subtitle="Signature Allocation by Cancer Types", y="Signature", x="CancerType"))

print(plot)

# aes(x = reorder(SampleId, -SampleCount), y = Count, fill = SigName)


data("midwest", package = "ggplot2")

View(midwest)




# plotting:
# SamplesWIthSig = 0 -> 1 (ie percentage of samples in cancer type with signature)
# 





# Positive correlations are displayed in blue and negative correlations in red color. 
# Color intensity and the size of the circle are proportional to the correlation coefficients.

# Matching PCAWG's plot:
# Size of circle: Proportion of samples in this cancer type with this signature - use SamplesWithSig
# Colour of circle: Median mutations/Mb due to signature (among samples with the signature)  - use MedianLoadPerMb


testM = matrix(c(1,4,9,5,5,4,5,7,6),nrow=3, n)
testMCorr <- cor(testM)
print(testMCorr)

# M <- cor(mtcars)
corrplot(testMCorr, method = "circle")

testM = matrix(c(0.25,0,0.75,-0.1,-1.0,0.8,1.0,-0.1,-0.5),nrow=3, ncol=3)

varNames = c('Sig1','Sig3','Sig3','Sig1','Sig3','Sig3')
cols = c('Sig1','Sig3','Sig3')
rows = c('Sig1','Sig3','Sig3')
print(varNames)
corrplot(testM, method = "circle")
corrplot(testM, method = "circle", col=cols)
corrplot(testM, varNames, method = "circle", abs=FALSE, n.col.legend=7)





library(ggplot2)
theme_set(theme_bw())  # pre-set the bw theme.
data("midwest", package = "ggplot2")
# midwest <- read.csv("http://goo.gl/G1K41K")  # bkup data source

# Scatterplot
gg <- ggplot(midwest, aes(x=area, y=poptotal)) + 
  geom_point(aes(col=state, size=popdensity)) + 
  geom_smooth(method="loess", se=F) + 
  xlim(c(0, 0.1)) + 
  ylim(c(0, 500000)) + 
  labs(subtitle="Area Vs Population", 
       y="Population", 
       x="Area", 
       title="Scatterplot", 
       caption = "Source: midwest")




#### DEBUG
runType = "SNV"
runId = "sking_lung_fit"
signatures = snvDpPcawgSigs
contribution = snvDpPcawgContribs
matrixData = snvDpMatrixData
summaryCounts = snvDpSampleCounts
bucketNames = snvBuckets
sigNamesNamed = snvPcawgSigsNamesStr
plotByCancerType = F
bgSigCount = 0
printAllPlots = F
sigColours = c()
sigInfo = data.frame()
sigAllocs = data.frame()

# data prep
sigCount = nrow(contribution)
hasBackgroundSigs = (bgSigCount > 0)

print(paste("evaluating run: type=", runType, ", id=", runId, ", sigCount=", sigCount, sep=''))

sampleNames = colnames(contribution)
length(sampleNames)
origSampleCounts = summaryCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))

sigNamesUnamed = get_signame_list(sigCount, F)
print(sigNamesUnamed)
colnames(signatures) = sigNamesUnamed

sigNamesCombined = cbind(sigNamesUnamed, sigNamesNamed)
colnames(sigNamesCombined) <- c("Signature", "SigName")
View(sigNamesCombined)

bucketIndex = data.frame(as.numeric(as.character(rownames(bucketNames))))
colnames(bucketIndex) <- c("BucketIndex")
bucketNamesIndexed = cbind(bucketNames, bucketIndex)
bucketNamesIndexed$BucketIndex = bucketNamesIndexed$BucketIndex-1
# View(bucketNamesIndexed)

sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
sigBucketData = merge(sigBucketData,bucketNamesIndexed,by="Bucket",all.x=T)
sigBucketData = merge(sigBucketData,sigNamesCombined,by="Signature",all.x=T)

sigBucketStats = get_sig_bucket_stats(sigBucketData)
sigBucketTopN = get_top_buckets(sigBucketData)
bucketSummaryData = get_bucket_stats(sigBucketData)
#View(sigBucketStats)
#View(sigBucketTopN)

sampleBucketData = get_sample_bucket_data(matrixData, origSampleCounts, bucketNames)

# A: routine for bucket analyser providing Unalloc and Excess counts
sigNamesCombined = rbind(sigNamesCombined, c(sigCount+1,"Unalloc"))
sigNamesCombined = rbind(sigNamesCombined, c(sigCount+2,"Excess"))
View(sigNamesCombined)
  
# ensure an entry in every sig or every sample
sampleSigFullSet = merge(sampleNames, sigNamesCombined)
colnames(sampleSigFullSet) = c("SampleId", "Signature", "SigName")

sampleSigData = merge(sigAllocs, sampleSigFullSet, by=c("SampleId", "Signature"), all=T)
sampleSigData[is.na(sampleSigData)] = 0
sampleSigData$PercBucket = round(sampleSigData$SigPercent/0.1)*0.1
sampleSigData = within(sampleSigData, rm(BgId))

# B: routine for when BucketAnalyser not used (eg NMF)
sampleSigData = get_sig_data(signatures, contribution, sigNamesNamed, sampleNames)
  
# calculate and factor in residuals
residuals = calc_contrib_sample_residuals(contribution, signatures, matrixData, bucketNames)
View(residuals)

sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)
View(sampleSigData)


# get cancer type and SV Count
sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)
# sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))
sampleSigData = merge(sampleSigData, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=T)
sampleSigDataNoResiduals = sampleSigData %>% filter(SigName!="Unalloc"&SigName!="Excess")

snvDpSampleSigData = sampleSigData
View(snvDpSampleSigData)

residualsSummary = (snvDpSampleSigData %>% group_by(SampleId,CancerType) 
     %>% summarise(SampleTotal=first(SampleCount), 
                   ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
     %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4)))

View(residualsSummary)

# Residuals Report

print(ggplot(data = residualsSummary, aes(ResidualTotal))
      + stat_ecdf()
      + scale_x_log10()
      # + facet_wrap(~ClusterDesc)
      )

print(ggplot(data = residualsSummary, aes(ResidualPerc))
      + stat_ecdf())

dev.off()












