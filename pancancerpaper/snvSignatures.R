library(devtools)
library(purple); # for multiplot
library(data.table)
library(dplyr)
library(tidyr)
library(stringi)
library(MutationalPatterns)
library(ggplot2)


snvDpSampleSigData = read.csv('~/data/sigs/dp/snvDpSampleSigData.csv')
snvDpSampleCounts = read.csv('~/data/sigs/dp/snvDpSampleCounts.csv')
nmfFitContributions = read.csv('~/data/sigs/dp/snvDpCosmicFitContributions.csv')
snvDpMatrixData = read.csv('~/data/sigs/dp/snvDpMatrixData.csv')
cosmicSigs = read.csv('~/data/sigs/dp/cosmicSigs.csv')

load('~/data/r_data/highestPurityCohortSummary.RData')
sampleCancerTypes = highestPurityCohortSummary %>% select(SampleId=sampleId,CancerType=cancerType)

snvDpSampleSigData = merge(snvDpSampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)

# merge in sample totals and cancer type
sampleTotals = snvDpSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=round(sum(Count),0))
snvDpSampleSigData = merge(snvDpSampleSigData,sampleTotals,by.x="SampleId",by.y="SampleId",all.x=T)

########################################
# 1. Signature Prevalence by Cancer Type

# filter for samples with a sig contribution above 5% and 300 count
sigPercentThreshold = 0.05
sigCountThreshold = 500
highMLSampleSigData = snvDpSampleSigData %>% filter(SigPercent>=sigPercentThreshold&Count>=sigCountThreshold)

cancerSampleCounts = sampleCancerTypes %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
rowIndex = data.frame(as.numeric(as.character(rownames(cancerSampleCounts))))
colnames(rowIndex) <- c("CancerIndex")
cancerSampleCounts = cbind(rowIndex, cancerSampleCounts)


# prepare summary for plot

# filter out residuals categories
cancerSampleSigSummary = (highMLSampleSigData %>% filter(SigName!='Excess'&SigName!='Unalloc') %>% group_by(CancerType,SigName) 
                          %>% summarise(CancerSigSampleCount=sum(SigPercent>=sigPercentThreshold),
                                        CancerSigSampleMutLoad=round(sum(Count),0)))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSampleCounts,by='CancerType',all.x=T)

# work counts of sample by cancer and signature
cancerTypesBySig = cancerSampleSigSummary %>% group_by(SigName) %>% summarise(SigCancerTypeCount=sum(CancerSigSampleCount>0))
rowIndex = data.frame(as.numeric(as.character(rownames(cancerTypesBySig))))
colnames(rowIndex) <- c("SigIndex")
cancerTypesBySig = cbind(rowIndex, cancerTypesBySig)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerTypesBySig,by='SigName',all.x=T)

# work out median mutational load per sig and cancer type for samples with the signature
cancerSigMedianCounts = highMLSampleSigData %>% filter(Count>0) %>% group_by(SigName,CancerType) %>% summarise(MedianSampleCount=round(median(Count),0))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSigMedianCounts,by=c('SigName','CancerType'),all.x=T)

# calculate counts per MB for each signature
mutLoadMin = 0.25
mutLoadMax = 16
cancerSampleSigSummary = cancerSampleSigSummary %>% mutate(SamplesWithSig=round(CancerSigSampleCount/CancerSampleCount,3),
                                                           MedianLoadPerMb=ifelse(MedianSampleCount>0,round(MedianSampleCount/3e3,3),0),
                                                           MedianLoadPerMb=ifelse(MedianLoadPerMb>mutLoadMax,mutLoadMax,
                                                                                  ifelse(MedianLoadPerMb<mutLoadMin,mutLoadMin,MedianLoadPerMb)))

plotBreaks = c(0.25, 0.5, 1, 2, 4, 8, 16)

plot = (ggplot(cancerSampleSigSummary %>% filter(SamplesWithSig>0), aes(x=CancerType, y=reorder(SigName, -SigIndex))) 
        + geom_point(aes(color=MedianLoadPerMb, size=SamplesWithSig)) 
        + scale_color_gradient(name='MedianLoadPerMb', trans='log2', low="lightblue", high="midnightblue", guide="colourbar",limits=c(mutLoadMin,mutLoadMax), breaks=plotBreaks, labels=plotBreaks)
        + theme_gray()
        + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
        + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0,size=10))
        + scale_x_discrete(position = 'top')
        + labs(subtitle="Signature Allocation by Cancer Types", y="Signature", x="CancerType"))

print(plot)


#####################
# 2. Residuals Report

residualsSummary = (snvDpSampleSigData %>% group_by(SampleId,CancerType) 
                    %>% summarise(SampleTotal=first(SampleTotal), 
                                  ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
                    %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4)))

resPlot = (ggplot(residualsSummary, aes(CancerType, ResidualPerc))
           + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
           + geom_boxplot(varwidth=T, fill="lightblue") +labs(title="Residuals % by Cancer Type", x="Cancer Type", y="Residuals %"))

print(resPlot)


######################################################
# 3. Bucket plots for samples with residuals above 50%

# these sample IDs need to be replaced with the HMF IDs
worstSamples = c('CPCT02010772T','CPCT02010662T') # these are the 2 SYD985 samples
worstSamplesMatrixData = snvDpSampleCounts %>% filter(SampleId %in% worstSamples) %>% spread(SampleId,Count)
worstSamplesMatrixData[is.na(worstSamplesMatrixData)] <- 0
worstSamplesMatrixData = within(worstSamplesMatrixData, rm(Bucket))

plot_96_profile(worstSamplesMatrixData)
