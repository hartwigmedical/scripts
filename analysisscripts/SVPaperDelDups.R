library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
library(grid)
library(gridExtra)
library(cowplot)

svData = read.csv('~/data/sv/CLUSTER.csv')
sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
svData$SampleId_CancerType = paste(svData$SampleId, svData$CancerType, sep='_')
View(head(svData,100))

svData = sv_set_common_fields(svData) # OR just run:
svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)

# extract only simple, unclustered DELs and DUPs
delsAndDups = svData %>% filter(Type=='DEL'|Type=='DUP') %>% filter(Length>0&ResolvedType=='SimpleSV')
delsAndDups$LenBucket = 2 ** round(log(delsAndDups$Length,2))

sampleDelDupCounts = delsAndDups %>% group_by(SampleId_CancerType,Type) %>% summarise(Count=n()) %>% arrange(-Count)
View(sampleDelDupCounts)

# PDF of DEL and DUP counts across all samples
delLengthPlot = (ggplot(data = sampleDelDupCounts %>% filter(Type=='DEL'&Count>=100), aes(x=reorder(SampleId_CancerType, -Count), y=Count))
      + geom_point()
      + labs(title = "DEL Counts All Samples") + xlab("Samples")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

print(delLengthPlot)

dupLengthPlot = (ggplot(data = sampleDelDupCounts %>% filter(Type=='DUP'&Count>=50), aes(x=reorder(SampleId_CancerType, -Count), y=Count))
      + geom_point()
      + labs(title = "DUP Counts All Samples") + xlab("Samples")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

print(dupLengthPlot)


getDeregulatedCutoff<-function(data,percentile,maxSamples=50)
{
  rowCount = nrow(data)
  if(rowCount == 0)
    return (0)
  
  topN = round(rowCount * percentile,0)
  topN = min(topN,maxSamples)
  topSamples = head(data %>% arrange(-Count), topN)
  cutoffRow = topSamples[nrow(topSamples),]
  return (cutoffRow$Count)
}

# delCountDeregulated = 250
# dupCountDeregulated = 250 # possibly higher
cutoffPercentile = 0.05
delCountDeregulated = getDeregulatedCutoff(sampleDelDupCounts %>% filter(Type=='DEL'), cutoffPercentile)
dupCountDeregulated = getDeregulatedCutoff(sampleDelDupCounts %>% filter(Type=='DUP'), cutoffPercentile)
print(paste("deregulated sample count cut-offs: DEL=", delCountDeregulated, ", DUP=", dupCountDeregulated, sep=''))

deregulatedSamples = (sampleDelDupCounts %>% filter((Type=='DEL'&Count>=delCountDeregulated)|(Type=='DUP'&Count>=dupCountDeregulated)) 
                      %>% group_by(SampleId_CancerType) %>% summarise(DelCount=sum(ifelse(Type=='DEL',Count,0)),
                                                                      DupCount=sum(ifelse(Type=='DUP',Count,0)),
                                                                      Total=sum(Count)))
View(deregulatedSamples)

deregDelDupData = (delsAndDups %>% filter(SampleId_CancerType %in% deregulatedSamples$SampleId_CancerType)
                   %>% group_by(SampleId_CancerType,Type,LenBucket) %>% summarise(Count=n()) %>% spread(Type,Count))

deregDelDupData = merge(deregDelDupData, deregulatedSamples, by="SampleId_CancerType", all.x=T)
deregDelDupData = deregDelDupData %>% arrange(-Total,LenBucket)

View(deregDelDupData)

deregDelDupData <- mutate(deregDelDupData,
                          SampleId_CancerType = reorder(SampleId_CancerType, -Total))

deregDelDupSamplesPlot = (ggplot(data = deregDelDupData, aes(x=LenBucket, y=Count))
      + geom_line(aes(y=DUP, colour='DUP'))
      + geom_line(aes(y=DEL, colour='DEL'))
      + scale_x_log10()
      + facet_wrap(~SampleId_CancerType)
      + labs(title = "Deregulated DEL and DUP Samples"))

print(deregDelDupSamplesPlot)

# then get all other samples and plot by cancer type
nonderegDelDupData = (delsAndDups %>% filter(!(SampleId_CancerType %in% deregulatedSamples$SampleId_CancerType))
                      %>% group_by(CancerType,Type,LenBucket) %>% summarise(Count=n()) %>% spread(Type,Count))

View(nonderegDelDupData)

nonderegDelDupSamplesByCancerTypePlot = (ggplot(data = nonderegDelDupData, aes(x=LenBucket, y=Count))
      + geom_line(aes(y=DUP, colour='DUP'))
      + geom_line(aes(y=DEL, colour='DEL'))
      + scale_x_log10()
      + facet_wrap(~CancerType)
      + labs(title = "Non-deregulated DEL and DUP lengths by cancer type"))

print(nonderegDelDupSamplesByCancerTypePlot)




# Synthetic DEL-DUP comparisons 

delExtTIs = svData %>% filter(grepl('DEL_Ext_TI',ResolvedType))
delExtTIs$LenBucket = ifelse(delExtTIs$SynDelDupLen<=32,32,2**round(log(delExtTIs$SynDelDupLen,2)))

View(delExtTIs)

dupExtTIs = svData %>% filter(grepl('DUP_Ext_TI',ResolvedType))
dupExtTIs$LenBucket = ifelse(dupExtTIs$SynDelDupLen<=32,32,2**round(log(dupExtTIs$SynDelDupLen,2)))

View(delExtTIs %>% select(SampleId,ClusterId,ResolvedType,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,
                        SynDelDupLen,SynDelDupTILen,LnkLenStart,LnkLenEnd,DBLenStart,DBLenEnd))

delExtTIs$LenBucket = ifelse(delExtTIs$SynDelDupLen<=32,32,2**round(log(delExtTIs$SynDelDupLen,2)))
delExtTIs$TILenBucket = ifelse(delExtTIs$SynDelDupTILen<=32,32,2**round(log(delExtTIs$SynDelDupTILen,2)))

View(delExtTIs %>% filter(SynDelDupTILen <= 1000) %>% group_by(LenBucket) %>% count())
View(delExtTIs %>% group_by(TILenBucket) %>% count())

# plot count of simple vs synthetic DELs by sample
delsAndDups$DelDupType = delsAndDups$Type
delExtTIs$DelDupType = "SYN_DEL"
dupExtTIs$DelDupType = "SYN_DUP"

combinedDelsAndDups = delsAndDups
combinedDelsAndDups = rbind(combinedDelsAndDups, delExtTIs %>% filter(ChainIndex==0))  # take only the first of the 2 SVs involved
combinedDelsAndDups = rbind(combinedDelsAndDups, dupExtTIs %>% filter(ChainIndex==0))  

delDupCountsPerSample = combinedDelsAndDups %>% group_by(SampleId,DelDupType) %>% summarise(Count=n()) %>% spread(DelDupType,Count)
delDupCountsPerSample[is.na(delDupCountsPerSample)] <- 0
View(delDupCountsPerSample)

delTypeComparisonPlot = (ggplot(data = delDupCountsPerSample, aes(x=DEL, y=SYN_DEL))
                 + geom_point()
                 + labs(title = "Simple vs Synthetic DELs Counts per Sample") 
                 + xlab("Simple DELs") + ylab("Synthetic DELs"))

print(delTypeComparisonPlot)

dupTypeComparisonPlot = (ggplot(data = delDupCountsPerSample, aes(x=DUP, y=SYN_DUP))
                         + geom_point()
                         + labs(title = "Simple vs Synthetic DUPs Counts per Sample") 
                         + xlab("Simple DUPs") + ylab("Synthetic DUPs"))

print(dupTypeComparisonPlot)




# PDF OUTPUT
outputFile = paste("~/logs/r_output/pdfs/SVA_DEL_DUP_Lengths.pdf", sep = "")

# DATA OUTPUT TO PDF
pdf(file=outputFile, height = 14, width = 20)

par(mar=c(1,1,1,1))

grid.arrange(delLengthPlot, dupLengthPlot, ncol = 1, nrow = 2, newpage = TRUE)

grid.arrange(deregDelDupSamplesPlot, ncol = 1, nrow = 1, newpage = TRUE)

grid.arrange(nonderegDelDupSamplesByCancerTypePlot, ncol = 1, nrow = 1, newpage = TRUE)

grid.arrange(delTypeComparisonPlot, dupTypeComparisonPlot, ncol = 2, nrow = 1, newpage = TRUE)

dev.off()





