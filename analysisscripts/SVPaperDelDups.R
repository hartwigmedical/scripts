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

svData = sv_load_and_prepare('~/data/sv/CLUSTER.csv')
svData$SampleId_CancerType = paste(svData$SampleId, svData$CancerType, sep='_')


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


# Synthetic DEL-DUPs

tiDirectData = read.csv("~/logs/SVA_LINKS.csv")
tiDirectData = tiDirectData %>% filter(TILength>=30)

synDelDups = tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'|ResolvedType=='DUP_Int_TI'|ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI')
synDelDups = synDelDups %>% filter(SynDelDupLen > 0) # keep these separate for now

View(synDelDups)

synDelDups$LengthBucket = 2**round(log(synDelDups$SynDelDupLen,2))

synDelDupsLengthSummary = synDelDups %>% group_by(LengthBucket,ResolvedType) %>% count() %>% spread(ResolvedType,n)
View(synDelDupsLengthSummary)

print(ggplot(data = synDelDups %>% group_by(LengthBucket,ResolvedType) %>% count(), aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ResolvedType))


synDelIntTIs = synDelDups %>% filter(ResolvedType=='DEL_Int_TI'&TraversedSVCount==0)
View(synDelIntTIs)
synDelIntTIs$IsReciprocalInv = (synDelIntTIs$ClusterDesc=='INV=2' & synDelIntTIs$TILength > 0.95 * synDelIntTIs$SynDelDupLen)
synDelIntTIs$LenBucket = ifelse(synDelIntTIs$SynDelDupLen<=32,32,2**round(log(synDelIntTIs$SynDelDupLen,2)))

print(ggplot(data = synDelIntTIs %>% group_by(LengthBucket,IsReciprocalInv) %>% count(), aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~IsReciprocalInv))

View(synDelIntTIs %>% filter(IsReciprocalInv))



# Comparison with simple DELs

delExtTIs = svData %>% filter(ResolvedType=='DEL_Ext_TI')
delExtTIs$LenBucket = ifelse(delExtTIs$SynDelDupLen<=32,32,2**round(log(delExtTIs$SynDelDupLen,2)))

View(delExtTIs)

View(delExtTIs %>% select(SampleId,ClusterId,ResolvedType,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,
                        SynDelDupLen,SynDelDupTILen,LnkLenStart,LnkLenEnd,DBLenStart,DBLenEnd))

delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainCount>0)
delIntTIs$LenBucket = ifelse(delIntTIs$SynDelDupLen<=32,32,2**round(log(delIntTIs$SynDelDupLen,2)))
nrow(delIntTIs)

synDelIntTIs$IsReciprocalInv = (synDelIntTIs$ClusterDesc=='INV=2' & synDelIntTIs$TILength > 0.99 * synDelIntTIs$SynDelDupLen)


# plot count of simple vs synthetic DELs by sample
delsAndDups$DelDupType = delsAndDups$Type
delsAndDups$LenBucket = ifelse(delsAndDups$SynDelDupLen<=32,32,2**round(log(delIntTIs$SynDelDupLen,2)))
delExtTIs$DelDupType = delExtTIs$ResolvedType
delIntTIs$DelDupType = ifelse(delIntTIs$ClusterDesc=='INV=2'&delIntTIs$SynDelDupTILen >= 0.99*delIntTIs$SynDelDupLen,'RecipInv',as.character(delIntTIs$ResolvedType))
View(delIntTIs %>% group_by(DelDupType)%>% count())

combinedDels = delsAndDups %>% filter(DelDupType=='DEL')
combinedDels = rbind(combinedDels, delExtTIs %>% filter(ChainIndex=='0s'))  # take only the first of the 2 SVs involved
combinedDels = rbind(combinedDels, delIntTIs %>% filter(ChainIndex=='0s'))  

View(combinedDels %>% group_by(DelDupType) %>% count())

plot_length_facetted(combinedDels, "DelDupType!='DEL'", 
                     "LenBucket,DelDupType", 
                     'LenBucket', 'DelDupType', "Simple and Synthetic DEL Lengths")

combinedDels$LengthGroup = ifelse(combinedDels$LenBucket<1e2,'Short',ifelse(combinedDels$LenBucket>2e3,'Medium','Long'))

delsPerSample = combinedDels %>% group_by(SampleId,DelDupType,LengthGroup) %>% summarise(Count=n()) %>% spread(DelDupType,Count)
delsPerSample[is.na(delsPerSample)] <- 0
View(delsPerSample)

print(ggplot(data = delsPerSample, aes(x=DEL, y=DEL_Ext_TI))
     + geom_point() + facet_wrap(~LengthGroup))

print(ggplot(data = delsPerSample, aes(x=DEL, y=DEL_Int_TI))
      + geom_point() + facet_wrap(~LengthGroup))

print(ggplot(data = delsPerSample, aes(x=DEL, y=RecipInv))
      + geom_point() + facet_wrap(~LengthGroup))




dupExtTIs$DelDupType = "SYN_DUP"


dupExtTIs = svData %>% filter(grepl('DUP_Ext_TI',ResolvedType))
dupExtTIs$LenBucket = ifelse(dupExtTIs$SynDelDupLen<=32,32,2**round(log(dupExtTIs$SynDelDupLen,2)))

dupTypeComparisonPlot = (ggplot(data = delDupCountsPerSample, aes(x=DUP, y=SYN_DUP))
                         + geom_point()
                         + labs(title = "Simple vs Synthetic DUPs Counts per Sample") 
                         + xlab("Simple DUPs") + ylab("Synthetic DUPs"))

print(dupTypeComparisonPlot)




# PDF OUTPUT
outputFile = paste("~/logs/r_output/pdfs/SVA_DEL_DUP_Lengths.pdf", sep = "")

pdf(file=outputFile, height = 14, width = 20)

par(mar=c(1,1,1,1))

grid.arrange(delLengthPlot, dupLengthPlot, ncol = 1, nrow = 2, newpage = TRUE)

grid.arrange(deregDelDupSamplesPlot, ncol = 1, nrow = 1, newpage = TRUE)

grid.arrange(nonderegDelDupSamplesByCancerTypePlot, ncol = 1, nrow = 1, newpage = TRUE)

grid.arrange(delTypeComparisonPlot, dupTypeComparisonPlot, ncol = 2, nrow = 1, newpage = TRUE)

dev.off()





