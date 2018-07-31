library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)


cohortSummary<-function(cluster,filterString = "",groupByString = "")
{
  summary = (cluster %>% s_filter(filterString) %>% s_group_by(groupByString)
             %>% summarise(count=n(),
                           DelCount=sum(Type=='DEL'),
                           DelPerc=round(sum(Type=='DEL')/n(),2),
                           DupCount=sum(Type=='DUP'),
                           DupPerc=round(sum(Type=='DUP')/n(),2),
                           InvCount=sum(Type=='INV'),
                           InvPerc=round(sum(Type=='INV')/n(),2),
                           BndCount=sum(Type=='BND'),
                           BndPerc=round(sum(Type=='BND')/n(),2),
                           FSPerc=round(sum(FSStart!='false'|FSEnd!='false')/n(),2),
                           LEPerc=round(sum(LEStart!='false'|LEEnd!='false')/n(),2),
                           QualityPerc=round(sum(!(Ploidy > 0.5 & (AdjCNChgStart < (Ploidy - 0.5) | AdjCNChgEnd < (Ploidy - 0.5))))/n(),2),
                           LenLT10KPerc=round(sum(Length > -1 & Length <= 1e4)/n(),2),
                           Len10Kto100KPerc=round(sum(Length > 1e4 & Length <= 1e5)/n(),2),
                           Len100Kto500KPerc=round(sum(Length > 1e5 & Length < 5e5)/n(),2),
                           Len500Kto5MPerc=round(sum(Length > 5e5 & Length < 5e6)/n(),2),
                           ClusteringPerc=round(sum(ClusterCount > 1)/n(),2),
                           SharedPerc=round(sum(MultipleBiopsy=='Shared')/n(),2),
                           HomologyPerc=round(sum(stri_length(Homology) > 0)/n(),4),
                           AvgClusterCount=round(sum(ClusterCount)/n(),0),
                           AvgProximity=round(sum(ifelse(NearestLen >= 0, NearestLen, 1000000))/n(),-2),
                           AvgPloidy=round(sum(Ploidy)/n(),2))
             %>% arrange(PatientId) %>% as.data.frame)
}


#########################
######## CODE STARTS HERE
#########################

# SV data file
cluster = read.csv('~/logs/CLUSTER_V12.csv')

cluster$Length = ifelse(cluster$Type=='BND'|cluster$Type=='INS'|cluster$ArmEnd!=cluster$ArmStart, -1, cluster$PosEnd-cluster$PosStart)
cluster$LengthBucket=ifelse(cluster$Length <= 0, 0, 2**(round(log(cluster$Length,2),0)))

cluster$ArmStressedStart =  ifelse(cluster$ArmCountStart / cluster$ArmExpStart >= 5,'true','false')
cluster$ArmStressedEnd =  ifelse(cluster$ArmCountEnd / cluster$ArmExpEnd >= 5,'true','false')
cluster$IsStressed = ifelse(cluster$ArmStressedStart=='true'|cluster$ArmStressedStart=='true','true','false')

multipleBiopsySampleIds = read.csv('~/data/multiple_biopsy_patientids.csv')
# cluster$SampleIdStripped = substr(cluster$SampleId,0,stri_locate_last(cluster$SampleId, regex="T")-1)
cluster = (merge(cluster,multipleBiopsySampleIds,by.x="SampleId",by.y="SampleId",all.x=TRUE))

multipleSVs = cluster %>% filter(MultipleBiopsy=='Shared'|MultipleBiopsy=='Private')

View(multipleSVs)

# standard sample summaries just for MB records, order by patientId
mbSampleSummary = cohortSummary(multipleSVs, "", "SampleId,PatientId")
View(mbSampleSummary)
