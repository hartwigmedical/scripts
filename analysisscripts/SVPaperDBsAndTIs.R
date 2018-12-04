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

svData = sv_set_common_fields(svData) 


# Deletion Bridge length profile
dbStartLens = svData %>% filter(DBLenStart>-31)
dbStartLens$DBLength = dbStartLens$DBLenStart
dbStartLens$Assembled = ifelse(dbStartLens$AsmbMatchStart=="MATCH","Assembled","NotAssembled")
dbEndLens = svData %>% filter(DBLenEnd>-31)
dbEndLens$DBLength = dbEndLens$DBLenEnd
dbEndLens$Assembled = ifelse(dbEndLens$AsmbMatchStart=="MATCH","Assembled","NotAssembled")

dbData = rbind(dbStartLens, dbEndLens)
dbData$DBLenBucket = ifelse(dbData$DBLength==0,0,ifelse(dbData$DBLength<0,-(2**round(log(-dbData$DBLength,2))),2**round(log(dbData$DBLength,2))))

dbSummaryData = dbData %>% group_by(DBLenBucket,Assembled) %>% summarise(Count=n()) %>% spread(Assembled,Count)
View(dbSummaryData)

dbData1 = dbData %>% filter(DBLength<=200) %>% group_by(DBLength,Assembled) %>% summarise(Count=n()) %>% spread(Assembled,Count)
dbData2 = dbData %>% filter(DBLength<=100&Type!='SGL'&Type!='NONE'&Type!='BND') %>% group_by(DBLength,Assembled) %>% summarise(Count=n()) %>% spread(Assembled,Count)
dbData3 = dbData %>% filter(DBLength<=100&ResolvedType!='SglPair_INS') %>% group_by(DBLength,Assembled) %>% summarise(Count=n()) %>% spread(Assembled,Count)
dbData4 = dbData %>% filter(DBLength<=100&ResolvedType=='SglPair_INS') %>% group_by(DBLength) %>% summarise(Count=n())

dbTypeData = dbData %>% filter(DBLength<=50,ResolvedType=='ComplexPartialChain') %>% group_by(DBLength,Type,ResolvedType) %>% summarise(Count=n()) %>% spread(Type,Count)

print(ggplot(data = dbTypeData, aes(x=DBLength, y=Count))
      + geom_line(aes(y=BND, colour='BND'))
      + geom_line(aes(y=DEL, colour='DEL'))
      + geom_line(aes(y=DUP, colour='DUP'))
      + geom_line(aes(y=INV, colour='INV'))
      + geom_line(aes(y=SGL, colour='SGL'))
      + theme(panel.grid.major = element_line(colour="grey", size=0.5))
      + facet_wrap(~ResolvedType))

View(dbData %>% filter(DBLength<=50,ResolvedType=='DEL_Int_TI')%>% group_by(perfect=DBLength==1) %>% 
        count())
View(dbData %>% filter(DBLength<=50,ResolvedType=='DEL_Int_TI')%>% group_by(DBLength) %>% 
       count())
View(svData %>% filter(DBLength<=50,ResolvedType=='Line'))

+ theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  scale_y_continuous(minor_breaks = seq(0 , 100, 5), breaks = seq(0, 100, 10))


print(ggplot(data = dbData1, aes(x=DBLength, y=Count))
      + geom_line(aes(y=Assembled, colour='Assembled'))
      + geom_line(aes(y=NotAssembled, colour='Not-Assembled')))



print(ggplot(data = dbData4, aes(x=DBLength, y=Count))
      + geom_line())

sglPairIns = dbData %>% filter(DBLength<=100&ResolvedType=='SglPair_INS')
sglPairIns$InsertSeqLen = stri_length(sglPairIns$InsertSeq) 
sglPairIns$InsertSeqLenPlusDB = sglPairIns$InsertSeqLen + sglPairIns$DBLength

sglPairIns$RepeatedSeq = ifelse(grepl('TTTTTTTT',sglPairIns$InsertSeq),'PolyT',ifelse(grepl('AAAAAAA',sglPairIns$InsertSeq),'PolyA','None'))

View(sglPairIns %>% group_by(RepeatedSeq,IsLINE) %>% count())

View(sglPairIns)

print(ggplot(data = sglPairIns %>% filter(DBLength<50) %>% group_by(InsertSeqLen) %>% count(), aes(x=InsertSeqLen, y=n))
      + geom_line() + labs(title = "INS Sequence Length Distribution"))

print(ggplot(data = sglPairIns %>% filter(DBLength<50) %>% group_by(InsertSeqLenPlusDB) %>% count(), aes(x=InsertSeqLenPlusDB, y=n))
      + geom_line() + labs(title = "INS Sequence Length + DB Length Distribution"))


View(dbData %>% filter(DBLength<50) %>% group_by(DBLength) %>% count())
View(dbData %>% filter(DBLength<50) %>% group_by(DBLength,Assembled) %>% count() %>% spread(Assembled,n))

print(ggplot(data = dbData %>% filter(DBLength<50) %>% group_by(DBLength) %>% count(), aes(x=DBLength, y=n))
                          + geom_line()
                          + labs(title = "Deletion Bridge Length Distribution"))


print(ggplot(data = dbSummaryData, aes(x=DBLenBucket, y=n))
                + geom_line()
                + scale_x_log10()
                # + scale_x_continuous(trans=reverselog_trans(base=10), labels=trans_format("identity", function(x) -x))
                + labs(title = "Deletion Bridge Length Distribution"))

# exploring the peak around -16 bases for non-assembled DBs
nonAsmbDbData = dbData %>% filter(Assembled=="NotAssembled"&DBLength<50)

View(nonAsmbDbData %>% group_by(Type,DBLength) %>% count() %>% spread(Type,n))
nonAsmbTypeSummary = nonAsmbDbData %>% group_by(Type,DBLength) %>% summarise(Count=n())

print(ggplot(data = nonAsmbTypeSummary, aes(x=DBLength, y=Count))
                      + geom_line()
                      + facet_wrap(~Type)
                      + labs(title = "Non-assembled DB Lengths by Type"))

nonAsmbDbData$ClusterCountBucket = 2**round(log(nonAsmbDbData$ClusterCount,2))

nonAsmbClusterSizeSummary = nonAsmbDbData %>% group_by(ClusterCountBucket,DBLength) %>% summarise(Count=n())

print(ggplot(data = nonAsmbClusterSizeSummary, aes(x=DBLength, y=Count))
      + geom_line()
      + facet_wrap(~ClusterCountBucket)
      + labs(title = "Non-assembled DB Lengths by Cluster Size"))

nonAsmbResolvedTypeSummary = nonAsmbDbData %>% group_by(ResolvedType,DBLength) %>% summarise(Count=n())

print(ggplot(data = nonAsmbResolvedTypeSummary, aes(x=DBLength, y=Count))
      + geom_line()
      + facet_wrap(~ResolvedType)
      + labs(title = "Non-assembled DB Lengths by Resolved Type"))



# Templated Insertion Analysis

svFile = '~/data/sv/CLUSTER.csv'
svData = sv_load_and_prepare(svFile)


# TI Direct File
tiDirectData = read.csv("~/logs/CLUSTER_LINKS.csv")

tiDirectData = tiDirectData %>% filter(TILength>=30)

View(tiDirectData)
View(tiDirectData %>% group_by(ResolvedType) %>% count())
nrow(tiDirectData %>% filter(IsLINE=='true'))

# criteria for evaluating TI lengths
# TI length but becomes less reliable for larger clusters or those with copy number uncertainty
# - remoteness of TI - 
# - DB length or none - short (<50 bases), proximite (<5K), long or none
# - assmebled chain length - single, short (2-4), long (5+)

tiDirectData$TiLenBucket = 2**round(log(tiDirectData$TILength,2))

tiDirectData$AssemblyType = ifelse(tiDirectData$IsAssembled=='true','ASSEMBLY','INFERRED')

tiDirectData$ChainCountSize = ifelse(tiDirectData$ChainCount==2,'1',ifelse(tiDirectData$ChainCount<=4,'2-3','4+'))
tiDirectData$AssembledCountSize = ifelse(tiDirectData$AssembledCount==0,0,ifelse(tiDirectData$AssembledCount==1,'1',ifelse(tiDirectData$ChainCount<=3,'2-3','4+')))

tiDirectData$ClusterType = tiDirectData$ResolvedType
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DUP_Ext_TI'|tiDirectData$ClusterType=='DEL_Ext_TI','DEL-DUP_EXT',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DEL_Int_TI','DEL_INT',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DUP_Int_TI','SIMPLE',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='ComplexChain'|tiDirectData$ClusterType=='ComplexPartialChain','COMPLEX',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='None','COMPLEX',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='SimpleChain'|tiDirectData$ClusterType=='SimplePartialChain','SIMPLE',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='SimpleSV','SIMPLE',as.character(tiDirectData$ClusterType))

View(tiDirectData %>% group_by(ClusterType) %>% count())

tiDirectData$Remoteness = ifelse(tiDirectData$NextSVLength==-1,'NONE',ifelse(tiDirectData$NextSVLength<=clusterDistance,'5K','REMOTE'))

minDBLength = -30
clusterDistance = 5000

tiDirectData$DBLengths = ifelse(tiDirectData$DBLenStart < minDBLength & tiDirectData$DBLenEnd < minDBLength,'NONE',
                            ifelse((tiDirectData$DBLenStart >= minDBLength & tiDirectData$DBLenStart <= 50)|(tiDirectData$DBLenEnd >= minDBLength & tiDirectData$DBLenEnd <= 50),'DB',
                            ifelse((tiDirectData$DBLenStart >= minDBLength & tiDirectData$DBLenStart <= clusterDistance)|(tiDirectData$DBLenEnd >= minDBLength & tiDirectData$DBLenEnd <= clusterDistance),'5K','REMOTE')))

View(tiDirectData %>% group_by(DBLengths) %>% count())

tiDirectData$Category = paste(tiDirectData$ClusterType,' Size=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$Remoteness,' DB=',tiDirectData$DBLengths, sep='')

print(ggplot(data = tiDirectData %>% group_by(Category,TiLenBucket) %>% summarise(Count=n()), aes(x=TiLenBucket, y=Count))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "TI Lengths by Combined Category"))

# 1. Synthetic DELs and DUPs
View(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'))

tiDirectData$Category = paste(tiDirectData$ClusterType,' DB=',tiDirectData$DBLengths, sep='')
#plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'), T)
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'), F)

View(svData %>% filter(SampleId=='CPCT02010003T'&ClusterId==209) 
     %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,ClusterId,LnkTypeStart,LnkTypeEnd,AsmbMatchStart,AsmbMatchEnd,DBLenStart,DBLenEnd))

# 2. DELs Internal TI - limit to INV pairs
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'))
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'&DBLengths=='NONE'))
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI') %>% group_by(ClusterDesc) %>% count())

tiDirectData$Category = paste(tiDirectData$ClusterType,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'&ClusterDesc=='INV=2'), F)


# 3. Simple Chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' Size=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$Remoteness,' DB=',tiDirectData$DBLengths, sep='')
#plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'), T)
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'), F)

# single-link chains
View(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize==1))
View(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize==1&Remoteness=='5K'&DBLengths=='DB'))
tiDirectData$Category = paste(tiDirectData$ClusterType,' Size=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$Remoteness,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize==1), F)

# 2-3 link chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' Size=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$Remoteness,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='2-3'), F)

View(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='2-3'))
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='2-3'), F)

# longer link chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' Size=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$Remoteness,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='4+'), F)

View(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='4+'))
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='4+'), F)

# 4. Complex clusters
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='COMPLEX'&ChainCountSize=='4+'), F)


plot_ti_by_category<-function(tiData, showAssembly)
{
  if(showAssembly)
  {
    plotData = tiData %>% group_by(Category,AssemblyType,TiLenBucket) %>% summarise(Count=n()) %>% spread(AssemblyType,Count)
      
    print(ggplot(data = plotData, aes(x=TiLenBucket, y=Count))
          + geom_line(aes(y=ASSEMBLY, colour='ASSEMBLY'))
          + geom_line(aes(y=INFERRED, colour='INFERRED'))
          + scale_x_log10()
          + facet_wrap(~Category)
          + labs(title = "TI Lengths by Combined Category"))
  }
  else
  {
    plotData = tiData %>% group_by(Category,TiLenBucket) %>% summarise(Count=n())

    print(ggplot(data = plotData, aes(x=TiLenBucket, y=Count))
          + geom_line()
          + scale_x_log10()
          + facet_wrap(~Category)
          + labs(title = "TI Lengths by Combined Category"))
  }
}


print(ggplot(data = tiDirectData %>% group_by(Category,AssemblyType,TiLenBucket) %>% summarise(Count=n()) %>% spread(AssemblyType,Count), aes(x=TiLenBucket, y=Count))
      + geom_line(aes(y=ASSEMBLY, colour='ASSEMBLY'))
      + geom_line(aes(y=INFERRED, colour='INFERRED'))
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "TI Lengths by Combined Category"))



# assembled only
View(tiDirectData %>% filter(AssemblyType=='ASSEMBLY') %>% group_by(TILength) %>% summarise(Count=n()))
print(ggplot(data = tiDirectData %>% filter(AssemblyType=='ASSEMBLY') %>% group_by(TILenBucket=round(TILength/2)*2) %>% summarise(Count=n()), 
                          aes(x=TILenBucket, y=Count))
                   + geom_line()
                   + theme(panel.grid.major = element_line(colour="grey", size=0.5)))


tiLengthSummary = tiDirectData %>% filter(TiLengthBucket<500) %>% group_by(TiLenSimBucket,AssemblyType) %>% summarise(Count=n()) %>% spread(AssemblyType,Count)
tiLengthSummary = tiDirectData %>% group_by(TiLengthBucket,IsAssembled) %>% summarise(Count=n()) %>% spread(IsAssembled,Count)

tiAssembledPlot = (ggplot(data = tiDirectData %>% filter(TiLengthBucket < 1e4) %>% group_by(TiLengthBucket,AssemblyType) %>% summarise(Count=n()) %>% spread(AssemblyType,Count), 
                          aes(x=TiLengthBucket, y=Count))
                      + geom_line(aes(y=ASSEMBLY, colour='ASSEMBLY'))
                      + geom_line(aes(y=INFERRED, colour='INFERRED'))
                      + scale_x_log10()
                      + theme(panel.grid.major = element_line(colour="grey", size=0.5))
                      + labs(title = "TI Lengths by Cluster Type"))

print(tiAssembledPlot)






# pre-direct TI file

View(svData5K %>% filter(IsLINE==T,ClusterCount==1,Type!='BND')) 
View(svData5K %>% filter(SampleId=='CPCT02020258T',ChrStart==13))
svData = sv_load_and_prepare(svFile)
summary10k=(svData %>% filter(grepl('Complex',ResolvedType)) %>% group_by(SampleId,Type) %>% count() %>% spread(Type,n))
summary5k=(svData5K %>% filter(grepl('Complex',ResolvedType)) %>% group_by(SampleId,Type) %>% count() %>% spread(Type,n))
View(merge(summary10k,summary5k,by='SampleId',all=T) %>% mutate(BNDDiff=BND.y-BND.x,DUPDiff=DUP.y-DUP.x,DELDiff=DEL.y-DEL.x,INVDiff=INV.y-INV.x) %>% 
       select(SampleId,BND.x,BND.y,INV.x,INV.y,DEL.x,DEL.y,DUP.x,DUP.y,BNDDiff,DUPDiff,INVDiff,DELDiff))
totalSVCount = nrow(svData)
View(svData %>% group_by(ResolvedType,ClusterSize) 
     %>% summarise(Clusters=n_distinct(paste(SampleId,ClusterId,sep='_')), TotalSVs=n(), AsPerc=round(n()/totalSVCount,2))
     %>% arrange(-AsPerc))


tiDataStart = svData %>% filter(LnkTypeStart=='TI'&Type!='SGL'&Type!='NONE') %>% filter(IsLINE==F)
tiDataStart$TiId1 = ifelse(tiDataStart$Id<tiDataStart$LnkSvStart,tiDataStart$Id,tiDataStart$LnkSvStart)
tiDataStart$TiId2 = ifelse(tiDataStart$Id>tiDataStart$LnkSvStart,tiDataStart$Id,tiDataStart$LnkSvStart)
tiDataStart$TiLength = tiDataStart$LnkLenStart
tiDataStart$Assembly = tiDataStart$AsmbMatchStart
tiDataStart$DBOnOther = (tiDataStart$DBLenEnd>-31 & tiDataStart$DBLenEnd<100)
tiDataEnd = svData %>% filter(LnkTypeEnd=='TI'&Type!='SGL'&Type!='NONE') %>% filter(IsLINE==F)
tiDataEnd$TiId1 = ifelse(tiDataEnd$Id<tiDataEnd$LnkSvEnd,tiDataEnd$Id,tiDataEnd$LnkSvEnd)
tiDataEnd$TiId2 = ifelse(tiDataEnd$Id>tiDataEnd$LnkSvEnd,tiDataEnd$Id,tiDataEnd$LnkSvEnd)
tiDataEnd$TiLength = tiDataEnd$LnkLenEnd
tiDataEnd$Assembly = tiDataEnd$AsmbMatchEnd
tiDataEnd$DBOnOther = (tiDataEnd$DBLenStart>-31 & tiDataEnd$DBLenStart<100)

tiData = rbind(tiDataStart,tiDataEnd)

tiDataPairs = (tiData %>% group_by(SampleId,ClusterId,TiId1,TiId2) 
               %>% summarise(Count=n(),
                             TiLength=first(TiLength),
                             Assembly=first(Assembly),
                             BndCount=sum(Type=='BND'),
                             CrossArmCount=sum(ArmStart!=ArmEnd),
                             ChrStart1=first(ChrStart),
                             ChrStart2=last(ChrStart),
                             ArmStart1=first(ArmStart),
                             ArmStart2=last(ArmStart),
                             DBOnOtherCount=sum(DBOnOther),
                             Assembly=last(Assembly),
                             ResolvedType=first(ResolvedType),
                             SynDelDupTILen=first(SynDelDupTILen),
                             ClusterSize=first(ClusterSize))
               %>% filter(Count==2))

# View(tiDataPairs %>% filter(Count==3))
View(tiDataPairs)
View(tiDataPairs %>% group_by(Assembly,DBOnOtherCount) %>% count())
View(tiDataPairs %>%filter(BndCount==0&CrossArmCount==2))

View(tiDataPairs %>%filter(BndCount==0&CrossArmCount==2) %>% group_by(SampleId,ChrStart1,ChrStart2,ArmStart1,ArmStart2) 
     %>% summarise(n_distinct(ClusterId)))


tiDataPairs$TiLenBucket = 2**round(log(tiDataPairs$TiLength,2))
tiDataPairs$ClusterType = tiDataPairs$ResolvedType
tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='ComplexPartialChain','ComplexChain',as.character(tiDataPairs$ClusterType))
tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='None','ComplexChain',as.character(tiDataPairs$ClusterType))
tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='SimplePartialChain','SimpleChain',as.character(tiDataPairs$ClusterType))
tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='SimpleSV','SimpleChain',as.character(tiDataPairs$ClusterType))
tiDataPairs$AssemblyType = ifelse(tiDataPairs$Assembly=='MATCH','ASSEMBLY','INFERRED')
tiDataPairs$DbFlanked = ifelse(tiDataPairs$DBOnOtherCount==2,'DBFlanked','Isolated')

nrow(tiDataPairs %>% filter(AssemblyType=='ASSEMBLY'))
nrow(svData)

View(svData %>% filter(AsmbTICount>0&IsLINE==F) %>% group_by(FoldbackCount>0,Type) %>% count())
View(svData %>% filter(AsmbTICount>0&IsLINE==F) %>% group_by(ResolvedType,Type) %>% count() %>% spread(Type,n))
nrow(svData %>% filter(AsmbTICount>0&Type=='BND'))

tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='SimpleChain','Chain',as.character(tiDataPairs$ClusterType))
tiDataPairs$ClusterType = ifelse(tiDataPairs$ClusterType=='ComplexChain','Chain',as.character(tiDataPairs$ClusterType))

tiDataPairs$TiLocation = ifelse(tiDataPairs$ResolvedType=='DEL_Ext_TI'|tiDataPairs$ResolvedType=='DUP_Ext_TI'|tiDataPairs$BndCount==2,'Remote','Unclear')

tiDataPairs$TiCategory = paste(tiDataPairs$ClusterType, tiDataPairs$TiLocation, tiDataPairs$DbFlanked, sep='_')
# tiDataPairs$TiCategory = paste(tiDataPairs$ClusterType, tiDataPairs$TiLocation, tiDataPairs$DbFlanked, sep='_')
# tiDataPairs$TiCategory = paste(tiDataPairs$ClusterType, tiDataPairs$AssemblyType,tiDataPairs$TiLocation,sep='_')

tiDataPairs = tiDataPairs %>% filter(ResolvedType!='Line')

#View(tiDataPairs %>% group_by(TiLenBucket,ResolvedType) %>% count() %>% spread(ResolvedType,n))
#View(tiDataPairs %>% group_by(TiLenBucket,ClusterType) %>% count() %>% spread(ClusterType,n))

View(tiDataPairs %>% filter(ClusterType=='Chain'&TiLength>2000&TiLength<=10000&TiLocation=='Remote'))

# tiClusterSummary = tiDataPairs %>% group_by(TiLenBucket,ClusterType) %>% summarise(Count=n())
tiClusterSummary = tiDataPairs %>% group_by(TiLenBucket,TiCategory) %>% summarise(Count=n())
tiDataPairs$TiLenSimBucket = round(tiDataPairs$TiLength,-1)
tiClusterSummary = tiDataPairs %>% filter(TiLenSimBucket<500) %>% group_by(TiLenSimBucket,AssemblyType) %>% summarise(Count=n()) %>% spread(AssemblyType,Count)


tiClusterTypesPlot = (ggplot(data = tiClusterSummary, aes(x=TiLenSimBucket, y=Count))
                      + geom_line(aes(y=ASSEMBLY, colour='ASSEMBLY'))
                      + geom_line(aes(y=INFERRED, colour='INFERRED'))
                      + theme(panel.grid.major = element_line(colour="grey", size=0.5))
                      + labs(title = "TI Lengths by Cluster Type"))

print(tiClusterTypesPlot)

tiClusterSummary2 = tiDataPairs %>% filter(TiLenSimBucket<2000) %>% group_by(TiLenSimBucket,ResolvedType) %>% summarise(Count=n())
View(tiClusterSummary2)
print(ggplot(data = tiClusterSummary2, aes(x=TiLenSimBucket, y=Count))
      + geom_line()
      + facet_wrap(~ResolvedType)
      + theme(panel.grid.major = element_line(colour="grey", size=0.5))
      + labs(title = "TI Lengths by Cluster Type"))



print(ggplot(data = tiClusterSummary, aes(x=TiLenBucket, y=Count))
      + geom_line()
      + geom_line()
      + scale_x_log10()
      + labs(title = "TI Lengths by Cluster Type")))

print(tiClusterTypesPlot)





