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

svData =  read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F)#gridssCohortVariantsread.csv('~/data/sv/CLUSTER.csv')
sampleCancerTypes = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/sample_cancer_types.csv')
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
svData$SampleId_CancerType = paste(svData$SampleId, svData$CancerType, sep='_')

svData = sv_set_common_fields(svData) 


# Calculate DB Lengths per breakend
dbStartLens = svData %>% filter(DBLenStart>-31)
dbStartLens$DBLength = dbStartLens$DBLenStart
dbStartLens$Assembled = ifelse(dbStartLens$AsmbMatchStart=="MATCH","Assembled","NotAssembled")
dbEndLens = svData %>% filter(DBLenEnd>-31)
dbEndLens$DBLength = dbEndLens$DBLenEnd
dbEndLens$Assembled = ifelse(dbEndLens$AsmbMatchStart=="MATCH","Assembled","NotAssembled")
dbData = rbind(dbStartLens, dbEndLens)
dbData$DBLenBucket = ifelse(dbData$DBLength==0,0,ifelse(dbData$DBLength<0,-(2**round(log(-dbData$DBLength,2))),2**round(log(dbData$DBLength,2))))

#1. DBLength by ResolvedType
print(ggplot(data = dbData %>% filter(DBLength<=50) %>% group_by(DBLength,Type,ResolvedType) %>% summarise(Count=n()) %>% spread(Type,Count),
      aes(x=DBLength, y=Count))
      + geom_line(aes(y=BND, colour='BND'))
      + geom_line(aes(y=DEL, colour='DEL'))
      + geom_line(aes(y=DUP, colour='DUP'))
      + geom_line(aes(y=INV, colour='INV'))
      + geom_line(aes(y=SGL, colour='SGL'))
      + theme(panel.grid.major = element_line(colour="grey", size=0.5))
      + facet_wrap(~ResolvedType))

#2. PolyA/T analysis across cluster type
View( dbData %>% group_by(IsLINE,ResolvedType,RepeatedSeq=ifelse(grepl('TTTTTTTT',InsertSeq)|grepl('AAAAAAAA',InsertSeq),'PolyA/T','None')) %>% 
        count() %>% spread(RepeatedSeq,n,fill=0))
#TO DO: why so many POLY A and T in unexpected cluster types

#3. The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific
View(dbData %>% filter(DBLength<=50,ResolvedType=='Line'|ResolvedType=='SglPair_INS') %>%
       group_by(CancerType,OLPeak=DBLength<(-7)) %>% count() %>% spread(OLPeak,n,fill=0))
View(dbData %>% filter(DBLength<=50,ResolvedType=='Line'|ResolvedType=='SglPair_INS') %>%
       group_by(SampleId,OLPeak=DBLength<(-7)) %>% count() %>% spread(OLPeak,n,fill=0))


#TO DO: could the 2 peaks stratify bg identiifcation of source line element in the ref genome?

########################################################
########## Templated Insertion Analysis ################
########################################################



svFile = '~/data/sv/CLUSTER.csv'
svData = sv_load_and_prepare(svFile)


# TI Direct File
rm(tiDirectData)

tiDirectData = read.csv("~/logs/SVA_LINKS.csv")
nrow(tiDirectData)
nrow(tiDirectData %>% filter(TILength>=30&TILength<=500))

View(tiDirectData)

nrow(tiDirectData %>% filter(TILength<30))
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
tiDirectData$CNGain = (tiDirectData$CopyNumberGain=='true')
tiDirectData$ArmOfOrigin = (tiDirectData$OnArmOfOrigin=='true')
tiDirectData$TraversesSvs = (tiDirectData$TraversedSVCount>0)
#View(tiDirectData %>% group_by(TraversedSVCount) %>% count())
#View(tiDirectData %>% group_by(TraversesSvs) %>% count())

tiDirectData$ChainCountSize = ifelse(tiDirectData$ChainCount==2,'1',ifelse(tiDirectData$ChainCount<=4,'2-3','4+'))

tiDirectData$ClusterType = tiDirectData$ResolvedType
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DUP_Ext_TI'|tiDirectData$ClusterType=='DEL_Ext_TI','DEL-DUP_EXT',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DEL_Int_TI','DEL_INT',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='DUP_Int_TI','SIMPLE',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='ComplexChain','COMPLEX',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='None','COMPLEX',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='SimpleChain','SIMPLE',as.character(tiDirectData$ClusterType))
tiDirectData$ClusterType = ifelse(tiDirectData$ClusterType=='SimpleSV','SIMPLE',as.character(tiDirectData$ClusterType))

# View(tiDirectData %>% group_by(ClusterType) %>% count())

tiDirectData$NextSV = ifelse(tiDirectData$NextSVDistance==-1,'NONE',ifelse(tiDirectData$NextSVDistance<=clusterDistance,'5K','REM'))
tiDirectData$NextSVTraversesSvs = (tiDirectData$NextSVTraversedCount>0)

minDBLength = -30
clusterDistance = 5000

tiDirectData$DBLengths = ifelse(tiDirectData$DBLenStart < minDBLength & tiDirectData$DBLenEnd < minDBLength,'NONE',
                            ifelse((tiDirectData$DBLenStart >= minDBLength & tiDirectData$DBLenStart <= 50)|(tiDirectData$DBLenEnd >= minDBLength & tiDirectData$DBLenEnd <= 50),'DB',
                            ifelse((tiDirectData$DBLenStart >= minDBLength & tiDirectData$DBLenStart <= clusterDistance)|(tiDirectData$DBLenEnd >= minDBLength & tiDirectData$DBLenEnd <= clusterDistance),'5K','NONE')))

# 1. Majority of templated insertions are short (250 - 1000 bases)

# split synthetic DELs and DUPs from longer clusters
plot_length_facetted(tiDirectData, "ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'", "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "TI Lengths for Synthetic DELs and DUPS")

plot_length_facetted(tiDirectData, "ResolvedType=='ComplexChain'|ResolvedType=='SimpleChain'", "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "TI Lengths for 3+ CLusters")


# the second peak (or more) beyond the short TIs is largely from shattering and repaired sections, not genuine inserted fragments
# show this by splitting by whether the TI has another SV in the same cluster on the same arm

plot_length_facetted(tiDirectData, "NextSVDistance<=0&(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI')", "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "Isolated TI Lengths for Synthetic DELs and DUPS")

plot_length_facetted(tiDirectData, "NextSVDistance<=0&(ResolvedType=='ComplexChain'|ResolvedType=='SimpleChain')", "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "Isolated TI Lengths for 3+ Clusters")

plot_length_facetted(tiDirectData, "NextSVDistance>0&(ResolvedType=='ComplexChain'|ResolvedType=='SimpleChain')&ArmOfOrigin&!CNGain&DBLengths!='NONE'", 
                     "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "Clustered TI Lengths for 3+ Clusters")

plot_length_facetted(tiDirectData, "NextSVDistance>0&(ResolvedType=='ComplexChain'|ResolvedType=='SimpleChain')", 
                     "TiLenBucket,ResolvedType", 
                     'TiLenBucket', 'ResolvedType', "Clustered TI Lengths for 3+ Clusters")

# 2. Distribution of chain link lengths based on proportion of short & assembled TI
shortTIData = (tiDirectData %>% group_by(SampleId,ClusterId,ChainId) 
               %>% summarise(ClusterCount=first(ClusterCount),
                             LinkCount=n(),
                             ShortTICount=sum(TILength<=1e3),
                             AssembledCount=sum(IsAssembled=='true')))

shortTIData$ShortTIPerc = shortTIData$ShortTICount/(shortTIData$ClusterCount-1)
shortTIData$ShortTIPercBucket = ifelse(shortTIData$ShortTIPerc>1,1,round(shortTIData$ShortTIPerc/0.2)*0.2)
shortTIData$AssembledTIPerc = shortTIData$AssembledCount/(shortTIData$ClusterCount-1)
shortTIData$AssembledTIPercBucket = ifelse(shortTIData$AssembledTIPerc>1,1,round(shortTIData$AssembledTIPerc/0.2)*0.2)

#View(shortTIData)
totalChains = nrow(shortTIData)
View(shortTIData %>% group_by(ShortTIPercBucket) %>% summarise(Count=n(),AsPerc=round(n()/totalChains,2)))
View(shortTIData %>% group_by(AssembledTIPercBucket) %>% summarise(Count=n(),AsPerc=round(n()/totalChains,2)))

print(ggplot(data = shortTIData %>% filter(ShortTIPerc>=0.9&LinkCount>1) %>% group_by(LinkCount) %>% count(), aes(x=LinkCount, y=n))
      + geom_line()
      + labs(title = "Chain Link Count for short consequetive TIs"))



############################
# working and experimental...

# short inferred TIs

shortTIData = tiDirectData %>% filter(TILength<250&ClusterCount<50&!grepl('Sgl',ResolvedType))
View(shortTIData)
View(shortTIData %>% group_by(ResolvedType) %>% count())

shortTISummary = shortTIData %>% group_by(TiLenBucket,AssemblyType) %>% count()

print(ggplot(data = shortTISummary, aes(x=TiLenBucket, y=n))
        + geom_line()
        + facet_wrap(~AssemblyType)
        + labs(title = 'TI Lengths assembled or not'))

shortTISum2 = shortTIData %>% filter(ResolvedType!='ComplexChain') %>% group_by(TILength,AssemblyType,NextSV=NextSVDistance>0) %>% count() %>% spread(AssemblyType,n)
View(shortTISum2)

print(ggplot(data = shortTISum2, aes(x=TILength, y=n))
      + geom_line(aes(y=INFERRED, colour='INFERRED'))
      + geom_line(aes(y=ASSEMBLY, colour='ASSEMBLY'))
      + facet_wrap(~NextSV)
      + labs(title = 'TI Lengths assembled or not'))


shortTIData = tiDirectData %>% filter(TILength<1e3&ClusterCount<50) %>% group_by(TiLenBucket,AssemblyType,NextSV=NextSVDistance>0) %>% count()


View(tiDirectData %>% filter(AssemblyType=='INFERRED'&TILength<=100))



# TIs be whether they traverse another SV or not
travsSvData = (tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI') 
               %>% group_by(ResolvedType,TraversesSvs,TiLenBucket) %>% count() %>% spread(ResolvedType,n))
View(travsSvData)

print(ggplot(data = travsSvData, aes(x=TiLenBucket, y=n))
      + geom_line(aes(y=DEL_Ext_TI, colour='DEL_Ext_TI'))
      + geom_line(aes(y=DUP_Ext_TI, colour='DUP_Ext_TI'))
      + scale_x_log10()
      + facet_wrap(~TraversesSvs)
      + labs(title = "TI Lengths for Traversing SVs"))

travsSvData2 = (tiDirectData %>% filter(ResolvedType=='DUP_Int_TI') 
                %>% group_by(ClusterDesc,TraversesSvs,TiLenBucket) %>% count() %>% spread(TraversesSvs,n))
travsSvData2 = (tiDirectData %>% filter(ClusterDesc=='DUP=2'&ResolvedType=='DUP_Ext_TI') 
               %>% group_by(ClusterDesc,TraversesSvs,TiLenBucket) %>% count() %>% spread(TraversesSvs,n))

print(ggplot(data = travsSvData2 , aes(x=TiLenBucket, y=n))
      + geom_line(aes(y=`TRUE`, colour='TraversesSVs'))
      + geom_line(aes(y=`FALSE`, colour='NoTraversal'))
      + scale_x_log10()
      + facet_wrap(~ClusterDesc)
      + labs(title = "TI Lengths for Traversing SVs"))


travsSvData3 = (tiDirectData %>% filter(ResolvedType=='SimpleChain'|ResolvedType=='ComplexChain') 
               %>% group_by(ResolvedType,TraversesSvs,TiLenBucket) %>% count() %>% spread(ResolvedType,n))

print(ggplot(data = travsSvData3, aes(x=TiLenBucket, y=n))
      + geom_line(aes(y=ComplexChain, colour='ComplexChain'))
      + geom_line(aes(y=SimpleChain, colour='SimpleChain'))
      + scale_x_log10()
      + facet_wrap(~TraversesSvs)
      + labs(title = "TI Lengths for Traversing SVs"))




tiDirectData$Connectivity = ifelse(tiDirectData$CNGain&tiDirectData$TraversesSvs==F&tiDirectData$ArmOfOrigin==F
                                   &(tiDirectData$DBLengths=='NONE'|tiDirectData$DBLengths=='REM')&(tiDirectData$NextSV=='REM'|tiDirectData$NextSV=='NONE'),'ISOLATED','LINKED')
View(tiDirectData %>% group_by(Connectivity) %>% count())

tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$NextSV,' DB=',tiDirectData$DBLengths, sep='')

print(ggplot(data = tiDirectData %>% group_by(Category,TiLenBucket) %>% summarise(Count=n()), aes(x=TiLenBucket, y=Count))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "TI Lengths by Combined Category"))


# 0. Connectivity

tiDirectData$Connectivity = ifelse(tiDirectData$CNGain&tiDirectData$TraversesSvs==F&tiDirectData$ArmOfOrigin==F
                                   &(tiDirectData$DBLengths=='NONE'|tiDirectData$DBLengths=='REM')&(tiDirectData$NextSV=='REM'|tiDirectData$NextSV=='NONE'),'ISOLATED','LINKED')

tiDirectData$TILengthGroup = ifelse(tiDirectData$TILength<1e3,'Short',ifelse(tiDirectData$TILength<=1e4,'LT5K','Long'))


View(tiDirectData %>% filter(ClusterDesc=='INV=2'&ResolvedType=='DEL_Ext_TI')
                                %>% group_by(TraversesSvs,DBLengths,NextSV,NextSVTraversesSvs,CNGain,TILengthGroup,ResolvedType) 
                                %>% count() %>% spread(TILengthGroup,n) %>% arrange(!TraversesSvs,!CNGain))

tiDirectData$Category = paste(tiDirectData$ClusterType,' CNGain=',tiDirectData$CNGain, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'&ClusterDesc=='INV=2'), F)
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'), F)
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='COMPLEX'), F)

View(tiDirectData %>% filter(ClusterDesc=='INV=2'&ResolvedType=='DEL_Ext_TI'&!TraversesSvs&!NextSVTraversesSvs))

tiDirectData$Category = paste(tiDirectData$ResolvedType,' TITravSV=',tiDirectData$TraversesSvs, ' NextSVTravSV=',tiDirectData$NextSVTraversesSvs, sep='')


tiLenSynDelDupByConnectivity = (tiDirectData %>% filter(ClusterDesc=='INV=2'&ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI')
                          %>% group_by(TraversesSvs,DBLengths,NextSV,NextSVTraversesSvs,CNGain,TILengthGroup,ResolvedType) 
                          %>% count() %>% spread(TILengthGroup,n) %>% arrange(!TraversesSvs,!CNGain))

View(tiLenSynDelDupByConnectivity)

tiLenChainsByConn = (tiDirectData 
       %>% filter((ClusterType=='SIMPLE'|ClusterType=='COMPLEX')&ResolvedType!='DUP_Int_TI')
       %>% group_by(TraversesSvs,DBLengths,NextSV,CNGain,TILengthGroup,ArmOfOrigin,ChainCountSize) 
       %>% count() %>% spread(TILengthGroup,n) %>% arrange(!TraversesSvs,!CNGain))

tiLenChainsByConn[is.na(tiLenChainsByConn)] = 0
tiLenChainsByConn$TotalCount = tiLenChainsByConn$Long+tiLenChainsByConn$Short+tiLenChainsByConn$LT5K
tiLenChainsByConn$ShortLongRatio = round(tiLenChainsByConn$Short/tiLenChainsByConn$TotalCount,2)

View(tiLenChainsByConn)
View(tiLenChainsByConn %>% filter(ChainCountSize=='4+'))


View(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'&NextSV=='NONE'))


tiLenSynLongChainsLinkByConnectivity = (tiDirectData 
                                          %>% filter(ClusterType=='SIMPLE'&ResolvedType!='DUP_Int_TI'&(ChainCountSize=='1'|ChainCountSize=='2-3'))
                                          %>% group_by(TraversesSvs,DBLengths,NextSV,CNGain,TILengthGroup,ArmOfOrigin) 
                                          %>% count() %>% spread(TILengthGroup,n) %>% arrange(!TraversesSvs,!CNGain) 
                                          %>% mutate(ShortLongRatio = ifelse(Long>0,round(Short/(Long+Short),2),1)))

View(tiLenSynSimpleChainsLinkByConnectivity)

# 1. Synthetic DELs and DUPs
View(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'))

tiDirectData$Category = paste(tiDirectData$ClusterType,' DB=',tiDirectData$DBLengths, sep='')
#plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'), T)
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI'), F)

tiLengthsByConnecivity = (tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI')
     %>% group_by(TraversesSvs,DBLengths,NextSV,CNGain,TILengthGroup,ResolvedType) 
     %>% count() %>% spread(TILengthGroup,n) %>% arrange(!TraversesSvs,!CNGain))

View(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI')
     %>% group_by(TraversesSvs,DBLengths,NextSV,CNGain,ShortTI=TILength<1e3) %>% count() %>% spread(ShortTI,n))

View(tiDirectData %>% filter((ResolvedType=='DEL_Ext_TI'|ResolvedType=='DUP_Ext_TI')
                             &CNGain&!TraversesSvs&DBLengths=='NONE'&(NextSV=='REM'|NextSV=='NONE')))


# 2. DELs Internal TI - limit to INV pairs
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'))
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'&DBLengths=='NONE'))
View(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI') %>% group_by(ClusterDesc) %>% count())

View(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'&ClusterDesc=='INV=2'))

plot_length_facetted(tiDirectData, "ResolvedType=='DEL_Ext_TI'&ClusterDesc=='INV=2'", "TiLenBucket", 
                     'TiLenBucket', '', "TI Lengths for DEL_Ext_TI from INV pairs")



tiDirectData$Category = paste(tiDirectData$ClusterType,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Int_TI'&ClusterDesc=='INV=2'), F)

print(ggplot(data = tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI',!NextSVTraversesSvs), aes(NextSVDistance))
           + stat_ecdf()
      + scale_x_log10()
      + facet_wrap(~ClusterDesc))

tiDirectData$NextSVDistBucket = ifelse(tiDirectData$NextSVDistance>0,10**round(log(tiDirectData$NextSVDistance,10)),-1)
tiDirectData$NextSVDistBucket = ifelse(tiDirectData$NextSVDistance<=0,'NONE',ifelse(tiDirectData$NextSVDistance<=5e3,'LT5K',
                                  ifelse(tiDirectData$NextSVDistance<1e5,'LT100K','REM')))

tiDirectData$Category = paste(' NextSVDistance=',tiDirectData$NextSVDistBucket, sep='')
tiDirectData$Category = paste(' NextSVDistance=',tiDirectData$ClusterType, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType!='aDUP_Ext_TI'), F)
tiDirectData$Category = paste(tiDirectData$ClusterType,' NextSVDistance=',tiDirectData$NextSVDistBucket, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterCount<4), F)

tiDirectData$Category = paste(tiDirectData$ClusterType,' DB Distance=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'), F)

tiDirectData$Category = paste(tiDirectData$ClusterType,' CopyNumberGain=',tiDirectData$CNGain, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'), F)

tiDirectData$Category = paste(tiDirectData$ClusterType,' TraversesSvs=',tiDirectData$TraversesSvs, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'), F)

tiDirectData$Category = paste(tiDirectData$ClusterType,' NextSVTraversesSvs=',tiDirectData$NextSVTraversesSvs, sep='')
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='DEL_Ext_TI'), F)

View(tiDirectData %>% group_by(NextSVDistBucket) %>% count())
View(tiDirectData %>% filter(NextSVDistBucket=='NONE'))



# 3. Simple Chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$NextSV,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'), F)

View(tiDirectData %>% filter(FullyChained=='true') %>% group_by(ResolvedType) %>% count())

# just looking at connectivity
tiDirectData$Category = paste('Connectivity=',tiDirectData$Connectivity, sep='')
plot_ti_by_category(tiDirectData %>% filter(FullyChained=='true'), F)

tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' CONN=',tiDirectData$Connectivity, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&FullyChained=='true'), F)
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize!='4+'), F)

# single-link chains
View(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize==1))
tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$NextSV,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize==1), F)

# 2-3 link chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$NextSV,' DB=',tiDirectData$DBLengths, ' CNG=', tiDirectData$CNChange, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='2-3'), F)

View(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='2-3'))

tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' DB=',tiDirectData$DBLengths, ' CNG=', tiDirectData$CNChange, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='2-3'), F)

View(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSZ=='2-3'))
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='2-3'), F)

# longer link chains
tiDirectData$Category = paste(tiDirectData$ClusterType,' SZ=',tiDirectData$ChainCountSize,' NextSV=',tiDirectData$NextSV,' DB=',tiDirectData$DBLengths, sep='')
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='SIMPLE'&ChainCountSize=='4+'), F)

View(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='4+'))
plot_ti_by_category(tiDirectData %>% filter(ResolvedType=='SimpleChain'&ChainCountSize=='4+'), F)

# 4. Complex clusters
plot_ti_by_category(tiDirectData %>% filter(ClusterType=='COMPLEX'&ChainCountSize=='4+'), F)

plot_ti_by_category(tiDirectData %>% filter(ClusterType=='COMPLEX'&FullyChained=='true'), F)


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


# Isolated TIs which are also long

View(tiDirectData %>% filter(FullyChained=='true') %>% group_by(Connectivity,Short=TILength<2e3) %>% count())
View(tiDirectData %>% filter(FullyChained=='true'&Connectivity=='ISOLATED'))
View(tiDirectData %>% filter(FullyChained=='true'&Connectivity=='ISOLATED') %>% group_by(TiLenBucket) %>% count())

View(tiDirectData %>% filter(FullyChained=='true'&Connectivity=='ISOLATED'))



view_cluster_sv_links('CPCT02010134T', 29)
view_chromosome_sv('CPCT02010134T', '5', '5')

View(svData %>% filter(SampleId=='CPCT02050182T'&ClusterId %in% c(2,22,18))
     %>% select(SampleId,ClusterId,Id,Type,ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd,LnkLenStart,LnkLenEnd,
                DBLenStart,DBLenEnd,AdjCNStart,AdjCNChgStart,Ploidy))  

view_cluster_sv_links('CPCT02050182T', 18)
view_cluster_sv_links('CPCT02050182T', 22)
view_cluster_sv_links('CPCT02050182T', 2)
view_chromosome_sv('CPCT02050182T', '1', '1')
view_chromosome_sv('CPCT02050182T', '3', '3')

view_clusters_sv('CPCT02050182T', c(2,18,22))



view_cluster_sv_links('CPCT02080125T', 244)
view_chromosome_sv('CPCT02080125T', '1', '1')
view_chromosome_sv('CPCT02080125T', '14', '14')



plot_cross_chr_sample_clusters(svData, c(1,3,5,17), 'CPCT02050182T', 1, F)


View(svData %>% filter(SampleId=='CPCT02070161T'&ClusterId==41) 
     %>% select(SampleId,ClusterId,Id,Type,ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd,LnkLenStart,LnkLenEnd,AsmbMatchStart,AsmbMatchEnd,
                DBLenStart,DBLenEnd,ChainIndex,Consistency,AdjCNStart,AdjCNChgStart,AdjCNEnd,AdjCNChgEnd,Ploidy))

View(tiDirectData %>% filter(Connectivity=='ISOLATED'&ResolvedType=='SimpleChain'&FullyChained=='true'))

view_chromosome_sv('CPCT02020250T', '1', '')

View(svData %>% filter(SampleId=='CPCT02020250T'&ChrStart==1) 
     %>% select(SampleId,ClusterId,ClusterDesc,ResolvedType,Id,Type,ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd,LnkLenStart,LnkLenEnd,AsmbMatchStart,AsmbMatchEnd,
                DBLenStart,DBLenEnd,ChainIndex,Consistency) %>% arrange(PosStart))

View(svData %>% filter(SampleId=='CPCT02020250T'&ChrStart==1) %>% group_by(ClusterDesc,ResolvedType) %>% count())




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






