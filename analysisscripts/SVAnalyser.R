library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)

#' Modified version of dplyr's filter that uses string arguments
#' @export
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}

#' Modified version of dplyr's select that uses string arguments
#' @export
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @export
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}
#' Internal function used by s_filter, s_select etc.
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df
}


cohortSummary<-function(cluster,filterString = "",groupByString = "")
{
  (cluster %>% s_filter(filterString) %>% s_group_by(groupByString)
             %>% summarise(count=n(),
                           countSGL=sum(Type=='SGL'),
                           countNONE=sum(Type=='NONE'),
                           countBND=sum(Type=='BND'),
                           countINV=sum(Type=='INV'),
                           countDEL=sum(Type=='DEL'),
                           countDUP=sum(Type=='DUP'))
             %>% arrange(-count) %>% as.data.frame)
}

plot_count_by_bucket_and_type<-function(countsData,bucket,facetWrap,titleString ="",useLogX = TRUE) {
  plot <- ggplot(data=countsData,aes_string(x=bucket))+geom_line(aes(y=countDEL,colour='DEL'))+
    geom_line(aes(y=countDUP,colour='DUP'))+geom_line(aes(y=countINV,colour='INV'))+geom_line(aes(y=countBND,colour='BND'))+geom_line(aes(y=countSGL, colour='SGL'))+
    facet_wrap(as.formula(paste("~", facetWrap)))+
    labs(title = titleString)+ theme(panel.grid.major = element_line(colour="grey", size=0.5))
  if (useLogX == TRUE) {
    plot<-plot+scale_x_log10()
  }
  print(plot)
}

logScatterPlot<-function(data,xVar,YVar) {
  print(ggplot(data=data,aes_string(xVar,YVar))+geom_point() +scale_x_log10()+scale_y_log10())
}

sv_set_common_fields<-function(cluster)
{ cluster %>% mutate( 
    IsLINE = ifelse(LEStart!='None'|LEEnd!='None',T,F),
    IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F),
    Length = ifelse(as.character(ChrStart)!=as.character(ChrEnd)|Type=='INS'|ArmEnd!=ArmStart, -1, PosEnd-PosStart),
    DoubleDupBE = ifelse(DupBEStart=='true'&DupBEEnd=='true',T,F),
    SingleDupBE = ifelse(DoubleDupBE==0&(DupBEStart=='true'|DupBEEnd=='true'),T,F),
    TICount = ifelse(LnkTypeStart=='TI',0.5,0)+ifelse(LnkTypeEnd=='TI',0.5,0),
    DBCount = ifelse(DBLenStart>=0,0.5,0)+ifelse(DBLenEnd>=0,0.5,0),
    IsSglTI = ifelse(LnkTypeStart=='SGL',0.5,0),
    AsmbTICount = ifelse(AsmbMatchStart=='MATCH',0.5,0)+ifelse(AsmbMatchEnd=='MATCH',0.5,0),
    InferTICount = TICount - AsmbTICount,
    ShortTICount=ifelse(LnkTypeStart=='TI'&LnkLenStart<=1000,0.5,0)+ifelse(LnkTypeEnd=='TI'&LnkLenEnd<=1000,0.5,0),
    ClusterSize = ifelse(ClusterCount==1,'Single',ifelse(ClusterCount<=4,'Small','Large')),
    IsConsistent = ifelse(Consistency==0,T,F),
    IsChained = (ChainCount>=1),
    FoldbackCount = ifelse(FoldbackLenStart>=0,0.5,0)+ifelse(FoldbackLenEnd>=0,0.5,0),
    RepeatedChainLink = (svData$ChainCount>0 & grepl(';',svData$ChainIndex)),
    isPolyA = grepl('TTTTTTTT',InsertSeq)|grepl('AAAAAAAA',InsertSeq),
    isLowQual = ResolvedType=='LowQual'
  )

}

createBuckets <- function(cluster) {
  cluster %>% mutate(
    PloidyBucket=2**(pmin(5,pmax(-3,round(log(Ploidy,2),0)))),
    CnChEndBucket=2**(pmin(5,pmax(-3,round(log(pmax(AdjCNChgEnd,0.01),2),0)))),
    CnChStartBucket=2**(pmin(5,pmax(-3,round(log(pmax(AdjCNChgStart,0.01),2),0)))),
    ClusterCountBucket=2**(pmin(5,pmax(-3,round(log(ClusterCount,2),0)))),
    LnkLenStartBucket=ifelse(!LnkLenStart>0,0,2**(pmin(20,pmax(0,round(log(LnkLenStart,2),0))))),
    FoldBackLenStartBucket=ifelse(!FoldbackLenStart>0,0,2**(pmin(20,pmax(0,round(log(FoldbackLenStart,2),0))))),
    LnkLenEndBucket=ifelse(!LnkLenEnd>0,0,2**(pmin(20,pmax(0,round(log(LnkLenEnd,2),0))))),
    FoldBackLenEndBucket=ifelse(!FoldbackLenEnd>0,0,2**(pmin(20,pmax(0,round(log(FoldbackLenEnd,2),0))))),
    LengthBucket=ifelse(Type=='BND'|Type=='INS'|PosEnd-PosStart==0|ArmEnd!=ArmStart,0,2**(pmin(25,pmax(0,round(log(PosEnd-PosStart,2),0))))),
    HomLenBucket=2**(round(log(nchar(Homology),2),0)),
    NearestLenBucket=2**(round(log((NearestLen),2),0)),
    InsLenBucket=2**(round(log(nchar(InsertSeq),2),0))
  )
}

createFoldbacks <- function(cluster) {
  rbind(cluster %>% filter(FoldbackLenStart>=0) %>% mutate(FoldbackLength=FoldbackLenStart,FoldbackLinkInfo=FoldbackLinkInfoStart),
        cluster %>% filter(FoldbackLenEnd>=0)  %>% mutate(FoldbackLength=FoldbackLenEnd,FoldbackLinkInfo=FoldbackLinkInfoEnd)) %>% 
  separate(FoldbackLinkInfo,c('ChainLinks','AssemblyLinks','FoldbackChainLength'),sep = ';') %>%
  mutate(FoldbackType = ifelse(!is.na(FoldbackLnkStart)&!is.na(FoldbackLnkEnd)&FoldbackLnkStart==FoldbackLnkEnd&Type=='INV','INV','Combo'),
         FoldbackLenBucket = 2**round(log(FoldbackLength,2)),
         ChainLinks = as.numeric(ChainLinks),
         ChainLinksBucket = ifelse(ChainLinks>0,2**round(log(ChainLinks,2)),0),
         AssemblyLinks = as.numeric(AssemblyLinks),
         FoldbackChainLength = as.numeric(FoldbackChainLength),
         FoldbackChainLengthBucket = ifelse(FoldbackChainLength>0,2**round(log(FoldbackChainLength,2)),0),
         FoldbackAsmbPercent = ifelse(ChainLinks>0,round(AssemblyLinks/ChainLinks/0.2)*0.2,0),
         ChainSize = ifelse(ChainLinks==0,'None',ifelse(ChainLinks==1,'Single',ifelse(ChainLinks<=3,'Small','Long')))
  )
}

createSampleSummary <- function(cluster) {
 cluster %>% 
   group_by(SampleId,ClusterId,ResolvedType,Type=ifelse(ClusterCount==1,as.character(Type),''),ShortTICount) %>% 
   summarise(n=n(),FSCount=sum(IsFS)) %>% group_by(SampleId) %>% 
   summarise(countCluster=n(),
              count=sum(n),
              countLine=sum(ResolvedType=='Line'),
              sumLine=sum(ifelse(ResolvedType=='Line',n,0)),
              countSGLPair_INS=sum(ResolvedType=='SglPair_INS'),
              countSGLPair_DEL=sum(ResolvedType=='SglPair_DEL'),
              countSGLPair_DUP=sum(ResolvedType=='SglPair_DUP'),
              countSimpleChain=sum(ResolvedType=='SimpleChain'|ResolvedType=='SimplePartialChain'),
              sumSimpleChain=sum(ifelse(ResolvedType=='SimpleChain'|ResolvedType=='SimplePartialChain',n,0)),
              countComplexChain=sum(ResolvedType=='ComplexChain'|ResolvedType=='ComplexPartialChain'),
              sumComplexChain=sum(ifelse(ResolvedType=='ComplexChain'|ResolvedType=='ComplexPartialChain',n,0)),
              countSimpleINS=sum(Type=='INS'),
              countSimpleDEL=sum(Type=='DEL'),
              countSimpleDUP=sum(Type=='DUP'),
              countSimpleBND=sum(Type=='BND'),
              countSimpleINV=sum(Type=='INV'),
              countSimpleSGL=sum(Type=='SGL'),
              countSimpleNONE=sum(Type=='NONE'),
              countDEL_Ext_TI=sum(ResolvedType=='DEL_Ext_TI'),
              countDEL_Int_TI=sum(ResolvedType=='DEL_Int_TI'),
              countDUP_Ext_TI=sum(ResolvedType=='DUP_Ext_TI'),
              countDUP_Int_TI=sum(ResolvedType=='DUP_Int_TI'),
              countRecipTrans=sum(ResolvedType=='RecipTrans'),
              countShortTI=sum(ShortTICount),
              countFS=sum(FSCount))
}
#############################################
################ LOADING ####################
#############################################

svData = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER.csv', header = T, stringsAsFactors = F)
svClusters = (read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/SVA_CLUSTERS.csv', header = T, stringsAsFactors = F))
sampleCancerTypesFileName='~/Dropbox/HMF Australia team folder/Structural Variant Analysis/sample_cancer_types.csv'
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)

#Annotations
svData=createBuckets(svData)
svData=sv_set_common_fields(svData)
svSampleSummary=createSampleSummary(svData)
foldbcks=createFoldbacks(svData)
dbData = rbind(svData %>% filter(DBLenStart>-31) %>% mutate(DBLength = DBLenStart,Assembled = ifelse(AsmbMatchStart=="MATCH","Assembled","NotAssembled")),
               svData %>% filter(DBLenEnd>-31) %>% mutate(DBLength = DBLenEnd, Assembled = ifelse(AsmbMatchEnd=="MATCH","Assembled","NotAssembled"))) %>%
               mutate(DBLenBucket = ifelse(DBLength==0,0,ifelse(DBLength<0,-(2**round(log(-DBLength,2))),2**round(log(DBLength,2)))))


########################################################
########## Deletion Bridge Analysis ####################
########################################################

#1.A DBLength by ResolvedType => less than 50 bases
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50",'DBLength,ResolvedType'),'DBLength','ResolvedType','DBLengthByResolvedType(<50bases)',useLogX = F)

#1.B DBLength by ResolvedType => large scale
plot_count_by_bucket_and_type(cohortSummary(dbData,"",'DBLenBucket,ResolvedType'),'DBLenBucket','ResolvedType','DBLengthByResolvedType',useLogX = T)

#2. PolyA/T analysis across cluster type  
View( dbData %>% group_by(IsLINE,ResolvedType,isPolyA) %>% count() %>% spread(isPolyA,n,fill=0))
#TO DO: why so many POLY A and T in unexpected cluster types

#3. The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific
dbDataLINE = dbData %>% filter(DBLength<=50,ResolvedType %in% c('Line','SglPair_INS'))
View(dbDataLINE %>% group_by(CancerType,OLPeak=DBLength<(-7)) %>% count() %>% spread(OLPeak,n,fill=0))
View(dbDataLINE %>% group_by(SampleId,OLPeak=DBLength<(-7)) %>% count() %>% spread(OLPeak,n,fill=0))

#TO DO: could the 2 peaks stratify bg identiifcation of source line element in the ref genome?

################################################
####### FOLDBACK ANALYSIS #######################
#################################################

#1.Simple + foldback length distribution for short combos
print(ggplot(data = foldbacks %>% filter(FoldbackType=="INV"|FoldbackChainLengthBucket<10000) %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=INV, colour='INV'))
      + geom_line(aes(y=Combo, colour='Combo > 50% assembled'))
      + scale_x_log10()
      + labs(title = "Foldback Length Distribution"))

#2. Anlaysis of combo foldbacks by Chainlinks
print(ggplot(data = foldbacks %>% filter(FoldbackType=="Combo") %>% group_by(FoldbackLenBucket,LongChain=FoldbackChainLengthBucket>10000,ChainLinksBucket) %>% summarise(Count=n()) %>% spread(LongChain,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=`TRUE`, colour='Chain Length > 10k'))
      + geom_line(aes(y=`FALSE`, colour='Chain Length < 10k'))
      + scale_x_log10()
      + facet_wrap(~ChainLinksBucket)
      + labs(title = "Combo Foldback Length Distribution by ChainLinksBucket"))      

#3. Anlaysis of combo foldbacks by ChainLength
print(ggplot(data = foldbacks %>% filter(FoldbackType=="Combo") %>% group_by(FoldbackLenBucket,LongChain=FoldbackChainLengthBucket>10000) %>% summarise(Count=n()) %>% spread(LongChain,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=`TRUE`, colour='Chain Length > 10k'))
      + geom_line(aes(y=`FALSE`, colour='Chain Length < 10k'))
      + scale_x_log10()
      + labs(title = "Combo Foldback Length Distribution by ChainLinksBucket"))    

#4. Count of Foldbacks per cluster
View(foldbacks %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(FBcount=n()/2,
          ClusterCount=first(ClusterCount),CNMax=max(pmax(AdjCNStart,AdjCNEnd)),FBPloidyMax=max(pmax(AdjCNChgStart,AdjCNChgEnd,Ploidy))))

###################################################
########### SAMPLE LEVEL ANALYSES #################
###################################################

#1. ScatterPlots of various counts
logScatterPlot(svSampleSummary,'countLine','sumLine')
logScatterPlot(svSampleSummary,'countLine','countSGLPair_INS')
logScatterPlot(svSampleSummary,'countSimpleDUP','countDUP_Ext_TI')
logScatterPlot(svSampleSummary,'countSimpleDEL','countDEL_Int_TI')
logScatterPlot(svSampleSummary,'countSimpleDEL','countDEL_Int_TI')
logScatterPlot(svSampleSummary,'countDEL_Int_TI','countRecipTrans')
logScatterPlot(svSampleSummary,'countFS','countSimpleDEL')

#2. Basic signature type plot top 50 samples
ResolvedTypeColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown", "darkgreen", "deepskyblue", "tan","black","orange","lightgreen")
filter = svSampleSummary %>% arrange (-count) %>% filter(row_number() <= 50) %>% .$SampleId

print(ggplot(data = svData %>% filter(SampleId %in% filter) %>% group_by(SampleId,ResolvedType=ifelse(grepl('Chain',ResolvedType),'Chain',ResolvedType)) %>% 
               summarise(SvCount=n()) %>% group_by(SampleId) %>% mutate(SampleCount=sum(SvCount)) %>% ungroup() %>% arrange(-SampleCount), 
       aes(x = reorder(SampleId, -SampleCount), y = SvCount, fill = ResolvedType))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = ResolvedTypeColours)
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))
      + ylab("SV Count") + xlab("Sample"))

#############################################
########## OVERVIEW OF COHORT ###############
#############################################

# 1. By RESOLOVED Type
View(svData %>% filter(!isLowQual) %>% group_by(ResolvedType,Type) %>% count() %>% spread(Type,n))

# 2. By CLUSTERCOUNTBUCKET
View(svData %>% filter(!isLowQual,ResolvedType!='Line') %>% group_by(ClusterCountBucket,Type) %>% count() %>% spread(Type,n))

#############################################
########## INDIVIDUAL SAMPLE ANALYSIS #######
#############################################

chartSample = 'CPCT02020670T'; chartChr = '17';

#1. Whole Chr View
View(svData %>% filter(sampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr,ResolvedType!='ALowQual') %>% 
       select(ResolvedType,LengthBucket,Ploidy,Type,ClusterCount,ClusterReason,everything()))
print(ggplot(data = svData %>% filter(!isLowQual,SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>% mutate(
  modStart=ifelse(ChrStart==chartChr,PosStart,as.numeric(ChrStart)*-1e6),modEnd=ifelse(ChrEnd==chartChr,PosEnd,as.numeric(ChrEnd)*-1e6)),aes(modStart,modEnd)) + 
    geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())

#1a. all chromosomes
print(ggplot(data = svData %>% filter(!isLowQual,SampleId==chartSample) %>% mutate(
  modStart=PosStart,modEnd=ifelse(ChrEnd==ChrStart,PosEnd,as.numeric(ChrEnd)*-1e6)),aes(modStart,modEnd)) + 
    geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw() +facet_wrap(~ChrStart))

#2.By CHR for large clusters
View(svData %>% filter(SampleId == chartSample,ClusterCount>8) %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% count() %>% spread(chr_arm_end,n,fill=""))

# 3. Single Cluster View
View(svData %>% filter(SampleId==chartSample,ClusterId==7)%>% arrange(ChrStart,PosStart) 
     %>% select(ResolvedType,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,Type,ClusterCount,ClusterReason,everything()))

# 4. Large cluster View
View(svData %>% filter(!isLowQual,SampleId==chartSample,ClusterCount>=1,ClusterCount>0,ResolvedType!='SimpleSV'))#%>% group_by(ClusterCountBucket,CnChStartBucket) %>% count() %>% spread(CnChStartBucket,n))
View(temp)

##############################################################
########## COMPLEX CLUSTERS AND CLUSTERS PER ARM #############
##############################################################

# 1. Per arm by ClusterCount
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>% filter(!isLowQual,Type!='NONE',ResolvedType!='Line') %>%
             group_by(SampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count(),svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>% filter(!isLowQual,Type!='NONE',ResolvedType!='Line') %>%
             group_by(SampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count()) %>% 
       group_by(SampleId,chr_arm,CC,id) %>% count() %>%
       group_by(SampleId,CC,chr_arm) %>% count() %>% spread(CC,n) %>%
       filter(chr_arm!='0_P'))

#2 LINE Clusters
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count(),svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count()) %>% 
       filter(ResolvedType=='Line',ClusterCount>8) %>%
       group_by(SampleId,id,chr_arm) %>% summarise(n=sum(n)/2) %>%
       group_by(SampleId,id) %>% summarise(arms=n(),svCount=sum(n)) )

#3 NON-LINE Clusters
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count(),svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count()) %>% 
       filter(ResolvedType!='Line',ClusterCount>8) %>%
       group_by(SampleId,id,chr_arm) %>% summarise(n=sum(n)/2) %>%
       group_by(SampleId,id) %>% summarise(arms=n(),svCount=sum(n)) )


#4. COUNT of NON_LINE PER SAMPLE
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count(),svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>%
             group_by(SampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count()) %>% 
       filter(ResolvedType!='Line',ClusterCount>8) %>%
       group_by(SampleId,id,chr_arm) %>% summarise(n=sum(n)/2) %>%
       group_by(SampleId,id) %>% summarise(arms=n(),svCount=sum(n)) %>% group_by(SampleId) %>% count())

#5. COUNT of NON_resolved
View(svData %>% filter(ResolvedType!='SimpleSV',!isLowQual,ResolvedType!='Line',!grepl('TI',ResolvedType),!grepl('Sgl',ResolvedType),ResolvedType!='RecipTrans',Type!='NONE',ClusterCount==1) %>% group_by(SampleId,Type) %>% count() %>% spread(Type,n))

#6. COUNT of clusters sharing inter-arms
View(svData %>% filter(ChrStart==ChrEnd,ArmStart!=ArmEnd,!isLowQual,ResolvedType!='Line') %>% group_by(SampleId,ChrStart,ClusterId) %>% count() %>% group_by(SampleId,ChrStart) %>% count())

##############################################################
########## ANALYSIS OF SIMPLE SVs ############################
##############################################################

#1. top50 simple
filter = cohortSummary(svData,"ResolvedType=='SimpleSV'",'SampleId')  %>% filter(row_number() <= 50) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"ResolvedType=='SimpleSV',!isLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','AllDel&Dup')

#2. By Cancer Type excluding top 50
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(!(SampleId %in% filter)),"ResolvedType=='SimpleSV'",'CancerType,LengthBucket'),'LengthBucket','CancerType','AllDel&Dup')

##############################################################
########## ANALYSIS OF LONE INCONSISTENT VARTIANTS ###########
##############################################################

#1. LONE INV
View(svData %>% filter(ClusterCount==1,Type =='INV',!isLowQual,AdjCNChgEnd>0.5,AdjCNChgStart>0.5)%>%  mutate(len=PosEnd-PosStart) %>% arrange(SampleId,ChrStart,PosStart) %>%
     select(ResolvedType,len,Ploidy,Type,ClusterCount,ClusterReason,everything()))

#2. LONE BND
View(svData %>% filter(ClusterCount==1,Type=='BND',!isLowQual,AdjCNChgStart>0.5,AdjCNChgEnd>0.5) %>% group_by(SampleId)%>%  count())

#3. LONE INS
View(svData %>% filter(ClusterCount==1,Type=='INS') %>% mutate(len=PosEnd-PosStart,insertLen=nchar(InsertSeq)) %>% select(len,insertLen,everything()))

#4. Excess SGL
View(svData %>% filter(AdjCNStart>0,!isLowQual) %>% group_by(SampleId) %>% summarise(count=n(),countSGL=sum(ifelse(Type=='SGL',1,0)),
      countLT0.2VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.2,1,0)),countLT0.15VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.15,1,0)),
      countLT0.1VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.1,1,0)),maxCN=max(AdjCNStart)) %>% arrange(-countLT0.15VafSGL))

View(svData %>% filter(AdjCNStart>0,!isLowQual) %>% group_by(ResolvedType,PolyA=grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)) %>% summarise(count=n(),countSGL=sum(ifelse(Type=='SGL',1,0)),
      countLT0.2VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.2,1,0)),countLT0.15VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.15,1,0)),
      countLT0.1VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.1,1,0)),maxCN=max(AdjCNStart)) %>% arrange(-countLT0.15VafSGL))

################################################################
###################### OTHER ###################################
################################################################

#1.  INS SEQUENCE DISTRIBTUTION

View(svData %>% filter(!isLowQual,ResolvedType!='Line') %>% group_by(InsLenBucket,Type) %>% count() %>% spread(InsLenBucket,n))

#2. What to do about POLYG?
View(svData %>% filter(grepl('CCCCCCCCCCCC',InsertSeq)|grepl('GGGGGGGGGGGG',InsertSeq),!isLowQual))
View(svData %>% filter(grepl('GGGGGGGGGGGG',InsertSeq)|grepl('CCCCCCCCCCCC',InsertSeq),!isLowQual) %>% group_by(round(AdjCNStart,-1),Type) %>% count() %>% spread(Type,n))

