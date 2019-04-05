library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)

myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")

base_complements <- function(bases) {
  complements = setNames(c("A", "C", "G","T"), c("T", "G", "C","A"))
  
  point_complement <- function(base) {
    paste(rev(sapply(strsplit(base, split = ""), function (x) {complements[x]})), collapse = "")
  }
  
  sapply(bases, point_complement)
}

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

insertAndHomologySummary<-function(cluster,filterString = "",groupByString = "") {
  (cluster %>% s_filter(filterString) %>% s_group_by(groupByString)
   %>% summarise(count=n(),
                 noInsnoHom=sum(HomLenBucket==0&InsLenBucket==0)/n(),
                 HomTotal=sum(HomLenBucket>0)/n(),
                 InsTotal=sum(InsLenBucket>0)/n(),
                 InsGTE32=sum(HomLenBucket==0&InsLenBucket>=32)/n(),
                 Ins16=sum(HomLenBucket==0&InsLenBucket==16)/n(),
                 Ins8=sum(HomLenBucket==0&InsLenBucket==8)/n(),
                 Ins4=sum(HomLenBucket==0&InsLenBucket==4)/n(),
                 Ins2=sum(HomLenBucket==0&InsLenBucket==2)/n(),
                 Ins1=sum(HomLenBucket==0&InsLenBucket==1)/n(),
                 Hom1=sum(HomLenBucket==1&InsLenBucket==0)/n(),
                 Hom2=sum(HomLenBucket==2&InsLenBucket==0)/n(),
                 Hom4=sum(HomLenBucket==4&InsLenBucket==0)/n(),
                 Hom8=sum(HomLenBucket==8&InsLenBucket==0)/n(),
                 Hom16=sum(HomLenBucket==16&InsLenBucket==0)/n(),
                 HomGTE32=sum(HomLenBucket>=32&InsLenBucket==0)/n())
  )
}


cohortSummary<-function(cluster,filterString = "",groupByString = "") {
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

plot_count_by_bucket_and_type<-function(countsData,bucket,facetWrap,titleString ="",useLogX = TRUE,useLogY = TRUE) {
  plot <- ggplot(data=countsData,aes_string(x=bucket))+geom_line(aes(y=countDEL,colour='DEL'))+
    geom_line(aes(y=countDUP,colour='DUP'))+geom_line(aes(y=countINV,colour='INV'))+geom_line(aes(y=countBND,colour='BND'))+geom_line(aes(y=countSGL, colour='SGL'))+
    facet_wrap(as.formula(paste("~", facetWrap)))+
    labs(title = titleString) +theme_bw() + theme(panel.grid.major = element_line(colour="grey", size=0.5))
  if (useLogX == TRUE) {
    plot<-plot+scale_x_log10()
  }
  if (useLogY == TRUE) {
    plot<-plot+scale_y_log10()
  }
  print(plot)
}

scatterPlot<-function(data,xVar,YVar,useLogX = TRUE,useLogY = TRUE,facetWrap='') {
  plot<-ggplot(data=data,aes_string(xVar,YVar))+geom_point() 
  if (useLogX == TRUE) {
    plot<-plot+scale_x_log10()
  }
  if (useLogY == TRUE) {
    plot<-plot+scale_y_log10()
  }
  if (facetWrap != '') {
    plot<-plot+facet_wrap(as.formula(paste("~", facetWrap,",scales='free'")))
  }
  print(plot)
}

sv_set_common_fields<-function(cluster){ cluster %>% mutate( 
    IsLINE = ifelse(LEStart!='None'|LEEnd!='None',T,F),
    IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F),
    IsGenicStart = ifelse(GeneStart!='',T,F),
    IsGenicEnd = ifelse(GeneEnd!='',T,F),
    Length = ifelse(as.character(ChrStart)!=as.character(ChrEnd)|Type=='INS'|ArmEnd!=ArmStart, -1, PosEnd-PosStart),
    #DoubleDupBE = ifelse(DupBEStart=='true'&DupBEEnd=='true',T,F),
    #SingleDupBE = ifelse(DoubleDupBE==0&(DupBEStart=='true'|DupBEEnd=='true'),T,F),
    #TICount = ifelse(LnkTypeStart=='TI',0.5,0)+ifelse(LnkTypeEnd=='TI',0.5,0),
    DBCount = ifelse(DBLenStart>=0,0.5,0)+ifelse(DBLenEnd>=0,0.5,0),
    #IsSglTI = ifelse(LnkTypeStart=='SGL',0.5,0),
    AsmbTICount = ifelse(AsmbMatchStart=='MATCH',0.5,0)+ifelse(AsmbMatchEnd=='MATCH',0.5,0),
    #InferTICount = TICount - AsmbTICount,
    #ShortTICount=ifelse(LnkTypeStart=='TI'&LnkLenStart<=1000,0.5,0)+ifelse(LnkTypeEnd=='TI'&LnkLenEnd<=1000,0.5,0),
    ClusterSize = ifelse(ClusterCount==1,'Single',ifelse(ClusterCount<=4,'Small','Large')),
    IsConsistent = ifelse(Consistency==0,T,F),
    IsChained = (ChainCount>=1),
    IsFoldBack = FoldbackLenStart>=0|FoldbackLenEnd>=0,
    RepeatedChainLink = (svData$ChainCount>0 & grepl(';',svData$ChainIndex)),
    IsPolyA = grepl('TTTTTTTTTT',InsertSeq)|grepl('AAAAAAAAAA',InsertSeq),
    IsLowQual = ResolvedType=='LowQual',
    CustomType=case_when(grepl('Chain',ResolvedType) ~ "Chain",grepl('Sgl',ResolvedType) ~ 'PairedSGL', grepl('TI',ResolvedType) ~ 'SyntheticDelDup',TRUE ~ ResolvedType)
  )

}

createBuckets <- function(cluster) {
  cluster %>% mutate(
    PloidyBucket=2**(pmin(7,pmax(-3,round(log(Ploidy,2),0)))),
    CnChEndBucket=2**(pmin(7,pmax(-3,round(log(pmax(AdjCNChgEnd,0.01),2),0)))),
    CnChStartBucket=2**(pmin(7,pmax(-3,round(log(pmax(AdjCNChgStart,0.01),2),0)))),
    ClusterCountBucket=2**(pmin(5,pmax(-3,round(log(ClusterCount,2),0)))),
    LnkLenStartBucket=ifelse(!LnkLenStart>0,0,2**(pmin(20,pmax(0,round(log(LnkLenStart,2),0))))),
    FoldBackLenStartBucket=ifelse(!FoldbackLenStart>0,0,2**(pmin(20,pmax(0,round(log(FoldbackLenStart,2),0))))),
    LnkLenEndBucket=ifelse(!LnkLenEnd>0,0,2**(pmin(20,pmax(0,round(log(LnkLenEnd,2),0))))),
    FoldBackLenEndBucket=ifelse(!FoldbackLenEnd>0,0,2**(pmin(20,pmax(0,round(log(FoldbackLenEnd,2),0))))),
    LengthBucket=ifelse(Type=='BND'|Type=='INS'|PosEnd-PosStart==0|ArmEnd!=ArmStart,0,2**(pmin(25,pmax(5,round(log(PosEnd-PosStart,2),0))))),
    HomLenBucket=2**(round(log(nchar(Homology),2),0)),
    NearestLenBucket=2**(round(log((NearestLen),2),0)),
    SynDelDupLenBucket=4**(round(log((SynDelDupLen),4),0)),
    SynDelDupTILenBucket=4**(round(log((SynDelDupTILen),4),0)),
    InsLenBucket=2**(round(log(nchar(InsertSeq),2),0)),
    QualScoreBucket=2**(round(log(nchar(QualScore),2),0))
  )
}

createsvLinksBuckets <- function(svLinks) {
  svLinks %>% mutate(
    TILengthBucket=ifelse(TILength==0,0,2**(pmin(25,pmax(5,round(log(TILength,2),0))))),
    synDelDupLengthBucket=ifelse(SynDelDupLen==0,0,2**(pmin(25,pmax(5,round(log(SynDelDupLen,2),0))))),
    DBLenEndBucket=ifelse(DBLenEnd==0,0,2**(pmin(25,pmax(5,round(log(DBLenEnd,2),0))))),
    DBLenStartBucket=ifelse(DBLenStart==0,0,2**(pmin(25,pmax(5,round(log(DBLenStart,2),0))))),
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

createClusterArmStats <- function(cluster) {
  rbind(cluster %>% mutate(Arm=ArmStart,Chr=ChrStart),svData %>% filter(ChrEnd!=0) %>% mutate(Arm=ArmEnd,Chr=ChrEnd)) %>%
          group_by(SampleId,ClusterId) %>% mutate(clusterSGLCount=sum(ChrEnd==0)) %>% ungroup() %>%
          group_by(SampleId,ClusterId,Arm,Chr,ResolvedType,ClusterBreakends=ClusterCount*2-clusterSGLCount) %>% 
          summarise(LocalCount=sum(ArmEnd==ArmStart&ChrEnd==ChrStart),CrossArmCount=sum(ArmEnd!=ArmStart&ChrEnd==ChrStart),SGLCount=sum(ChrEnd==0),RemoteCount=sum((ArmEnd!=ArmStart|ChrEnd!=ChrStart)&ChrEnd!=0)) %>%
          mutate(expectedLocalProportion = round((LocalCount+RemoteCount+SGLCount)/ClusterBreakends,4),actualLocalProportion=(round((LocalCount+SGLCount)/(LocalCount+RemoteCount+SGLCount),4)),diff=actualLocalProportion-expectedLocalProportion) %>%
          mutate(pValue=pmin(ppois(LocalCount+SGLCount,(LocalCount+RemoteCount+SGLCount)*expectedLocalProportion,FALSE,FALSE),ppois(LocalCount+SGLCount,(LocalCount+RemoteCount+SGLCount)*expectedLocalProportion,TRUE,FALSE)))
}

createSampleSummary <- function(cluster) {
 cluster %>% 
   group_by(SampleId,ClusterId,ResolvedType,Type=ifelse(ClusterCount==1,as.character(Type),'')) %>% 
   summarise(n=n(),countFS=sum(IsFS),
             countSimpleDUPShort=sum(Type=='DUP' & PosEnd-PosStart<1e5),
             countSimpleDUPLong=sum(Type=='DUP' & PosEnd-PosStart>=1e5),
             countDUP_Ext_TIShort=sum(ResolvedType=='DUP_Ext_TI'&SynDelDupLen<1e5),
             countDUP_Ext_TILong=sum(ResolvedType=='DUP_Ext_TI'&SynDelDupLen>=1e5)) %>% group_by(SampleId) %>% 
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
              countSimpleDUPShort=sum(countSimpleDUPShort),
              countSimpleDUPLong=sum(countSimpleDUPLong) ,
              countSimpleBND=sum(Type=='BND'),
              countSimpleINV=sum(Type=='INV'),
              countSimpleSGL=sum(Type=='SGL'),
              countSimpleNONE=sum(Type=='NONE'),
              countDEL_Ext_TI=sum(ResolvedType=='DEL_Ext_TI'),
              countDEL_Int_TI=sum(ResolvedType=='DEL_Int_TI'),
              countDUP_Ext_TI=sum(ResolvedType=='DUP_Ext_TI'),
              countDUP_Ext_TIShort=sum(countDUP_Ext_TIShort),
              countDUP_Ext_TILong=sum(countDUP_Ext_TILong),
              countDUP_Int_TI=sum(ResolvedType=='DUP_Int_TI'),
              countRecipInv=sum(ResolvedType=='ReciprocalInversion'),
              countRecipTrans=sum(ResolvedType=='RecipTrans'),
              #countShortTI=sum(ShortTICount),
              countFS=sum(countFS))
}

#############################################
################ LOADING ####################
#############################################

PATH='~/Dropbox/HMF Australia team folder/Structural Variant Analysis/'
PATH='~/hmf/analyses/SVAnalysis/'

svData = read.csv(paste(PATH,'SVA_SVS.csv',sep=''), header = T, stringsAsFactors = F)
svClusters = (read.csv(paste(PATH,'SVA_CLUSTERS.csv',sep=''), header = T, stringsAsFactors = F))
svDrivers = read.csv(paste(PATH,'SVA_DRIVERS.csv',sep=''), header = T, stringsAsFactors = F)
svLinks = read.csv(paste(PATH,'SVA_LINKS.csv',sep=''), header = T, stringsAsFactors = F)
sampleCancerTypes= (read.csv(paste(PATH,'sample_cancer_types.csv',sep=''), header = T, stringsAsFactors = F))
svChordStatus = read.csv(paste(PATH,'SVChordStatus.csv',sep=''), header = T, stringsAsFactors = F)
svDriverCatalog = read.csv(paste(PATH,'SVDriverCatalog.csv',sep=''), header = T, stringsAsFactors = F)
svGermline = read.csv(paste(PATH,'SVGermline.csv',sep=''), header = T, stringsAsFactors = F)
svDriverAndGermline=(rbind(svDriverCatalog %>% select(gene,sampleId,driver,driverLikelihood), svGermline %>% mutate(driverLikelihood=1,driver=ifelse(biallelic==1,'Germline Biallleic','Germline Monoallleic')) %>% select(gene,sampleId,driverLikelihood,driver)))
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
svDrivers = merge(svDrivers, sampleCancerTypes, by='SampleId', all.x=T)
sampleList = svData %>% distinct(SampleId,CancerType)
#Annotations
svData=createBuckets(svData)
svData = svData %>% mutate(ResolvedType=ifelse(ResolvedType=='DEL_Int_TI'&Type=='INV'&SynDelDupTILen>=0.5*SynDelDupLen,'ReciprocalInversion',ResolvedType))  # TEMP DEFINITION - CAN DO BETTER
svData=sv_set_common_fields(svData)
svSampleSummary=createSampleSummary(svData)
foldbacks=createFoldbacks(svData)
dbData = rbind(svData %>% filter(DBLenStart>-31) %>% mutate(DBLength = DBLenStart,Assembled = ifelse(AsmbMatchStart=="MATCH","Assembled","NotAssembled"),LE=LEStart,RefContext=RefContextStart,Orient=OrientStart),
               svData %>% filter(DBLenEnd>-31) %>% mutate(DBLength = DBLenEnd, Assembled = ifelse(AsmbMatchEnd=="MATCH","Assembled","NotAssembled"),LE=LEEnd,RefContext=RefContextEnd,Orient=OrientEnd)) %>%
               mutate(DBLenBucket = ifelse(DBLength==0,0,ifelse(DBLength<0,-(2**round(log(-DBLength,2))),2**round(log(DBLength,2)))))
svLinks=createsvLinksBuckets(svLinks)

#############################################
########## OVERVIEW OF COHORT ###############
#############################################

# 1. By RESOLVED Type
View(svData %>% filter(!IsLowQual) %>% group_by(ResolvedType,Type) %>% tally() %>% spread(Type,n))

# 2. By CLUSTERCOUNTBUCKET
View(svData %>% filter(!IsLowQual) %>% group_by(ClusterCountBucket,Type) %>% tally() %>% spread(Type,n))

#3. Simple SV make up 70% of clusters but only 30% of SVs
View(svData %>% filter(!IsLowQual) %>% group_by(CustomType,ClusterId,SampleId) %>% tally() %>% group_by(CustomType) %>% summarise(variants=sum(n),clusters=n()))


##############################################################
########## Short Deletion Bridge Analysis ####################
#############################################################


#TEMP...
plot_count_by_bucket_and_type(cohortSummary(dbData %>% filter(ResolvedType=='ReciprocalInversion'),
                                            "SynDelDupLen<=10000,!IsLowQual,!(CustomType %in% c('SimpleSV','PairedSGL')),!IsPolyA",'SynDelDupLen,ClusterDesc'),
                              'SynDelDupLen','ClusterDesc','SyntheticDelDupLength',useLogX = F,useLogY = F)
plot_count_by_bucket_and_type(cohortSummary(dbData %>% filter(ResolvedType=='DUP_Ext_TI'),
                                            "SynDelDupLen<=30,!IsLowQual,!(CustomType %in% c('SimpleSV','PairedSGL'))",'SynDelDupLen,IsPolyA'),
                              'SynDelDupLen','IsPolyA','SyntheticDelDupLength',useLogX = F,useLogY = F)

#1.DBLength by ResolvedType => less than 50 bases (NB - DBLength of 1 means exact break - shoould correct this)
# Note the LINE double peak and the sharp feature for Reciprocal Inversion, Reciprocal Translocations and None
plot_count_by_bucket_and_type(cohortSummary(dbData %>% mutate(CustomType=ifelse(CustomType=='Chain'&Length<1e3&Type=='INV','ShortINVinChain',CustomType)),
                            "DBLength<=30,!IsLowQual,!(CustomType %in% c('SimpleSV','PairedSGL'))",'DBLength,CustomType'),
                            'DBLength','CustomType','DBLengthByCustomType(<50bases)',useLogX = F,useLogY = F)


#2. DB Length - Exact break short INV feature also exists in chains but only when INV length is <1e3   
# Occurs prinicipally in the most complex clusters.  Short Inversions are much more likely than other variant types than other short variants.
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length>1e3,CustomType=='Chain'",'DBLength,ClusterCountBucket'),'DBLength','ClusterCountBucket','DBLength for short variants in chains by cluster count',useLogX = F,useLogY = F)
# Is not primarily a foldback feature
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length<1e3,CustomType=='Chain'",'DBLength,IsFoldBack'),'DBLength','IsFoldBack','DBLength for short variants in chains by isFoldback',useLogX = F,useLogY = F)
# Mostly around 60-500 bases.
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length<1e3,CustomType=='Chain'",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for short variants in chains by length bucket',useLogX = F,useLogY = F)

#3. 'NONE' clusters - appear to have characteristics of short INVs and likely LINE elements for BND + SGL
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,ResolvedType=='None'",'DBLength,LengthBucket=Length<0|Length<1e3'),'DBLength','LengthBucket','DBLength for Complex Chains ByLengthBucket>1e3(<50bases)',useLogX = F,useLogY = F)

#4. LINE Elements - no clear indication of why there are 2 peaks.
#The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific.   Also appars to be the same regardless of source element
print(ggplot(data = dbData %>% filter(DBLength<=50,ResolvedType %in% c('Line')) %>% group_by(CancerType,OLPeak=DBLength<(-7)) %>% tally(), 
             aes(x = reorder(CancerType, -n), y = n, fill =OLPeak))
      + geom_bar(stat = "identity", colour = "black") + ylab("DB Count") + xlab("Tumor Type") + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15)) + labs(title = "DB Counts for LINE Elements by cancerType and Overlapping Peak"))
# Also does not seem to depend significantly on refContext.  
View(dbData %>% filter(ResolvedType=='Line',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context,longOverlap=DBLength<(-7)) %>% tally() %>% arrange(-n) %>% spread(longOverlap,n))

#######################################################
########## LINE specific properties ###################
#######################################################

# 1. INSERTION MOTIF: ~15% are A-TTTTT with another 30% or so a similar variant
View(dbData %>% filter(ResolvedType=='Line',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context) %>% tally() %>% arrange(-n))

View(dbData %>% filter(ResolvedType=='Line',LE=='None',IsPolyA,RefContext!='',Orient==-1) %>% select(InsertSeq,RefContext,DBLength,everything()))

# 2. CLUSTER COUNT vs TOTAL COUNT Count of clusters per sample is very related to the total number of variants (averaging at 3-4 per cluster)
scatterPlot(svClusters %>% filter(ResolvedType=='Line') %>% group_by(SampleId) %>% summarise(countLine=n(),sumLine=sum(ClusterCount)),'countLine','sumLine',F,F)

# 3. LINE Clusters affect a predictable number of chromosome arms that grows with the svCount in the cluster
scatterPlot(svClusters %>% filter(ResolvedType=='Line'),'ClusterCount','ArmCount',F,F)
View(svClusters %>% filter(ResolvedType=='Line'))
# 4. LINE elements mmuch more common in Esophagus and Stomach cancers;  
ggplot(data=(merge(sampleList,svClusters %>% filter(ResolvedType=='Line',!(ClusterDesc %in% c('SGL','SGL=2','INS'))) %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>%
        group_by(CancerType) %>% mutate(TypeCount=n()) %>% filter(TypeCount>20) %>% ungroup())) +
        stat_ecdf(aes(n), geom = "step", pad = FALSE) + facet_wrap(~CancerType) + labs(title = 'CDF # of LINE Clusters per Sample') +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.grid.minor.x = element_blank(),
         legend.position="none", legend.title = element_blank()) +  coord_flip()
# VIOLIN equivalent
#ggplot(merge(sampleList,svClusters %>% filter(ResolvedType=='Line',!(ClusterDesc %in% c('SGL','SGL=2','INS'))) %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>%
#         group_by(CancerType) %>% mutate(TypeCount=n()) %>% filter(TypeCount>20) %>% ungroup() , aes(CancerType, n)) + 
#  geom_violin(scale="width",fill="#6baed6") +   theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) + ggtitle("") + xlab("") + ylab("Line Elemnts Per Sample") + 
#  theme(axis.ticks = element_blank(), legend.position="bottom") + scale_y_continuous(expand=c(0.01, 0.01)) + scale_fill_manual(values = cancerTypeColours) + theme(legend.position="none") +coord_flip() 

# 5. Known and suspected LINE element source locations => many potentially novel source locations that feature
View(rbind(svData %>% filter(ResolvedType=='Line') %>% filter(LEStart!='None') %>% group_by(chr=ChrStart,pos=round(PosStart,-4),SampleId) %>% summarise(n=n(),countSuspected=sum(LEStart=='Suspect'),
           countKnown=sum(LEStart=='Known'),countPOLYA=sum(IsPolyA&Type=='BND'),countSuspectKnown=sum(LEStart=='Known;Suspect')),
           svData %>% filter(ResolvedType=='Line') %>% filter(LEEnd!='None') %>% group_by(chr=ChrEnd,pos=round(PosEnd,-4),SampleId) %>% summarise(n=n(),countSuspected=sum(LEEnd=='Suspect'),
           countKnown=sum(LEEnd=='Known'),countPOLYA=sum(IsPolyA&Type=='BND'),countSuspectKnown=sum(LEEnd=='Known;Suspect'))) %>%
       group_by(chr,pos,SampleId) %>% summarise(n=sum(n),countSuspected=sum(countSuspected),countSuspectKnown=sum(countSuspectKnown),countKnown=sum(countKnown),countPOLYA=sum(countPOLYA)) %>% 
       group_by(chr,pos) %>% summarise(n=sum(n),countSuspected=sum(countSuspected),countSuspectKnown=sum(countSuspectKnown),countKnown=sum(countKnown),countPOLYA=sum(countPOLYA),numSamples=n()) )

#######################################################
########## Pseudogene insertions ###################
#######################################################

svPseudoGeneLinks = svLinks %>% filter(ExonMatch!='') %>% separate(ExonMatch,c('transcriptId','exonRank','exonLength'))

View(svLinks %>% group_by(Id1,Id2) %>% tally() %>% filter(n>1) )

#1. 440 events
View(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd,transcriptId) %>% summarise(count=n(),minExon=min(as.numeric(exonRank)),maxExon=max(as.numeric(exonRank)),range=maxExon-minExon))

#2. found in 6% of samples - nearly all with highly activated LINE elements
View(merge(svLinks %>% filter(ExonMatch!='') %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd) %>% count() %>% group_by(SampleId) %>% count(),svSampleSummary %>% select(SampleId,count,countLine,sumLine),by='SampleId',all=T,fill=0))

# TO DO:  Create circos charts with gene transcripts for further analysis

################################################
############## FOLDBACKS #######################
################################################

#1.Foldback length distribution for simple inversions short chained foldbacks
print(ggplot(data = foldbacks %>%  filter(!IsLowQual)  %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=INV, colour='INV')) + geom_line(aes(y=Combo, colour='Combo ( FB Chain Length < 5k)' )) + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))

#2. Combo foldbacks counts increase with simple foldback counts
scatterPlot(foldbacks %>%  filter(!IsLowQual) %>% group_by(SampleId) %>% summarise(simpleFoldbackCount=sum(FoldbackType=='INV'),comboFoldbackCount=sum(FoldbackType!='INV'&FoldbackChainLengthBucket<5000)),
            'simpleFoldbackCount','comboFoldbackCount',F,F)

############################################################
########### Reciprocal Inversions ##########################
############################################################

# #1. TO DO:  show the length distribution and come up with a standardised definition
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='ReciprocalInversion'",'CancerType,LengthBucket') 
                              ,'LengthBucket','CancerType','Reciprical INV',useLogY =F )

#2. Reciprocal INV NOT enriched in typical DEL lengths OR FS!
plot_count_by_bucket_and_type(cohortSummary(svData,"!IsLowQual,ResolvedType %in% c('SimpleSV','ReciprocalInversion')",'IsFS,LengthBucket'),'LengthBucket','IsFS','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )


############################################################
########### Reciprocal Translocations ######################
############################################################

# TO DO:  show the length distribution and come up with a standardised definition
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='RecipTrans'",'CancerType,NearestLenBucket') 
                              ,'NearestLenBucket','CancerType','Reciprical Trans',useLogY =F )
View(svData %>% filter(ClusterDesc=="BND=2") %>% group_by(ResolvedType) %>% tally()) 

# TO DO: 
View(svData %>% filter(ClusterDesc=="BND=2",ArmCount==2,Consistency==0,ResolvedType=='None')) 

############################################################
########### Short TIs #####################################
############################################################
#0. Numbers per sample
ggplot(data=merge(svData %>% group_by(SampleId) %>% tally(),svLinks %>% 
    group_by(SampleId,shortTI=ifelse(TILength<1000,'ShortTICount','LongTICount')) %>% tally() %>% spread(shortTI,n),by='SampleId')) +
      geom_point(aes(n,ShortTICount))

#1. TIs length distibution, showing 2 peaks for most types ('short', ie <1k and long -> presumed to be cause by random length of distance between breakages)
ggplot(data=svLinks %>% group_by(TILengthBucket,ResolvedType) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_x_log10()+ facet_wrap(~ResolvedType)

ggplot(data=svLinks  %>% filter(TILength<30000) %>% mutate(TILengthBucket=round(TILength+0.001,-2)) %>% group_by(TILengthBucket) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_y_log10()
svLinks  %>% filter(TILength<500) %>% mutate(TILengthBucket=round(TILength,-1)) %>% filter(TILengthBucket==100) %>% summarise(count=n())
svLinks %>% group_by(SampleId,shortTI=TILength<1000) %>% summarise(count=n()) %>% spread(shortTI,n)
ggplot(data=svLinks %>% group_by(SampleId,shortTI=TILength<1000) %>% summarise(count=n()) %>% spread(shortTI,count))  +
    geom_point(aes(`FALSE`,`TRUE`))
ggplot(data=svLinks %>% group_by(SampleId,shortTI=TILength<1000) %>% summarise(count=n())) + 
  stat_ecdf(aes(`count`,color='count'),geom = "step", pad = FALSE) + facet_wrap(~shortTI)
    

ggplot(data=svData %>% filter(ResolvedType=='SimpleSV',PosEnd-PosStart>20000,PosEnd-PosStart<500000,Type=='DEL') %>% group_by(SampleId,IsFS) %>% tally() %>% spread(IsFS,n,fill=0)) + 
  stat_ecdf(aes(`TRUE`,color='isFS'),geom = "step", pad = FALSE) + stat_ecdf(aes(`FALSE`,color='FSLength non-FS'),geom = "step", pad = FALSE) + 
  ylim(0,1) + labs(title = 'CDF # of FS length DELs')

#2. Short TIs create 'synthetic' events that match the simpleSV event profile per sample
scatterPlot(svSampleSummary,'countSimpleDEL','countDEL_Ext_TI',F,F) # DEL correlated with synthetic DEL
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_Ext_TI',F,F)  # DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_Int_TI',F,F)  # Int DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPShort','countDUP_Ext_TIShort',F,F)  # Short DUP correlated with short synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPLong','countDUP_Ext_TILong',F,F)  # Long DUP correlated with long synthetic DUP

##############################################################
########## FRAGILE SITE ######################################
##############################################################

#1. Fragile Site Length Distribution
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType %in% c('SimpleSV','ReciprocalInversion')","LengthBucket,IsFS"),'LengthBucket','IsFS','Length Distribution of SimpleSV & Reciprocal Inversion by IsFS',useLogY =F)

#2. FS like length dels outside FS are correlated, but not all 20k-500k are FS Dels
scatterPlot(svData %>% filter(Type=='DEL'&ResolvedType=='SimpleSV') %>% group_by(SampleId) %>% summarise(countFS=sum(IsFS),countFSLike=sum(!IsFS&PosEnd-PosStart>2e4&PosEnd-PosStart<5e5)),'countFS','countFSLike',F,F) 

#3. Extreme tail possible for FS length DELS => Unrelated process
ggplot(data=svData %>% filter(ResolvedType=='SimpleSV',PosEnd-PosStart>20000,PosEnd-PosStart<500000,Type=='DEL') %>% group_by(SampleId,IsFS) %>% tally() %>% spread(IsFS,n,fill=0)) + 
  stat_ecdf(aes(`TRUE`,color='isFS'),geom = "step", pad = FALSE) + stat_ecdf(aes(`FALSE`,color='FSLength non-FS'),geom = "step", pad = FALSE) + 
  ylim(0,1) + labs(title = 'CDF # of FS length DELs')

#4. Location of FS + potentially novel FS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL'|Type=='DUP',PosEnd-PosStart>2e4,PosEnd-PosStart<5e5) %>% group_by(chr=ChrStart,pos=round(PosStart,-6), IsFS) %>% summarise(n=n()) %>% spread(IsFS,n,fill=0) %>% mutate(Total=`TRUE`+`FALSE`))

#5. Types of events in FS
View(svData %>% group_by(IsFS,ResolvedType,Type) %>% tally() %>% spread(IsFS,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))

#6. TO DO: CHECK GENIC vs non-GENIC and length of GENES
ggplot(data=svData %>% mutate(gene=pmax(GeneStart,GeneEnd)) %>% filter(IsFS,Type=='DEL',PosEnd-PosStart<1e6) %>% group_by(gene) %>% mutate(count=n()) 
       %>% ungroup() %>% filter(count>20),aes(PosStart,PosEnd))+geom_point()+facet_wrap(~gene,scales ='free')

View(svData %>% group_by(driver) %>% count())
#################################################
################# COMPLEX CLUSTER STATS #########
#################################################

#1 LOCAL - REMOTE CONSISTENCY
clusterArmStats = createClusterArmStats(svData)
View(clusterArmStats %>% filter(LocalCount+RemoteCount+SGLCount>10))

# 2. Per arm by ClusterCount
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='Line') %>% group_by(SampleId,chr_arm,CC=pmin(ClusterCount,8),id=ClusterId) %>% tally(),
           svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='Line') %>%group_by(SampleId,chr_arm,CC=pmin(ClusterCount,8),id=ClusterId) %>% tally()) %>% 
       group_by(SampleId,chr_arm,CC,id) %>% tally() %>%
       group_by(SampleId,CC,chr_arm) %>% tally() %>% spread(CC,n) %>%
       filter(chr_arm!='0_P'))

#3. LARGE Clusters - # arms not super correlated with # of SVs
scatterPlot(svClusters %>% filter(ClusterCount>5,ResolvedType!='LINE') %>% mutate(realArmCount=ArmCount-FragmentArms),'realArmCount','ClusterCount',F,F)

#4. CDF of large non-line clusters per sample
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>10,ResolvedType!='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>10 Per Sample')
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>50,ResolvedType!='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>50 Per Sample')

##############################################################
########## SIMPLE SV SIGNATURES ##############################
##############################################################

#0. Counts per sample 
# TO DO: MOVE TO VIOLIN PLOT
ggplot(data=svData %>% filter(ResolvedType=='SimpleSV',Type %in% c('DEL','DUP')) %>% group_by(SampleId,Type,CancerType) %>% tally() %>% spread(Type,n,fill=0)) + 
  stat_ecdf(aes(DUP,color='DUP'),geom = "step", pad = FALSE) + ylim(0,1) + stat_ecdf(aes(DEL,color='DEL'),geom = "step", pad = FALSE) +
      labs(title = 'CDF # of SimpleSV per Sample') + scale_x_log10() + facet_wrap(~CancerType)

#1. top50 simple
filter = cohortSummary(svData,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV')",'SampleId') %>% top_n(50,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','AllDel&Dup',useLogY =F)

#2. SHORT DELS: <100 bases length is related 1k to 15k bucket
scatterPlot(svData %>% filter(Type=='DEL',ResolvedType=='SimpleSV') %>% group_by(SampleId) %>% summarise(DEl_LT_100_bases=sum(PosEnd-PosStart<100),DEL_1k_to_5k_bases=sum(PosEnd-PosStart>1e3&PosEnd-PosStart<0.5e4)),
            'DEl_LT_100_bases','DEL_1k_to_5k_bases',F,F) 

#3. ULTRA SHORT DUP: Top 20:
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length<100",'SampleId') %>% top_n(20,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Ultra Short Dup',useLogY =F)
#NB - Colorectal cancers dominate this list (treatment?)
View(svData %>% filter(Type=='DUP',Length<100,Length>0,ResolvedType=='SimpleSV') %>% group_by(SampleId,CancerType) %>% tally() %>% arrange(-n))
ggplot(data=svData %>% filter(Type=='DUP',Length<100,Length>0,ResolvedType=='SimpleSV') %>% group_by(SampleId,CancerType) %>% tally(),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of Foldback Clusters per Sample') + facet_wrap(~CancerType)


#4.ULTRA SHORT DEL: Top 20:
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length<100",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Ultra Short Del',useLogY =F)


#5. SHORT DEL: Top 30:
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length>1000,Length<10000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Short Del',useLogY =F)

#5. BRCA Length DUP (5k-100k): Top 30:
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>5000,Length<100000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','BRCA LENGTH DUP',useLogY =F)

#5. CDK12 Length DUP (100k+): Top 30:
## NB - there is a 
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>100000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','CDK12 LENGTH DUP',useLogY =F)


# GENE ENRICHMENT AND SIMPLE SIGNATURES TP53 enriched in all cases!

#CDK12 AND CCNE1 associated with 2 different Long DUP types.  TP53 also enriched
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>100000",'SampleId,CancerType') %>% top_n(40) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('CDK12','CCNE1','ANKRD11')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))# %>% group_by(gene) %>% tally())

# BRCA1 AND ANKRD11 associated with midrange DUPS associated with 2 different Med DUP types
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length<100000,Length>1000",'SampleId,CancerType') %>% top_n(100,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('BRCA1','ANKRD11')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))# %>% group_by(gene) %>% tally())

# BRCA2,BRCA1 AND ANKRD11 strongly associated with very short dels
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length<1000",'SampleId,CancerType') %>% top_n(100,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('BRCA2','ANKRD11','BRCA1')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))

# TP53 and BRCA2 STRONGLY positively associated with 1k-10k dels ()
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length>1000,Length<10000",'SampleId,CancerType') %>% top_n(50,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('TP53','BRCA2')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))
View(merge((filter),svDriverAndGermline %>% filter(),by.x='SampleId',by.y='sampleId',all.x=T) %>% group_by(gene) %>% summarise(n=n(),meanDL=mean(driverLikelihood)) %>% arrange(-n))

# TP53 STRONGLY positively associated with 10k-500k dels() (especially if IsFS but also regardless)
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length>20000,Length<500000,IsFS",'SampleId,CancerType') %>% top_n(200,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('TP53')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length>20000,Length<500000,!IsFS",'SampleId,CancerType') %>% top_n(200,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('TP53')),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))

# ANKRD11 drivers have with a diverse range of effects, is mutually exclusive to BRCA1 & BRCA2 (but may just be a weak site of some description?)
filter = cohortSummary(merge(svData,svDriverAndGermline %>% filter(gene=='ANKRD11',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(50,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','ANKRD11 samples',useLogY =F)

# SPOP causes 1 1-10k peak but no strong LT100 peak in prostate cancer only
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV',CancerType=='Prostate'),svDriverAndGermline %>% filter(gene=='SPOP',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(40,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'|ResolvedType=='ReciprocalInversion')",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','SPOP Top 30',useLogY =F)

# PTPRN2 causes BRCA like signature in Breast cancers only
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV'),svDriverAndGermline %>% filter(gene=='PTPRN2',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'|ResolvedType=='ReciprocalInversion')",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','PTPRN2 Top 30',useLogY =F)

# CDK12 causes 2 peaks whereas CCNE1 only causes 1
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV'),svDriverAndGermline %>% filter(gene=='CDK12',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','CDK12 Top 30',useLogY =F)
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV'),svDriverAndGermline %>% filter(gene=='CCNE1',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','CCNE1 Top 30',useLogY =F)

# BRCA1 to p 30
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV'),svDriverAndGermline %>% filter(gene=='BRCA1',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','BRCA1 Top 30',useLogY =F)

# BRCA2 top 30
filter = cohortSummary(merge(svData %>% filter(ResolvedType=='SimpleSV'),svDriverAndGermline %>% filter(gene=='BRCA2',driverLikelihood>0.8) %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','BRCA2 Top 30',useLogY =F)


# PTCH1 in Basal Cell Carcinoma /LINE correlation
View(merge(svData %>% filter(CancerType=='Skin',ResolvedType=='Line') %>% group_by(SampleId,cancerSubtype) %>% count(),svDrivers %>% filter(Gene=='PTCH1') %>% distinct(SampleId,DriverType) %>% mutate(PTCH1=TRUE),by='SampleId',all.x=T))


# TP53 highly enriched in LINE
filter = cohortSummary(svData,"IsFS,CancerType=='Breast'",'SampleId,CancerType') %>% top_n(100,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDriverAndGermline %>% filter(gene %in% c('TP53',"MDM2")),by.x='SampleId',by.y='sampleId',all.x=T) %>% arrange(-count))

############################################################
#############Driver vs Signature correlation ###############
############################################################
driverSignatureCoocurrence = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/coc_sv_results_extended.csv') 
View(driverSignatureCoocurrence %>% arrange(FDR) %>% filter(CancerType!='All',Category!='SIMPLE_SV',FDR<0.1))
write.csv(driverSignatureCoocurrence %>% arrange(FDR) %>% filter(CancerType!='All',Category!='SIMPLE_SV',FDR<0.1),'~/driverSignatureCoocurrence.csv')

cancerTypeLOHCoocurrence = read.csv('~/hmf/analyses/SVAnalysis/coc_sv_loh_cancer.csv') 
View(cancerTypeLOHCoocurrence %>% group_by(CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))
write.csv(cancerTypeLOHCoocurrence %>% group_by(CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR),'~/cancerTypeLOHCoocurrence.csv')
driverGeneLOHCoocurrence = read.csv('~/hmf/analyses/SVAnalysis/coc_gene_loh_results_extended.csv') 
View(driverGeneLOHCoocurrence %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,positive=CountGtExp,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))
write.csv(driverGeneLOHCoocurrence %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR),'~/driverGeneLOHCoocurrence.csv')

##########################################################
########### Enrichment by location #######################
##########################################################

#1. Enriched segments of short DEL / DUP
#### NB - Replication timing average needs to remove 0s from the mean
#SHORT DELS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% summarise(n=n(),avgRep=mean(RepOriginStart)) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))

temp=(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4,!IsFS) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot(data=temp,aes(pos,n)) + geom_point() + facet_wrap(~chr)

#FS DELS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4,PosEnd-PosStart<5e5) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows)%>% summarise(n=n(),fs=sum(ifelse(IsFS,1,0))) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
temp=(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL',PosEnd-PosStart<2e4,PosEnd-PosStart<5e5) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot(data=temp,aes(pos,n)) + geom_point() + facet_wrap(~chr)


#SHORT DUP
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart<1e3,!IsFS) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
temp=(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart<1e3,!IsFS) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot(data=temp,aes(pos,n)) + geom_point() + facet_wrap(~chr)
#BRCA DUPS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e3,PosEnd-PosStart<1e5) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-5),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e5,FALSE),q=p.adjust(p,"BH",2.8e9/1e5)))
temp=(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>5e3,PosEnd-PosStart<8e4) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-5),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e5,FALSE),q=p.adjust(p,"BH",2.8e9/1e5)))
ggplot(data=temp,aes(pos,n)) + geom_point() + facet_wrap(~chr)
#CDK12 DUPS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e5,PosEnd-PosStart<1e6) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)) %>% filter(q<0.01))
temp=(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e5,PosEnd-PosStart<1e6) %>% mutate(nrows=n()) %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% tally() %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot(data=temp,aes(pos,n)) + geom_point() + facet_wrap(~chr)
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DUP',PosEnd-PosStart>1e5,!IsFS,ChrStart==11,round(PosStart,-6)==6.5e7))

###################################################
########### Driver Variants #######################
###################################################
driverVariants=merge(svDrivers,svData %>% select(SampleId,SvId=Id,ResolvedType,ClusterCount,Length,LengthBucket,Ploidy),by=c('SvId',"SampleId"),all.x=T) %>%  unite(DriverCombined,DriverType,MatchInfo) %>% mutate(DriverCombined=ifelse(grepl('AMP',DriverCombined),'AMP',DriverCombined)) 

View(svDrivers %>% filter(DriverType=='BIALLELIC',MatchInfo==''))
# 1. By resolvedType
View(driverVariants %>% group_by(Gene,GeneType,SampleId,ClusterId,DriverCombined,ResolvedType) %>% count() %>% group_by(Gene,GeneType,ResolvedType) %>% count() %>% spread(ResolvedType,nn))
View(driverVariants %>% group_by(Gene,GeneType,SampleId,ClusterId,DriverCombined,Chromosome) %>% count() %>% group_by(Chromosome,Gene,GeneType,DriverCombined) %>% count())

View(driverVariants %>% group_by(Gene,GeneType,SampleId,ClusterId,DriverCombined,CancerType) %>% count() %>% 
       group_by(Gene,GeneType,DriverCombined,CancerType) %>% count() %>% spread(CancerType,nn))
View(driverVariants %>% group_by(Gene,GeneType,SampleId,ClusterId,DriverCombined,CancerType) %>% count() %>% 
       group_by(GeneType,DriverCombined,CancerType) %>% count() %>% spread(CancerType,nn))

View(driverVariants %>% filter(DriverCombined!='DNDS_') %>% group_by(Gene,GeneType,SampleId,ClusterId,DriverCombined) %>% count() %>% 
       group_by(Gene,GeneType) %>% count())

View(svDrivers %>% filter(GeneType=='TSG') %>% group_by(Gene,DriverType) %>% count() %>%  spread(DriverType,n))
       group_by(Gene,GeneType,DriverCombined,CancerType) %>% count() %>% spread(DriverCombined,nn))

View(driverVariants %>% filter(ResolvedType=='SimpleSV') %>% group_by(Gene,DriverCombined,LengthBucket) %>% count() %>% spread(LengthBucket,n))

##### TESTING
View(driverVariants %>% filter(Gene=='MYC') %>% group_by(SampleId,ClusterCount,ClusterId) %>% summarise(n=n(),maxPloidy=max(Ploidy)))
View(driverVariants %>% filter(DriverCombined=='AMP') %>% group_by(SampleId,ClusterCount,ClusterId) %>% summarise(n=n(),maxPloidy=max(Ploidy)))

#View(driverVariants %>% filter(DriverGene=='MYC'))
#View(svData %>% filter(!(DriverStart==''),!(DriverEnd==''),DriverStart!=DriverEnd) %>% select(DriverStart,DriverEnd))
#View(svData %>% filter(grepl('TP53',DriverStart)|grepl('TP53',DriverEnd)))

###################################################
########### HIGH LEVEL AMPLIFICTATION #############
###################################################

# 566/3582 (16%) of  samples have a very high copy nummber
View(svData %>% filter(Ploidy>20,AdjCNChgStart>20,Type!='SGL',Type!='NONE') %>% group_by(SampleId,isSimpleSV=ResolvedType=='SimpleSV',ClusterId) %>% tally() %>% group_by(SampleId,isSimpleSV) %>% count() %>% arrange(-nn) %>% spread(isSimpleSV,nn))

#TEMP
View(svData %>% filter(Ploidy>20,AdjCNChgStart>20,Type!='SGL',Type!='NONE') %>% filter(SampleId==chartSample))


# ~50%  of those involve a single CHR
View(svData %>% filter(Type!='SGL',Type!='NONE') %>% group_by(SampleId,ClusterId) %>% mutate(maxClusterPloidy=max(Ploidy),minChr=min(ChrStart),maxChr=max(ChrEnd)) %>% ungroup() %>% filter(maxClusterPloidy>20,AdjCNChgStart>20,minChr==maxChr) 
     %>% group_by(SampleId,CancerType,minChr) %>% count() %>% spread(minChr,n))

# By ploidy bucket       
View(svData %>% group_by(SampleId,ClusterId) %>% mutate(maxClusterPloidy=max(Ploidy),minChr=min(ChrStart),maxChr=max(ChrEnd)) %>% ungroup() %>% filter(maxClusterPloidy>20) %>% 
       group_by(SampleId,CancerType,ClusterId,minChr,maxChr,sameChr=minChr==maxChr,ResolvedType,ClusterDesc,PloidyBucket) %>% tally() %>% spread(PloidyBucket,n))

###################################################
########### DOUBLE MINUTE ########################
###################################################
View(svClusters %>% filter(grepl('DM',Annotations)))

View(merge(svClusters %>% filter(grepl('DM',Annotations),!grepl('DM_Unclear',Annotations)) %>% select(SampleId,ClusterId,Annotations),svData,by=c('SampleId','ClusterId')) %>% filter(ClusterCount==2))           
View(merge(svClusters %>% filter(grepl('DM',Annotations),!grepl('DM_Unclear',Annotations)) %>% select(SampleId,ClusterId,Annotations),svData,by=c('SampleId','ClusterId')) %>% group_by(Type,SampleId,ClusterId) %>% count() %>% spread(Type,n))

########################################
########## GENIC #######################
########################################

# TO DO: REVIVE THIS ANALYSIS AFTER CHARLES UPDATES WITH GENE ANNOTATIONS
View(svData %>% group_by(Type,IsGenicStart) %>% tally() %>% spread(IsGenicStart,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsFS","LengthBucket,IsGenic=IsGenicStart"),'LengthBucket','IsGenic','Length Distribution by isGenicStart',useLogY =F)

# INDIVIDUAL SIGNATURE CHECK
filter = cohortSummary(merge(svData,svDriverAndGermline %>% filter(gene=='CCNE1') %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(20,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),SampleId %in% filter","LengthBucket,IsGenic=IsGenicStart"),'LengthBucket','IsGenic','Length Distribution by isGenicEnd',useLogY =F)

#3. Zoom in on FragileSites
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),IsFS","LengthBucket,GeneStart"),'LengthBucket','GeneStart','Length Distribution by isGenicStart',useLogY =F)

#4.
View(svLinks %>% group_by(ResolvedType,isGenic=GeneStart!=''&GeneEnd!='') %>% count() %>%spread(isGenic,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))
##############################################################
############# INSERT and HOM Seq #############################
##############################################################

#1. INS SEQUENCE DISTRIBTUTION - 
plot_count_by_bucket_and_type(cohortSummary(svData,"nchar(InsertSeq)<50,nchar(InsertSeq)>=0,nchar(Homology)==0,!IsPolyA,Type!='SGL'","InsLen=nchar(InsertSeq),CustomType"),'InsLen','CustomType','InsertSequenceLength(excl PolyA)',F,T)



#2. Frequency of short distribution is relatively random - rate of As and Ts is proportional to ref genome rate.
View(svData %>% filter(!IsLowQual,ResolvedType!='Line',nchar(InsertSeq)<=3) %>% group_by((substring(InsertSeq,1,1))) %>% tally())

#3. Hom and INs length By resolved Type
# LINE, Chains & reciprocal events have less homology then simpleSV (~58% compared to ~68%)
# Chains and TIs have relative more INS. 
# Chains and LINE have more long INS in particulare (32+)
View(insertAndHomologySummary(svData,"!(Type %in% c('SGL','NONE','INS')),!IsPolyA,!IsLowQual",'ResolvedType'))

#4. Hom and INs length  for SimpleSV by Type and Length Bucket
### Short and mid-range DELs & DUPS are all around high 60%s with homology. Dels slightly lower than DUPs
## Very Short Dels have weird lack of INS SEQ
## Very Short DUPs have weird lack of Homology
View(insertAndHomologySummary(svData,"!(Type %in% c('SGL','NONE','INS')),!IsPolyA,!IsLowQual,ResolvedType=='SimpleSV'",'Type,LengthBucket'))

#5. POLYA in non SGL
#TO DO: why so many POLY A and T in unexpected cluster types
View(insertAndHomologySummary(svData,"!(Type %in% c('SGL','NONE','INS')),IsPolyA,!IsLowQual","ResolvedType") %>% arrange(-count))

#3. Source of long insert sequences
# TO DO: How much is local vs remote?  
View(svData %>% filter(Type!='SGL',Type!='NONE',!IsPolyA) %>% group_by(HasAln= InsSeqAlignments!='',nchar(InsertSeq)) %>% count() %>% spread(HasAln,n))

#4. VIRAL INSERSTIONS
# TO DO: Count frequency per viral genome cancer type. PORT from SQL

#5. SGL BE Insertion Sequence
# TO DO: Breakdown INS Sequence into known 
#• POLYA (presumed LINE)    ~40$
#• Viral    - 
#• Telomeric Sequence (towards telomere)  - 2%
#• Telomeric Sequence (away from telomere) - 1%
#• Centromeric Sequence (towards centromere) - 10%???
#• Centromeric Sequence (away from centromere) - 10%??
#• Other satellite DNA - ?
#• Other repeat sequences - ?
#• Unknown - ?

# eg. TELOMERE - 2.6k
View(svData %>% filter(Type=='SGL',(InsSeqAlignments)!='',!IsPolyA,(grepl('TTAGGGTTAGGG',InsertSeq)|grepl('AACCCTAACCCT',InsertSeq))) %>% group_by(OrientStart,fwdTelomereStrand=grepl('TTAGGGTTAGGG',InsertSeq)) %>% count())

# eg. CENTROMERE - at least 6k..
View(svData %>% filter(Type=='SGL',(InsSeqAlignments)!='',!IsPolyA,
                       (grepl('ACTTGAAACA',InsertSeq)|grepl('AACTTCCTG',InsertSeq)|grepl('CAGAGTTGAA',InsertSeq)|
                          grepl('GAAATATCTT',InsertSeq)|grepl('AGTGGAGATT',InsertSeq))) %>% count())

##############################################################
############# SUPSICIOUS Findings ############################
##############################################################

#4. TO DO: Investigate events in general with Qual score < 450
View(svData %>% mutate(lowQS=QualScore<450) %>% unite(lowQ,lowQS,IsLowQual) %>% group_by(SampleId,lowQ) %>% count() %>% spread(lowQ,n))

#5. Recurrent Variants
# TO DO: ensure we filter out duplicate samples and contamination
temp =(svData %>% filter(PosEnd != -1,SampleId %in% highestPurityCohortSummary$sampleId) %>% group_by(ChrStart,ChrEnd,PosStart,PosEnd,Type,Length=PosEnd-PosStart,OrientStart,OrientEnd) %>% 
        summarise(n=n(),min(Ploidy),max(Ploidy),mean(Ploidy),min(AdjCNChgStart),max(AdjCNChgStart),minSampleId=min(SampleId),maxSampleId=max(SampleId),mean(AdjCNChgStart),min(AdjCNChgEnd),max(AdjCNChgEnd),mean(AdjCNChgEnd)) %>% filter(n>1))
View(temp)
View(temp %>% group_by(maxSampleId,minSampleId) %>% count())
View(temp %>% filter(!(maxSampleId %in% c('CPCT02140045T','DRUP01140005T','DRUP01020026T','CPCT02020516T','DRUP01020027T','CPCT02020536T','CPCT02190031T','CPCT02190029T','DRUP01060010T','CPCT02370001T','DRUP01060008T','CPCT02060059T'))) %>% 
       mutate(PloidyRatio=`max(Ploidy)`/`min(Ploidy)`,sameSample=minSampleId==maxSampleId) %>% select(PloidyRatio,sameSample,everything()))

#6. Overlapping DB are associated with much shorter INV foldback events - these may be some other type of event????
print(ggplot(data = foldbacks %>% filter(FoldbackType=="INV",!IsLowQual) %>% group_by(FoldbackLenBucket,DBOverlap=((DBLenStart>-30&DBLenStart==0)|(DBLenEnd>-30&DBLenEnd<=0))) %>% summarise(Count=n()) %>% spread(DBOverlap,Count), 
             aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=`TRUE`, colour='Short Overlapping DB'))
      + geom_line(aes(y=`FALSE`, colour='Other' ))
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + scale_y_log10()
      + labs(title = "Foldback Length Distribution"))

#7. LONE INV
# 2/3 are Foldbacks and match the length disrtibution of foldbacks
View(svData %>% filter(ClusterCount==1,Type =='INV',!IsLowQual,ArmStart==ArmEnd) %>% arrange(SampleId,ChrStart,PosStart) %>%
       group_by(IsFoldBack,LengthBucket) %>% tally() %>% spread(IsFoldBack,n))
# surplus facing centromere
View(svData %>% filter(ClusterCount==1,Type =='INV',!IsLowQual,ArmStart==ArmEnd,IsFoldBack) %>% mutate(CentromereOriented=(OrientStart*ifelse(ArmStart=='P',-1,1))) %>%
       group_by(CentromereOriented,LengthBucket) %>% tally() %>% spread(CentromereOriented,n))
#select(ResolvedType,len,Ploidy,Type,ClusterCount,ClusterReason,everything()))

#8. LONE BND
#Around 3500 / 13000 lone BND appear to be geuine unrecpirocated translocations.  
# Rest likely missing a variant - should we try to join to SGL / NONE or INV
View(svData %>% filter(!IsLowQual,Type=='BND',ClusterCount==1) %>% group_by(CustomType,Consistency) %>% tally() %>% spread(Consistency,n))

#9. Should we still have imprecise?
# TO DO: check if recovered
View(svData %>% group_by(Imprecise,Type) %>% tally())

#11. Failure to cluster high ploidy INVs (not overlapping and not foldaback)
### FIND EXAMPLES


#############################################
########## Over Clustering ##################
#############################################

#3. Count of Foldbacks per cluster
# TO Do - Why do so many samples have 5+ separate clusters with foldbacks???
foldbackClusters=(foldbacks %>% filter(!IsLowQual) %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(FBcount=n()/2,
                                                                                                               ClusterCount=dplyr::first(ClusterCount),CNMax=max(pmax(AdjCNStart,AdjCNEnd)),FBPloidyMax=max(pmax(AdjCNChgStart,AdjCNChgEnd,Ploidy))))

ggplot(data=merge(sampleList,foldbackClusters %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of Foldback Clusters per Sample')
ggplot(data=foldbackClusters ,aes(FBcount)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of Foldback Count per Cluster') +scale_x_log10()
View(merge(sampleList,foldbackClusters %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>% arrange(-n))
View(foldbackClusters %>% group_by(SampleId,CC=pmin(5,ClusterCount)) %>% tally() %>% arrange(-n) %>% spread(CC,n)) 

View(foldbackClusters %>% group_by(FBcount,CC=pmin(5,ClusterCount)) %>% tally() %>% arrange(-n) %>% spread(CC,n)) 

View(foldbackClusters %>% group_by(ceiling(FBcount)) %>% tally())

#############################################
########## INDIVIDUAL SAMPLE ANALYSIS #######
#############################################

chartSample = 'CPCT02010690T'; chartChr = '11';
View(svData %>% filter(SampleId==chartSample,ResolvedType!='LowQual') %>% 
       select(IsFoldBack,ResolvedType,Length,ClusterId,Ploidy,Type,ClusterCount,ChrStart,ChrEnd,PosStart,PosEnd,ClusterReason,everything()) %>% arrange(-Ploidy))

#1. Whole Chr View
View(svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr,ResolvedType!='ALowQual') %>% 
       select(CancerType,IsFoldBack,ResolvedType,Length,SynDelDupLen,SynDelDupTILen,Ploidy,Type,ClusterCount,ClusterReason,everything()))
print(ggplot(data = svData %>% filter(!IsLowQual,SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>% mutate(
  modEnd=ifelse(ChrEnd==chartChr,ifelse(ChrStart==chartChr,PosEnd,as.numeric(ChrStart)*-1e7),as.numeric(ChrEnd)*-1e7),modStart=ifelse(ChrStart==chartChr,PosStart,PosEnd)),aes(modStart,modEnd)) + 
    geom_point(aes(size = Ploidy,colour=ifelse(ResolvedType == 'Line'|ClusterCount==2,ResolvedType,ClusterDesc))) + theme_bw() + scale_color_manual( values= myCOLORS))

#1a. all chromosomes
print(ggplot(data = rbind(svData %>% filter(!IsLowQual,SampleId==chartSample) %>% mutate(modStart=PosStart*1e7,modEnd=ifelse(ChrEnd==ChrStart,PosEnd,as.numeric(ChrEnd)*-1e7)) %>% select(ChrStart,Ploidy,ResolvedType,ClusterCount,ClusterDesc,modStart,modEnd),
                          svData %>% filter(!IsLowQual,ChrEnd!=ChrStart,SampleId==chartSample,ChrEnd!=0) %>% mutate(modStart=PosEnd*1e7,modEnd=as.numeric(ChrStart)*-1e7) %>% select(ChrStart,Ploidy,ResolvedType,ClusterCount,ClusterDesc,modStart,modEnd)),
             aes(modStart,modEnd)) + geom_point(aes(size = Ploidy,colour=ifelse(ResolvedType == 'Line'|ClusterCount==2,ResolvedType,ClusterDesc))) + theme_bw() +facet_wrap(~ChrStart))

#2.By CHR for large clusters
View(svData %>% filter(SampleId == chartSample,ClusterCount>8) %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% tally() %>% spread(chr_arm_end,n,fill=""))

View(svData %>% filter(SampleId == chartSample,IsFoldBack) %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% tally() %>% spread(chr_arm_end,n,fill=""))

# 3. Single Cluster View
View(svData %>% filter(SampleId==chartSample,ClusterId==59)%>% arrange(ChrStart,PosStart) 
     %>% select(ResolvedType,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,Type,ClusterCount,ClusterReason,everything()))

# 4. Large cluster View
View(svData %>% filter(!IsLowQual,SampleId==chartSample,ClusterCount>=1,ClusterCount>0,ResolvedType!='ASimpleSV')%>% 
       select(ResolvedType,Length,Ploidy,Type,ClusterCount,ClusterReason,everything()) %>% arrange(-Length))#%>% group_by(ClusterCountBucket,CnChStartBucket) %>% tally() %>% spread(CnChStartBucket,n))










############################################### EXPERIMENTAL BELOW HERE ##############################################



##################################################
############### SV CLUSTERS ######################
##################################################

View(svClusters %>% filter(ClusterCount==3) %>% group_by(AcComplexFoldback,AcSimpleFoldback,AcComplexOther,AcMultipleDsb,AcDsb,AcRemoteTI,AcSoloSv) %>% tally())
scatterPlot(svClusters %>% filter(ClusterCount>5),'AcComplexFoldback','AcSimpleFoldback')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcDsb')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcRemoteTI')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcSoloSv')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','ShortTIRemotes')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','UnlinkedRemotes')

View(svClusters %>% filter(ClusterCount>0,ResolvedType!='Line',AcSimpleFoldback>0) %>% group_by(SampleId) %>% tally())

##################################################
############### CONSECUTRIVE BRREAKENDS ##########
##################################################

consecutiveBEs = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/consecutive_be.csv')
consecutiveBEs$HasBEBetween = (consecutiveBEs$TraversedSvCount>0)
consecutiveBEs$NextNonClusteredSVOrientMatch = (consecutiveBEs$Orientation!=consecutiveBEs$OrientNext)

# 1. All 81k
nrow(consecutiveBEs)

# 2. INCONSISTENT CN Change / Likely FP - 19k
nrow(consecutiveBEs %>% filter(CNChgFront<=0.5|CNChgBack<=0.5))

# 3. FB - 13k
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,IsFoldback=='true'))

View(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,IsFoldback=='true',FtoNextCLLength!=0) %>% select(CNChgBack,CNChgFront,CNChgNextCL,everything()))
### 3.A) FB with no next Variant - 2k
View(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,IsFoldback=='true',FtoNextCLLength==0,NextNonClusteredSVOrientMatch) %>% select(CNChgBack,CNChgFront,CNChgNext,everything()))
ggplot(data=consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,IsFoldback=='true',FtoNextCLLength==0,NextNonClusteredSVOrientMatch),aes(FtoNextLength)) + 
  stat_ecdf(geom = "step", pad = FALSE) + scale_x_log10() + ylim(0,1) + labs(title = 'CDF of next unclustered SV in foldback')

# Has BE - 9k
View(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==T) %>% select(CNChgBack,CNChgFront,CNChgNextCL,everything()))

# DB - 2k
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextCLLength<30,FtoNextCLLength>0))

# No next clustered variant 2k
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextCLLength==0,FtoNextLength>0))

# No next Variant at all 1k
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextLength==0))

# Possible DUP BE - 2k
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextCLLength>=30,BtoFLength<30))

# No CN Change
nrow(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextCLLength>=30,BtoFLength>=30,abs(CNChgBack-CNChgNextCL)<0.5|abs(CNChgNextCL-CNChgBack)/pmax(CNChgNextCL,CNChgBack) < 0.1,
                               (CNChgBack-CNChgFront) < -0.5,(CNChgNextCL-CNChgFront) < -0.5))

#Other
temp = (consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextLength>=5000,BtoFLength>=30,TypeBack!='SGL',TypeFront!='SGL') %>% group_by(SampleId,ClusterId,ClusterCount,Long=BtoFLength>5000) %>% summarise(n=n()) %>% spread(Long,n))
scatterPlot(temp,'ClusterCount','n',T,T)
View(temp)

temp=(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,HasBEBetween==F,IsFoldback=='false',FtoNextCLLength>=30,BtoFLength>=3) %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(n=n(),maxCN=max(CNChgFront)) )
scatterPlot(temp,'maxCN','n',T,T)
View(consecutiveBEs %>% filter(CNChgFront>0.5,CNChgBack>0.5,FtoNextCLLength>=30,SampleId=='CPCT02010869T',ClusterId==233))


############# TESTING #############

View(svData %>% filter(grepl('ATATATAT',InsertSeq)) %>% group_by(Type,ResolvedType,PloidyBucket) %>% count() %>% spread(Type,n))


View(svDrivers %>% filter(Gene=='PTCH1',CancerType=='Skin') %>% group_by(SampleId) %>% tally())

View(svData %>% filter(SampleId =='DRUP01010099T',ClusterId==77))

View(svData %>% group_by(ResolvedType) %>% tally())
# 1. INSERTION MOTIF: ~15% are A-TTTTT with another 30% or so a similar variant
View(dbData %>% filter(RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context,ResolvedType) %>% tally() %>% arrange(-n) %>% spread(ResolvedType,n))


repTiming = read.csv(paste(PATH,'heli_rep_origins.bed',sep=''), header = F, stringsAsFactors = F,sep = '\t',col.names = c('chr','posStart','posEnd','repTiming'))
(repTiming %>% group_by(round(repTiming,-1)) %>% count())
ggplot(data=repTiming %>% filter(chr=='chrX'),aes(as.numeric(posStart),as.numeric(repTiming)))+geom_point() + xlim(2e7,4e7)
head(repTiming)hea
     
View(svData %>% filter(ChrEnd!=0,ResolvedType=='SimpleSV',RepOriginStart!=0,RepOriginEnd!=0) %>% 
       group_by(Type,Length>5000,RepStart=round(RepOriginStart,1),RepEnd=round(RepOriginEnd,1)) %>% tally() %>% spread(RepEnd,n))

View(svData %>% filter(ChrEnd!=0,ResolvedType=='SimpleSV',RepOriginStart!=0,RepOriginEnd!=0) %>% 
       group_by(Type,ChrStart,RepEnd=round(RepOriginEnd,1)) %>% tally() %>% spread(RepEnd,n))

View(svData %>% filter(ChrEnd!=0,RepOriginStart!=0,RepOriginEnd!=0) %>% 
       group_by(ResolvedType,ClusterCountBucket,RepEnd=round(RepOriginEnd,1)) %>% tally() %>% spread(RepEnd,n))

plot_count_by_bucket_and_type(cohortSummary(dbData %>% mutate(CustomType=ifelse(CustomType=='Chain'&Length<1e3&Type=='INV','ShortINVinChain',CustomType)),
                                            "DBLength<=30,!IsLowQual,!(CustomType %in% c('SimpleSV','PairedSGL')),RepOriginStart<0.6|RepOriginEnd<0.6",'DBLength,CustomType'),
                              'DBLength','CustomType','DBLengthByCustomType(<50bases)',useLogX = F,useLogY = F)

View(svData %>% filter(Type=='DUP',ResolvedType=='SimpleSV',Length<100) %>% 
       group_by(PosStart,PosEnd,ChrStart,ChrEnd) %>% mutate(count=n()) %>% ungroup() %>% group_by(SampleId,count) %>% count() %>% spread(count,n))

