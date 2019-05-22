library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)

base_complements <- function(bases) {
  complements = setNames(c("A", "C", "G","T"), c("T", "G", "C","A"))
  
  point_complement <- function(base) {
    paste(rev(sapply(strsplit(base, split = ""), function (x) {complements[x]})), collapse = "")
  }
  
  sapply(bases, point_complement)
}

#' Modified version of dplyr's filter that uses string arguments
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df
}

# plotting tools
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
  plot <- ggplot(data=countsData,aes_string(x=bucket))+
      geom_line(aes(y=countDEL,colour='DEL'))+
      geom_line(aes(y=countDUP,colour='DUP'))+
      geom_line(aes(y=countINV,colour='INV'))+
      geom_line(aes(y=countBND,colour='BND'))+
      geom_line(aes(y=countSGL, colour='SGL'))+
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
  IsFoldBack = FoldbackLenStart>=0|FoldbackLenEnd>=0,
  IsPolyA = grepl('TTTTTTTTTT',InsertSeq)|grepl('AAAAAAAAAA',InsertSeq)
)
}

createBuckets <- function(cluster) {
  cluster %>% mutate(
    PloidyBucket=2**(pmin(7,pmax(-3,round(log(Ploidy,2),0)))),
    CnChEndBucket=2**(pmin(7,pmax(-3,round(log(pmax(CNChgEnd,0.01),2),0)))),
    CnChStartBucket=2**(pmin(7,pmax(-3,round(log(pmax(CNChgStart,0.01),2),0)))),
    ClusterCountBucket=2**(pmin(5,pmax(-3,round(log(ClusterCount,2),0)))),
    LnkLenStartBucket=ifelse(!LnkLenStart>0,0,2**(pmin(20,pmax(0,round(log(pmax(LnkLenStart,0.01),2),0))))),
    FoldBackLenStartBucket=ifelse(!FoldbackLenStart>0,0,2**(pmin(20,pmax(0,round(log(pmax(FoldbackLenStart,0.01),2),0))))),
    LnkLenEndBucket=ifelse(!LnkLenEnd>0,0,2**(pmin(20,pmax(0,round(log(pmax(LnkLenEnd,0.01),2),0))))),
    FoldBackLenEndBucket=ifelse(!FoldbackLenEnd>0,0,2**(pmin(20,pmax(0,round(log(pmax(FoldbackLenEnd,0.01),2),0))))),
    LengthBucket=ifelse(Type=='BND'|Type=='INS'|PosEnd-PosStart==0|ArmEnd!=ArmStart,0,2**(pmin(25,pmax(5,round(log(pmax(0.01,PosEnd-PosStart),2),0))))),
    HomLenBucket=2**(round(log(nchar(Homology),2),0)),
    NearestLenBucket=2**(round(log(pmax(0.01,NearestLen),2),0)),
    SynDelDupLenBucket=2**(round(log((SynDelDupLen),2),0)),
    SynDelDupTILenBucket=2**(round(log((SynDelDupAvgTILen),2),0)),
    InsLenBucket=2**(round(log(nchar(InsertSeq),2),0)),
    QualScoreBucket=2**(round(log(nchar(QualScore),2),0))
  )
}

createsvLinksBuckets <- function(svLinks) {
  svLinks %>% mutate(
    TILengthBucket=ifelse(TILength==0,0,2**(pmin(25,pmax(5,round(log(TILength,2),0))))),
    synDelDupLengthBucket=ifelse(SynDelDupLen==0,0,2**(pmin(25,pmax(5,round(log(pmax(0.01,SynDelDupLen),2),0))))),
    DBLenEndBucket=ifelse(DBLenEnd==0,0,2**(pmin(25,pmax(5,round(log(pmax(0.01,DBLenEnd),2),0))))),
    DBLenStartBucket=ifelse(DBLenStart==0,0,2**(pmin(25,pmax(5,round(log(pmax(0.01,DBLenStart),2),0))))),
    ClusterCountBucket=2**(pmin(5,pmax(-3,round(log(ClusterCount,2),0))))
  )
}

createFoldbacks <- function(cluster) {
  rbind(cluster %>% filter(FoldbackLenStart>=0) %>% mutate(FoldbackLength=FoldbackLenStart,FoldbackLinkInfo=FoldbackInfoStart),
        cluster %>% filter(FoldbackLenEnd>=0,!FoldbackLenStart>=0)  %>% mutate(FoldbackLength=FoldbackLenEnd,FoldbackLinkInfo=FoldbackInfoEnd)) %>% 
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
              countDUP_EXT_TIShort=sum(ResolvedType=='DUP_EXT_TI'&SynDelDupLen<1e5),
              countDUP_EXT_TILong=sum(ResolvedType=='DUP_EXT_TI'&SynDelDupLen>=1e5)) %>% group_by(SampleId) %>% 
    summarise(countCluster=n(),
              count=sum(n),
              countLine=sum(ResolvedType=='LINE'),
              sumLine=sum(ifelse(ResolvedType=='LINE',n,0)),
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
              countDEL_EXT_TI=sum(ResolvedType=='DEL_EXT_TI'),
              countDEL_INT_TI=sum(ResolvedType=='DEL_INT_TI'),
              countDUP_EXT_TI=sum(ResolvedType=='DUP_EXT_TI'),
              countDUP_EXT_TIShort=sum(countDUP_EXT_TIShort),
              countDUP_EXT_TILong=sum(countDUP_EXT_TILong),
              countDUP_INT_TI=sum(ResolvedType=='DUP_INT_TI'),
              countRecipInv=sum(ResolvedType=='ReciprocalInversion'),
              countRecipTrans=sum(ResolvedType=='RecipTrans'),
              #countShortTI=sum(ShortTICount),
              countFS=sum(countFS))
}

#############################################
################ LOADING ####################
#############################################

PATH='~/hmf/analyses/SVAnalysis/'

svData = read.csv(paste(PATH,'SVA_SVS.csv',sep=''), header = T, stringsAsFactors = F)
svClusters = (read.csv(paste(PATH,'SVA_CLUSTERS.csv',sep=''), header = T, stringsAsFactors = F))
svDrivers = read.csv(paste(PATH,'SVA_DRIVERS.csv',sep=''), header = T, stringsAsFactors = F)
svLinks = read.csv(paste(PATH,'SVA_LINKS.csv',sep=''), header = T, stringsAsFactors = F)
sampleCancerTypes= (read.csv(paste(PATH,'sample_cancer_types.csv',sep=''), header = T, stringsAsFactors = F))
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
svData=merge(svData,svClusters %>% select(SampleId,ClusterId,Subclonal,SynDelDupLen,SynDelDupAvgTILen,Foldbacks),by=c('SampleId','ClusterId'))
svDrivers = merge(svDrivers, sampleCancerTypes, by='SampleId', all.x=T)
sampleList = svData %>% distinct(SampleId,CancerType)

#svChordStatus = read.csv(paste(PATH,'SVChordStatus.csv',sep=''), header = T, stringsAsFactors = F)
#svDriverCatalog = read.csv(paste(PATH,'SVDriverCatalog.csv',sep=''), header = T, stringsAsFactors = F)
#svGermline = read.csv(paste(PATH,'SVGermline.csv',sep=''), header = T, stringsAsFactors = F)
#svDriverAndGermline=(rbind(svDriverCatalog %>% select(gene,sampleId,driver,driverLikelihood), svGermline %>% mutate(driverLikelihood=1,driver=ifelse(biallelic==1,'Germline Biallleic','Germline Monoallleic')) %>% select(gene,sampleId,driverLikelihood,driver)))

#Annotations
svData=createBuckets(svData)
svData=sv_set_common_fields(svData)
svLinks=createsvLinksBuckets(svLinks)
foldbacks=createFoldbacks(svData)
svSampleSummary=createSampleSummary(svData)

#RI temp fix
#svLinks = svLinks %>% mutate(ResolvedType=ifelse(ResolvedType=='DEL_INT_TI'&SynDelDupAvgTILen>SynDelDupLen*0.5,'RECIP_INV',ResolvedType))
#svLinks = svLinks %>% mutate(ResolvedType=ifelse(ResolvedType=='DEL_INT_TI'&TILength>SynDelDupLen*0.5,'RECIP_INV',ResolvedType))

# BE data
beData = rbind(svData %>% mutate(LocTopTI=LocTopTIStart,LocTopType=LocTopTypeStart,chr=ChrStart,pos=PosStart,Orient=OrientStart,RefContext=RefContextStart,LE=LEStart,DBLength = DBLenStart,Assembled = ifelse(AsmbMatchStart=="MATCH","Assembled","NotAssembled")),
               svData %>% mutate(LocTopTI=LocTopTIEnd,LocTopType=LocTopTypeEnd,chr=ChrEnd,pos=PosEnd,Orient=OrientEnd,RefContext=RefContextEnd,LE=LEEnd,DBLength = DBLenEnd, Assembled = ifelse(AsmbMatchEnd=="MATCH","Assembled","NotAssembled"))) 
dbData = beData %>% filter(DBLength>-1000) %>% mutate(DBLenBucket = ifelse(DBLength==0,0,ifelse(DBLength<0,-(2**round(log(-DBLength,2))),2**round(log(DBLength,2)))))

#############################################
########## OVERVIEW OF COHORT ###############
#############################################

# 1. By RESOLVED Type
View(svData %>% filter(Subclonal=='false') %>% group_by(ResolvedType,Type) %>% tally() %>% spread(Type,n))

# 2. By CLUSTERCOUNTBUCKET
View(svData %>% filter(Subclonal=='false',ResolvedType!="Line") %>% group_by(ClusterCountBucket,Type) %>% tally() %>% spread(Type,n))

#3. Simple SV make up 70% of clusters but only 30% of SVs
View(svClusters %>% filter(Subclonal=='false') %>% group_by(ResolvedType) %>% summarise(variants=sum(ClusterCount),clusters=n()))

############################################################
########### Short TIs #####################################
############################################################
#0. Numbers per sample
ggplot(data=merge(svData %>% group_by(SampleId) %>% tally(),svLinks %>% 
                    group_by(SampleId,shortTI=ifelse(TILength<1000,'ShortTICount','LongTICount')) %>% tally() %>% spread(shortTI,n),by='SampleId')) +
  geom_point(aes(n,ShortTICount))

#1. TIs length distibution, showing 2 peaks for most types ('short', ie <1k and long -> presumed to be cause by random length of distance between breakages)
ggplot(data=svLinks %>% filter(ResolvedType!='COMPLEX') %>% group_by(TILengthBucket,ResolvedType) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_x_log10()+ facet_wrap(~ResolvedType)

#2. LONG LINE CLUSTER LENGTHS Issue
ggplot(data=svLinks %>% filter(ResolvedType=='LINE') %>% group_by(TILengthBucket,ClusterCountBucket) %>% summarise(count=n()),
aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_x_log10()+ facet_wrap(~ClusterCountBucket)

#3. SHORT TI actual shape
ggplot(data=svLinks  %>% filter(TILength<5000) %>% mutate(TILengthBucket=round(TILength+0.001,-1)) %>% group_by(TILengthBucket) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_y_log10()

#4. Short TIs create 'synthetic' events that match the simpleSV event profile per sample
scatterPlot(svSampleSummary,'countSimpleDEL','countDEL_EXT_TI',F,F) # DEL correlated with synthetic DEL
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_EXT_TI',F,F)  # DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_INT_TI',F,F)  # Int DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPShort','countDUP_EXT_TIShort',F,F)  # Short DUP correlated with short synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPShort','countDUP_EXT_TILong',F,F)  # Short DUP correlated with short synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPLong','countDUP_EXT_TILong',F,F)  # Long DUP correlated with long synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPLong','countDUP_EXT_TIShort',F,F)  # Short DUP correlated with short synthetic DUP

#5. Synthetic Length Distributions
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='DEL_INT_TI'",'SynDelDupLenBucket,SynDelDupTILenBucket'),'SynDelDupTILenBucket','SynDelDupLenBucket','DEL_INT_TI',useLogY =F)
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='DEL_EXT_TI'",'SynDelDupLenBucket,SynDelDupTILenBucket'),'SynDelDupTILenBucket','SynDelDupLenBucket','DEL_EXT_TI',useLogY =F)
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='DUP_EXT_TI'",'SynDelDupLenBucket,SynDelDupTILenBucket'),'SynDelDupTILenBucket','SynDelDupLenBucket','DUP_EXT_TI',useLogY =F)
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='DUP_INT_TI'",'SynDelDupLenBucket,SynDelDupTILenBucket'),'SynDelDupTILenBucket','SynDelDupLenBucket','DUP_INT_TI',useLogY =F)


##############################################################
########## Short Deletion Bridge Analysis ####################
#############################################################

##### TEMP

#### WHY???
plot_count_by_bucket_and_type(cohortSummary(dbData, "DBLength<=100,DBLength>=-100,Subclonal=='false',ResolvedType=='COMPLEX',LocTopType=='DSB',LocTopTI<3",'DBLength,LocTopTI'),
                              'DBLength','LocTopTI','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = F)

View(svData %>% filter(DBLenStart==0,DBLenEnd==0,Subclonal=='false',ResolvedType=='COMPLEX',LocTopTypeStart=='DSB',LocTopTypeEnd=='DSB',ClusterDesc =='INV=4')
     %>%select(Length,everything()))#


View(svData %>% filter(DBLenStart==0,DBLenEnd==0,Subclonal=='false',ResolvedType=='COMPLEX')
     %>%select(Length,everything()))# %>% group_by(ClusterDesc,Type) %>% count() %>% spread(Type,n)) 
plot_count_by_bucket_and_type(cohortSummary(dbData, "DBLength<=30,DBLength>=-30,Subclonal=='false',ResolvedType=='COMPLEX',Type!='BND'",'DBLength,LengthBucket'),
                              'DBLength','LengthBucket','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = F)

plot_count_by_bucket_and_type(cohortSummary(dbData, "DBLength<=30,DBLength>=-30,Subclonal=='false',ResolvedType=='COMPLEX',!IsFoldBack",'DBLength,ClusterCountBucket'),
                              'DBLength','ClusterCountBucket','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = F)



#1.DBLength by ResolvedType => less than 50 bases (NB - DBLength of 1 means exact break - shoould correct this)
# Note the LINE double peak and the sharp feature for Reciprocal Inversion, Reciprocal Translocations and None
# NB LE=='None' excludes line source elements with complex topology
plot_count_by_bucket_and_type(cohortSummary(dbData, "DBLength<=100,DBLength>=-100,Subclonal=='false',LE=='None'",'DBLength,ResolvedType'),
                              'DBLength','ResolvedType','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = T)

#2. EXACT BREAKS
# RI have a sharp feature at 0 bases reagrdless of RI length!
plot_count_by_bucket_and_type(cohortSummary(dbData,"Type!='BND',DBLength>-50,DBLength<=50,Subclonal=='false',ResolvedType=='RECIP_INV'",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for INV in RI by length bucket',useLogX = F,useLogY = F)
# Sharp RI feature exists for Complex short INV also. Excess compared to DEL & DUP and Mostly around 60-500 bases (and not a foldback feature)
plot_count_by_bucket_and_type(cohortSummary(dbData,"Type!='BND',DBLength>-50,DBLength<=50,Subclonal=='false',Length<1e7,ResolvedType=='COMPLEX'",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for short variants in COMPLEX clusters by length bucket',useLogX = F,useLogY = F)

plot_count_by_bucket_and_type(cohortSummary(dbData,"Type!='BND',DBLength>-50,DBLength<=50,Subclonal=='false',Length<1e7,ResolvedType=='COMPLEX'",'DBLength,LocTopType'),'DBLength','LocTopType','DBLength for short variants in COMPLEX clusters by length bucket',useLogX = F,useLogY = F)


#3 DUP EXT_TI which are likely LINE
plot_count_by_bucket_and_type(cohortSummary(dbData %>% filter(ResolvedType=='DUP_EXT_TI'),
                                            "SynDelDupLen<=30,Subclonal=='false'",'SynDelDupLen,IsPolyA'),
                              'SynDelDupLen','IsPolyA','SyntheticDelDupLength',useLogX = F,useLogY = F)

#4. SOME SINGLES in COMPLEX 2 clusters and above LOOK LIKE LINE
plot_count_by_bucket_and_type(cohortSummary(dbData, "DBLength<=100,DBLength>=-100,Subclonal=='false',LocTopType=='DSB',LE=='None',ClusterCount==2,ResolvedType=='COMPLEX'",'DBLength,ClusterDesc'),
                              'DBLength','ClusterDesc','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = F)


#5. LINE Elements - no clear indication of why there are 2 DB peaks.
#The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific.   Also appars to be the same regardless of source element
print(ggplot(data = dbData %>% filter(DBLength<=50,ResolvedType %in% c('LINE')) %>% group_by(CancerType,OLPeak=DBLength<(-7)) %>% tally(), aes(x = reorder(CancerType, -n), y = n, fill =OLPeak))
      + geom_bar(stat = "identity", colour = "black") + ylab("DB Count") + xlab("Tumor Type") + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15)) + labs(title = "DB Counts for LINE Elements by cancerType and Overlapping Peak"))
# Also does not seem to depend significantly on refContext.  
View(dbData %>% filter(ResolvedType=='LINE',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context,longOverlap=DBLength<(-7)) %>% tally() %>% arrange(-n) %>% spread(longOverlap,n))


################################################
############## FOLDBACKS #######################
################################################

foldbacks=createFoldbacks(svData)
View(svClusters %>% filter(Foldbacks>0) %>% group_by(Subclonal,Foldbacks) %>% tally %>% spread(Subclonal,n))

#1.Foldback length distribution for simple inversions short chained foldbacks
print(ggplot(data = foldbacks %>% filter(Subclonal=='false') %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
             aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=INV, colour='INV')) + geom_line(aes(y=Combo, colour='Chain' )) 
      + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))


# TO DO: WHAT DO WE WANT TO DO WITH OUR DETAILED FB analysis?
# B. When we have 2 FB how often are they facing each other?
# C. How often is the most telomeric FB facing to the centromere

#2. Per sample counts (needs work)
temp=(svData %>% filter(Subclonal=='false',IsFoldBack,ClusterCount>0) %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(FBCount=sum(ifelse(IsFoldBack,1,0)),maxPloidy=max(Ploidy)))
View(temp %>% group_by(SampleId) %>% count)
scatterPlot(temp,'FBCount','maxPloidy',F,T)

# LOCATION
telomereOrientedFB=(svData %>% filter(CancerType!='ABreast',Subclonal=='false',IsFoldBack,Type=='INV',PosEnd-PosStart<1e4) %>% mutate(nrows=n(),FBOrient=ifelse(OrientStart*ifelse(ArmStart=='P',1,-1)==1,'TELOMERE','CENTROMERE')) %>% filter(FBOrient=="TELOMERE") %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% summarise(n=n(),avgRep=mean(RepOriginStart)) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
centromereOrientedFB=(svData %>% filter(CancerType!='ABreast',Subclonal=='false',IsFoldBack,Type=='INV',PosEnd-PosStart<1e4) %>% mutate(nrows=n(),FBOrient=ifelse(OrientStart*ifelse(ArmStart=='P',1,-1)==1,'TELOMERE','CENTROMERE')) %>% filter(FBOrient=="CENTROMERE") %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% summarise(n=n(),avgRep=mean(RepOriginStart)) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot() + geom_point(data=telomereOrientedFB,aes(pos,n,color="TELOMERE"),shape="x") + geom_point(data=centromereOrientedFB,aes(pos,n,color="CENTROMERE"),shape="x") + facet_wrap(~chr)

#######################################################
#################### LINE #############################
#######################################################

# 1. INSERTION MOTIF: ~15% are A-TTTTT with another 30% or so a similar variant
View(dbData %>% filter(Subclonal=='false',ResolvedType=='LINE',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context) %>% tally() %>% arrange(-n))

# 2. LINE Clusters affect a predictable number of chromosome arms that grows with the svCount in the cluster
scatterPlot(svClusters %>% filter(Subclonal=='false',ResolvedType=='LINE'),'ClusterCount','ArmCount',F,F)

# 3. LINE elements mmuch more common in Esophagus and Stomach cancers;  
#A. Source Elements
ggplot(merge(sampleList,svClusters %>% filter(Subclonal=='false',ResolvedType=='LINE',!(ClusterDesc %in% c('SGL','SGL=2','INS'))) %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>%
         group_by(CancerType) %>% mutate(TypeCount=n()) %>% filter(TypeCount>20) %>% ungroup() , aes(CancerType, n)) + 
  geom_violin(scale="width",fill="#6baed6") +   theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) + ggtitle("") + xlab("") + ylab("Line Source Elements Per Sample") + 
  theme(axis.ticks = element_blank(), legend.position="bottom") + scale_y_continuous(expand=c(0.01, 0.01)) + scale_fill_manual(values = cancerTypeColours) + theme(legend.position="none") +coord_flip() 
#B. Total LINE related SVs
ggplot(merge(sampleList,svData%>% filter(Subclonal=='false',ResolvedType=='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>%
         group_by(CancerType) %>% mutate(TypeCount=n()) %>% filter(TypeCount>20) %>% ungroup() , aes(CancerType, n)) + 
  geom_violin(scale="width",fill="#6baed6") +   theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) + ggtitle("") + xlab("") + ylab("Line Total SV Per Sample") + 
  theme(axis.ticks = element_blank(), legend.position="bottom") + scale_y_continuous(expand=c(0.01, 0.01)) + scale_fill_manual(values = cancerTypeColours) + theme(legend.position="none") +coord_flip() 

# 5. Known and suspected LINE element source locations => many potentially novel source locations that feature
View(beData %>% filter(Subclonal=='false',LE!='None',ResolvedType=='LINE') %>% group_by(chr,pos=round(pos,-4),SampleId) %>% 
     summarise(n=n(),countSuspected=sum(LE=='Suspect'),countKnown=sum(LE=='Known'),countPOLYA=sum(IsPolyA&Type=='BND'),countSuspectKnown=sum(LE=='Known;Suspect')) %>%
     group_by(chr,pos) %>% summarise(n=sum(n),countSuspected=sum(countSuspected),countSuspectKnown=sum(countSuspectKnown),countKnown=sum(countKnown),countPOLYA=sum(countPOLYA),numSamples=n()) )

LineClusters=(beData %>% filter(Subclonal=='false',ResolvedType=='LINE') %>% group_by(SampleId,ClusterId) %>% 
       summarise(n=n()/2,countAllSuspected=sum(LE=='Suspect'|LE=='Known;Suspect'),countKnown=sum(LE=='Known'),countPOLYA=sum(IsPolyA)/2,
                 countSuspectKnown=sum(LE=='Known;Suspect'),
                 countSGL=sum(Type=='SGL')/2,countINV=sum(Type=='INV')/2,countBND=sum(Type=='BND')/2,countDEL=sum(Type=='DEL')/2,countDUP=sum(Type=='DUP')/2))

View(LineClusters %>% filter(n==4))
scatterPlot(LineClusters%>% filter(countAllSuspected==0),'countPOLYA','countKnown',F,F)
scatterPlot(LineClusters %>% mutate(n=pmin(n,100),countAllSuspected=pmin(100,countAllSuspected+countKnown)),'n','countAllSuspected',F,F)
scatterPlot(LineClusters,'countAllSuspected','countPOLYA',F,F)
#######################################################
########## Pseudogene insertions ###################
#######################################################

svPseudoGeneLinks = svLinks %>% filter(ExonMatch!='') %>% separate(ExonMatch,c('transcriptId','exonRank','exonLength','other')) #what is the last field in the ExonMatch?

#1. 440 events
View(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd,transcriptId) %>% summarise(count=n(),minExon=min(as.numeric(exonRank)),maxExon=max(as.numeric(exonRank)),range=maxExon-minExon+1))

#2. found in 6% of samples - nearly all with highly activated LINE elements
View(merge(merge(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd) %>% count() %>% group_by(SampleId) %>% count(),sampleList,by='SampleId',all=T),
           svClusters %>% filter(ResolvedType=='LINE') %>% group_by(SampleId) %>% summarise(numLine=n()),by='SampleId',all.x = T) %>% arrange(-nn))


##TO DO: Create circos charts with gene transcripts for further analysis

############################################################
########### Reciprocal Inversions ##########################
############################################################

# #1. Reciprocal inv
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType %in% c('DEL_INT_TI','RECIP_INV'),Subclonal=='false'",'ResolvedType,LengthBucket') 
                              ,'LengthBucket','ResolvedType','Reciprical INV',useLogY =F )

#2. Reciprocal INV NOT enriched in typical DEL lengths OR FS!
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='SIMPLE'|ResolvedType=='RECIP_INV'",'IsFS,LengthBucket'),'LengthBucket','IsFS','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )


############################################################
########### Reciprocal Translocations ######################
############################################################

# TO DO:  show the length distribution and come up with a standardised definition
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='RECIP_TRANS'",'CancerType,NearestLenBucket') 
                              ,'NearestLenBucket','CancerType','Reciprical Trans',useLogY =F )

# TO DO: How to find synthetic RTs?

##############################################################
########## FRAGILE SITE ######################################
##############################################################

#1. Fragile Site Length Distribution
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType %in% c('SIMPLE','RECIP_INV')","LengthBucket,IsFS"),'LengthBucket','IsFS','Length Distribution of SimpleSV & Reciprocal Inversion by IsFS',useLogY =F)

#2. FS like length dels outside FS are correlated, but not all 20k-500k are FS Dels eg. CPCT02300036T and CPCT02190030T
# TODO: COLOUR this plot by cancer type
scatterPlot(svData %>% filter(Type=='DEL'&ResolvedType=='SIMPLE') %>% group_by(SampleId) %>% 
              summarise(countFS=sum(IsFS),countFSLike=sum(!IsFS&PosEnd-PosStart>2e4&PosEnd-PosStart<5e5)),'countFS','countFSLike',F,F) 

#3. Extreme tail possible for FS length DELS => Unrelated process
ggplot(data=svData %>% filter(ResolvedType=='SIMPLE',PosEnd-PosStart>20000,PosEnd-PosStart<500000,Type=='DEL') %>% group_by(CancerType,SampleId,IsFS) %>% tally() %>% spread(IsFS,n,fill=0)) + 
  stat_ecdf(aes(`TRUE`,color='isFS'),geom = "step", pad = FALSE) + stat_ecdf(aes(`FALSE`,color='FSLength non-FS'),geom = "step", pad = FALSE) + 
  ylim(0,1) + labs(title = 'CDF # of FS length DELs') + facet_wrap(~CancerType)

#4. Types of events in FS -  TO DO: Clean up and divided Simple events by length category
View(svData %>% group_by(IsFS,ResolvedType,Type) %>% tally() %>% spread(IsFS,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))

############################################################
############# HOMOLOGY & INSERT SEQUENCE ###################
############################################################

#1. INS SEQUENCE DISTRIBTUTION - 
View(svData  %>% filter(Type!='NONE',Type=='BND'|Length>100,Subclonal=='false',!IsPolyA) %>%  group_by(ResolvedType,InsLenBucket) %>% count() %>% spread(InsLenBucket,n))
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',nchar(InsertSeq)<50,nchar(InsertSeq)>=0,!IsPolyA,nchar(Homology)==0,Type!='SGL',Type!='NONE'","InsLen=nchar(InsertSeq),ResolvedType"),'InsLen','ResolvedType','InsertSequenceLength(excl PolyA)',F,T)

#2. HOM SEQUENCE DISTRIBTUTION (NO homology more common in complex & reciprocal breaks)
View(svData  %>% filter(Type!='NONE',InsertSeq=="",Type=='BND'|Length>100,Subclonal=='false') %>%  group_by(ResolvedType,HomLenBucket) %>% count() %>% spread(HomLenBucket,n))
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',nchar(InsertSeq)==0,Type!='SGL',Type!='NONE'","HomLen=nchar(Homology),ResolvedType"),'HomLen','ResolvedType','HomSeq',F,T)

### TO DO: Short and mid-range DELs & DUPS are all around high 60%s with homology. Dels slightly lower than DUPs.   Check previous short variant resultes

#3. Frequency of short distribution is relatively random - rate of As and Ts is proportional to ref genome rate.
View(svData %>% filter(Subclonal=='false',nchar(InsertSeq)>0,Type!='SGL',!IsPolyA,Type!='NONE',nchar(InsertSeq)<=3) %>% group_by((substring(InsertSeq,1,1))) %>% tally())

#### TO DO: SGL BE Insertion Sequence - DANIEL











#### TO CLEAN UP########

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


########################################
########## GENIC #######################
########################################

# TO DO: REVIVE THIS ANALYSIS AFTER CHARLES UPDATES WITH GENE ANNOTATIONS
View(svData %>% group_by(Type,IsGenicStart) %>% tally() %>% spread(IsGenicStart,n))# %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SIMPLE')|(ResolvedType=='DEL_INT_TI'&Type=='INV'),!IsFS","LengthBucket,IsGenic=IsGenicStart"),'LengthBucket','IsGenic','Length Distribution by isGenicStart',useLogY =F)

# INDIVIDUAL SIGNATURE CHECK
filter = cohortSummary(merge(svData,svDriverAndGermline %>% filter(gene=='CCNE1') %>% select(SampleId=sampleId),by='SampleId'),"",'SampleId') %>% top_n(20,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_INT_TI'&Type=='INV'),SampleId %in% filter","LengthBucket,IsGenic=IsGenicStart"),'LengthBucket','IsGenic','Length Distribution by isGenicEnd',useLogY =F)

#3. Zoom in on FragileSites
plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_INT_TI'&Type=='INV'),IsFS","LengthBucket,GeneStart"),'LengthBucket','GeneStart','Length Distribution by isGenicStart',useLogY =F)

#4.
View(svLinks %>% group_by(ResolvedType,isGenic=GeneStart!=''&GeneEnd!='') %>% count() %>%spread(isGenic,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))

############################################################
#############Driver vs Signature correlation ###############
############################################################
#### TO DO RERUN!!!!!
driverSignatureCoocurrence = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/coc_sv_results_extended.csv') 
View(driverSignatureCoocurrence %>% arrange(FDR) %>% filter(CancerType!='All',Category!='SIMPLE_SV',FDR<0.1))
write.csv(driverSignatureCoocurrence %>% arrange(FDR) %>% filter(CancerType!='All',Category!='SIMPLE_SV',FDR<0.1),'~/driverSignatureCoocurrence.csv')

cancerTypeLOHCoocurrence = read.csv('~/hmf/analyses/SVAnalysis/coc_sv_loh_cancer.csv') 
View(cancerTypeLOHCoocurrence %>% group_by(CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))
write.csv(cancerTypeLOHCoocurrence %>% group_by(CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR),'~/cancerTypeLOHCoocurrence.csv')
driverGeneLOHCoocurrence = read.csv('~/hmf/analyses/SVAnalysis/coc_gene_loh_results_extended.csv') 
View(driverGeneLOHCoocurrence %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,positive=CountGtExp,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))
write.csv(driverGeneLOHCoocurrence %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR),'~/driverGeneLOHCoocurrence.csv')


####################
# TO DO

View(svLinks %>% group_by(Id1,Id2) %>% tally() %>% filter(n>1) )   ##### CHECK THIS WITH CHALRES?

ggplot(data=svData %>% filter(ResolvedType=='SIMPLE',PosEnd-PosStart>20000,PosEnd-PosStart<500000,Type=='DEL') %>% group_by(SampleId,IsFS) %>% tally() %>% spread(IsFS,n,fill=0)) + 
  stat_ecdf(aes(`TRUE`,color='isFS'),geom = "step", pad = FALSE) + stat_ecdf(aes(`FALSE`,color='FSLength non-FS'),geom = "step", pad = FALSE) + 
  ylim(0,1) + labs(title = 'CDF # of FS length DELs')

#CHECK
scatterPlot(svData %>% filter(ResolvedType=='DEL_EXT_TI',Type=='BND') %>% mutate(diff=SynDelDupLen-SynDelDupAvgTILen),'SynDelDupLen','diff',F,F) 
scatterPlot(svData %>% filter(ResolvedType=='DEL_EXT_TI',Type!='BND'),'SynDelDupLen','SynDelDupAvgTILen',F,F) 
scatterPlot(svData %>% filter(ResolvedType=='DUP_EXT_TI',Type=='BND'),'SynDelDupLen','SynDelDupAvgTILen',F,F) 


