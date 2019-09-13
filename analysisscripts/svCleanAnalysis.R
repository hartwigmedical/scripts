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
    InsLenBucket=2**(round(log(nchar(InsertSeq),2),0)),
    QualScoreBucket=2**(round(log(nchar(QualScore),2),0))
  )
}

createsvLinksBuckets <- function(svLinks) {
  svLinks %>% mutate(
    TILengthBucket=ifelse(!TILength>0,TILength,2**(pmin(25,pmax(5,round(log(TILength,2),0))))),
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


#############################################
################ LOADING ####################
#############################################

PATH='~/hmf/analyses/SVAnalysis/'

svData = read.csv(paste(PATH,'SVA_SVS.csv',sep=''), header = T, stringsAsFactors = F)
#svClusterHistory = (read.csv(paste(PATH,'SVA_CLUSTERING_HISTORY.csv',sep=''), header = T, stringsAsFactors = F))
svClusters = (read.csv(paste(PATH,'SVA_CLUSTERS.csv',sep=''), header = T, stringsAsFactors = F))
svDisruptions = (read.csv(paste(PATH,'SVA_DISRUPTIONS.csv',sep=''), header = T, stringsAsFactors = F))
svFusions = (read.csv(paste(PATH,'SVA_FUSIONS.csv',sep=''), header = T, stringsAsFactors = F))
svDrivers = read.csv(paste(PATH,'SVA_DRIVERS.csv',sep=''), header = T, stringsAsFactors = F) #%>% distinct(Gene,SampleId,ClusterId,Category,EventType,ResolvedType,Chromosome,Arm,MatchInfo,GeneMinCN,) 
svDoubleMinute = read.csv(paste(PATH,'SVA_DM.csv',sep=''), header = T, stringsAsFactors = F)# %>% distinct(SampleId,ClusterId,DMSvTypes,DMSvCount,FullyChained,SvIds,ChainLength,Chromosomes,DupPosStart,DupPosEnd,MaxCopyNumber,AmpGenes)
svLinks = read.csv(paste(PATH,'SVA_LINKS.csv',sep=''), header = T, stringsAsFactors = F)
sampleCancerTypes= (read.csv(paste(PATH,'sample_cancer_types.csv',sep=''), header = T, stringsAsFactors = F))
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
svData=merge(svData,svClusters %>% separate('Annotations',c('SynLength1','SynLength2','SynGapLength'),sep=';') %>% select(SampleId,ClusterId,Subclonal,Foldbacks,Synthetic,SynLength1,SynLength2),by=c('SampleId','ClusterId'))   #TO DO: improve synLength annotations
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

# BE data
beData = rbind(svData %>% mutate(Arm=ArmStart,LnkLen=LnkLenStart,LocTopTI=LocTopTIStart,LocTopType=LocTopTypeStart,LocTopId=LocTopIdStart,Chr=ChrStart,Pos=PosStart,Orient=OrientStart,IsStart=T,Anchor=AnchorStart,
              RefContext=RefContextStart,LE=LEStart,DBLength = DBLenStart,Assembled = ifelse(grepl('asm',AsmbStart),"Assembled","NotAssembled")),
              svData %>% mutate(Arm=ArmEnd,LnkLen=LnkLenEnd,LocTopTI=LocTopTIEnd,LocTopType=LocTopTypeEnd,LocTopId=LocTopIdEnd,Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,IsStart=F,Anchor=AnchorEnd,
              RefContext=RefContextEnd,LE=LEEnd,DBLength = DBLenEnd, Assembled = ifelse(grepl('asm',AsmbStart),"Assembled","NotAssembled"))) %>%
              select(SampleId,Id,IsStart,ClusterId,CancerType,Type,ResolvedType,Subclonal,Arm,LnkLen,LocTopTI,LocTopType,LocTopId,Chr,Pos,Orient,RefContext,LE,DBLength,Assembled,Anchor,IsPolyA,ClusterCount,PloidyMin,PloidyMax,ClusterCountBucket,ClusterDesc,Synthetic,LengthBucket,IsFoldBack) %>% 
              mutate(DBLenBucket = ifelse(DBLength==0,0,ifelse(DBLength<0,-(2**round(log(-DBLength,2))),2**round(log(DBLength,2)))))

#############################################
########## OVERVIEW OF COHORT ###############
#############################################

# 1. By RESOLVED Type
View(svData %>% group_by(ResolvedType,Type) %>% tally() %>% spread(Type,n))

View(svData %>% filter(Type=='INF',grepl('Prox',ClusterReason)) %>% group_by(ClusterDesc) %>% count) 


View(svData %>% filter(Subclonal=='false') %>% group_by(ResolvedType,LengthBucket=LengthBucket>5E6) %>% tally() %>% spread(LengthBucket,n))

#2. Synthetic vs non synthetic
View(merge(svData %>% filter(Subclonal=='false',!ResolvedType %in% c('INF','POLY_G_C'),!grepl('INF',ClusterDesc)|Synthetic=='false') %>% group_by(ResolvedType,Synthetic) %>% tally() %>% spread(Synthetic,n),
           svClusters %>% filter(Subclonal=='false') %>% group_by(ResolvedType,Synthetic) %>% tally() %>% spread(Synthetic,n),
           by='ResolvedType',suffixes=c('.variants','.clusters',all.x=T)))

View(svData %>% filter(ResolvedType=='DEL') %>% group_by(ClusterDesc) %>% count)

# 2. By CLUSTERCOUNTBUCKET
View(svData %>% filter(Subclonal=='false') %>% group_by(ClusterCountBucket,ResolvedType) %>% tally() %>% spread(ClusterCountBucket,n))

#3. Simple SV make up 70% of clusters but only 30% of SVs
View(svClusters %>% filter(Subclonal=='false') %>% group_by(ResolvedType) %>% summarise(variants=sum(ClusterCount),clusters=n()))

###########################################################
############### SIMPLE vs COMPLEX #########################
###########################################################

plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',Type %in% c('DEL','DUP','INV'),!IsFoldBack,ResolvedType=='COMPLEX'|(ClusterCount==1&(ResolvedType=='DEL'|ResolvedType=='DUP')),ArmStart==ArmEnd","LengthBucket,ResolvedType"),'LengthBucket','ResolvedType','Length Distribution of local variants by resolved type (excluding FB INV)',useLogY =F)

#############################################
############### DUP #########################
#############################################

#0. Length distribution
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType %in% c('DUP'),ClusterCount==1","LengthBucket,ClusterCount"),'LengthBucket','ClusterCount','Length Distribution of Tandem DUP',useLogY =F)
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',Type %in% c('DEL','DUP','INV'),!IsFoldBack,ResolvedType=='COMPLEX'|(ClusterCount==1&(ResolvedType=='DEL'|ResolvedType=='DUP')),ArmStart==ArmEnd","LengthBucket,ResolvedType"),'LengthBucket','ResolvedType','Length Distribution of local variants by resolved type (excluding FB INV)',useLogY =F)

topSamplesCount=100
#1. NOVEL SHORT DUP SIG
top100shortDUP = (svData %>% filter(ResolvedType=='DUP',ClusterCount==1,Length<1e3) %>% group_by(SampleId,CancerType) %>% count() %>% as.data.frame()
                  %>% top_n(topSamplesCount,n) %>% group_by(CancerType) %>% count)
View(merge(top100shortDUP,sampleList %>% group_by(CancerType) %>% summarise(total=n()),by='CancerType',all=T) %>% 
       mutate(rate=round(nn/total,3),p=ppois(nn,topSamplesCount*total/nrow(sampleList),FALSE)) %>% filter(total>0) %>% arrange(p,-total))

#1. MEDIUM DUP
top100medDUP = (svData %>% filter(ResolvedType=='DUP',ClusterCount==1,Length>1e3,Length<1e5) %>% group_by(SampleId,CancerType) %>% count() %>% as.data.frame()
                  %>% top_n(topSamplesCount,n) %>% group_by(CancerType) %>% count)
View(merge(top100medDUP,sampleList %>% group_by(CancerType) %>% summarise(total=n()),by='CancerType',all=T) %>% 
       mutate(rate=round(nn/total,3),p=ppois(nn,topSamplesCount*total/nrow(sampleList),FALSE)) %>% filter(total>40) %>% arrange(p,-total))
#1. LONG DUP
top100LongDUP = (svData %>% filter(ResolvedType=='DUP',ClusterCount==1,Length>1e5) %>% group_by(SampleId,CancerType) %>% count() %>% as.data.frame()
          %>% top_n(topSamplesCount,n) %>% group_by(CancerType) %>% count)
View(merge(top100LongDUP,sampleList %>% group_by(CancerType) %>% summarise(total=n()),by='CancerType',all=T) %>% 
       mutate(rate=round(nn/total,3),p=ppois(nn,topSamplesCount*total/nrow(sampleList),FALSE)) %>% filter(total>40) %>% arrange(p,-total))


############################################################
########### Short TIs #####################################
############################################################
#0. Numbers per sample

ggplot(data=merge(svData %>% group_by(SampleId) %>% tally(),svLinks %>% group_by(SampleId,shortTI=ifelse(TILength<1000,'ShortTICount','LongTICount')) %>% tally() %>% spread(shortTI,n),by='SampleId'),aes(n,ShortTICount)) + geom_point()

#1. TIs length distibution, showing 2 peaks for most types ('short', ie <1k and long -> presumed to be cause by random length of distance between breakages)
ggplot(data=svLinks %>% filter(ResolvedType=='COMPLEX') %>% group_by(TILengthBucket,ResolvedType) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_x_log10()+ facet_wrap(~ResolvedType)
#2. SHORT TI actual shape
ggplot(data=svLinks  %>% filter(TILength<1000) %>% mutate(TILengthBucket=round(TILength+0.001,-1)) %>% group_by(TILengthBucket) %>% summarise(count=n()),
       aes(x=TILengthBucket))+geom_line(aes(y=count,colour='count'))+scale_y_log10()

#4. Short TIs create 'synthetic' events that match the simpleSV event profile per sample
scatterPlot(svData %>% filter(ResolvedType=='DEL',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # DEL correlated with synthetic DEL
scatterPlot(svData %>% filter(ResolvedType=='DEL',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc),(Synthetic=='false'&Length<5e3)|(Synthetic=='true'&SynLength1<5e3)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # Short DEL correlated with short synthetic DEL
scatterPlot(svData %>% filter(ResolvedType=='DEL',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc),(Synthetic=='false'&Length>5e3)|(Synthetic=='true'&SynLength1>5e3)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # Short DEL correlated with short synthetic DEL
scatterPlot(svData %>% filter(ResolvedType=='DUP',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # DUP correlated with synthetic DUP
scatterPlot(svData %>% filter(ResolvedType=='DUP',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc),(Synthetic=='false'&Length<8e4)|(Synthetic=='true'&SynLength2<8e4)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # Short DEL correlated with short synthetic DUP
scatterPlot(svData %>% filter(ResolvedType=='DUP',!grepl('SGL',ClusterDesc),!grepl('INF',ClusterDesc),(Synthetic=='false'&Length>8e4)|(Synthetic=='true'&SynLength2>8e4)) %>% group_by(SampleId,Synthetic)
            %>% count() %>% spread(Synthetic,n),'false','true',F,F) # Short DEL correlated with short synthetic DEL
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DEL','DUP',F,F) # Simple DEL and DUP not highly correlated
scatterPlot(svData %>% filter(Synthetic=='true') %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DEL','DUP',F,F)  # Synth DEL and DUP also not


plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ClusterDesc=='INV=2',ClusterCount==2",'ResolvedType,LengthBucket'),'LengthBucket','ResolvedType','Synthetic DEL and DUP raw lengths',useLogY =F )


plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='DEL'|ResolvedType=='DUP',ClusterCount==2,Type!='BND',(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd))",'ClusterDesc,LengthBucket'),'LengthBucket','ClusterDesc','Synthetic DEL and DUP raw lengths',useLogY =F )
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='DEL'|ResolvedType=='DUP',ClusterCount==2,Type!='BND',((DBLenStart<0&DBLenStart>-1000)|(DBLenEnd<0&DBLenEnd>-1000))",'ClusterDesc,LengthBucket'),'LengthBucket','ClusterDesc','Synthetic DEL and DUP raw lengths',useLogY =F )
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='DEL'|ResolvedType=='DUP',ClusterCount==2,Type!='BND',((DBLenStart<=-1000)|(DBLenEnd<=-1000))",'ClusterDesc,LengthBucket'),'LengthBucket','ClusterDesc','Synthetic DEL and DUP raw lengths',useLogY =F )


##############################################################
########## RECIPROCAL EVEBTS ####################
#############################################################
#1. Reciprocal DUPS Correlated with DUPS
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DUP','RECIP_INV_DUPS',F,F)  # Recip DUP correlated with DUP
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DUP','RECIP_TRANS_DUPS',F,F)  # Recip DUP correlated with DUP
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'RECIP_INV_DUPS','RECIP_TRANS_DUPS',F,F)  # Recip DUP correlated with DUP

#2. Reciprocal DUPSS Not correlated with other reciprocal types
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DEL','RECIP_INV_DUPS',F,F)  # Recip DUP not correlated with RECIP_INV
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'DEL','RECIP_TRANS_DUPS',F,F)  # Recip DUP not correlated with RECIP_INV
scatterPlot(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n),'RECIP_INV','RECIP_INV_DUPS',F,F)  # Recip DUP not correlated with RECIP_INV

#3. Reciprocal INV NOT enriched in typical DEL lengths OR FS!
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='DEL'|ResolvedType=='RECIP_INV'",'IsFS,LengthBucket'),'LengthBucket','IsFS','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='DEL'|ResolvedType=='RECIP_INV',ClusterCount<=2",'ClusterDesc,LengthBucket'),'LengthBucket','ClusterDesc','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )
plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',Synthetic=='true',ResolvedType=='DEL'|ResolvedType=='RECIP_INV'",'ResolvedType,LengthBucket'),'LengthBucket','ResolvedType','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )
#scatterPlot(svSampleSummary,'countSimpleDELShortMed','countRECIP_TRANS',F,F)  # some short DELs can be correlated with RECIP_TRANS
#scatterPlot(svSampleSummary,'countSimpleDELShortMed','countRECIP_INV',F,F)  # some short DELs can be correlated with RECIP_TRANS
#scatterPlot(svSampleSummary,'countSimpleDELLong','countRECIP_TRANS',F,F)  # Long DELs NOT correlated with RECIP_TRANS

plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',ResolvedType=='RECIP_TRANS'",'IsFS,LengthBucket'),'LengthBucket','IsFS','Reciprocal Inversions and SimpleSV by IsFS (log scale)',useLogY =T )


# TO DO:  show the length distribution and come up with a standardised definition
ggplot(data=svData %>% filter(Subclonal=='false',Synthetic=='false',grepl('RECIP',ResolvedType)) %>% group_by(NearestLenBucket,ResolvedType) %>% count %>% spread(ResolvedType,n,fill=0),aes(x=NearestLenBucket))+
  geom_line(aes(y=RECIP_TRANS,colour="RECIP_TRANS"))+geom_line(aes(y=RECIP_TRANS_DEL_DUP,colour="RECIP_TRANS_DEL_DUP"))+geom_line(aes(y=RECIP_TRANS_DUPS,colour="RECIP_TRANS_DUPS"))+
  geom_line(aes(y=RECIP_INV,colour="RECIP_INV"))+geom_line(aes(y=RECIP_INV_DEL_DUP,colour="RECIP_INV_DEL_DUP"))+geom_line(aes(y=RECIP_INV_DUPS,colour="RECIP_INV_DUPS"))+
  scale_x_log10()
View(svData %>% filter(Subclonal=='false',Synthetic=='false',grepl('RECIP',ResolvedType)) %>% group_by(NearestLenBucket,ResolvedType) %>% count %>% spread(ResolvedType,n,fill=0))

##############################################################
########## Short Deletion Bridge Analysis ####################
#############################################################

## NB - Might have to reapply this filter (filter(DBLength>-1000))
     
#1.DBLength by ResolvedType => less than 50 bases (NB - DBLength of 1 means exact break - shoould correct this)
# Note the LINE double peak and the sharp feature for Reciprocal Inversion, Reciprocal Translocations and None
# NB LE=='None' excludes line source elements with complex topology
plot_count_by_bucket_and_type(cohortSummary(beData, "DBLength<=50,DBLength>=-50,Subclonal=='false',LE=='None'",'DBLength,ResolvedType'),
                              'DBLength','ResolvedType','DB Length for selected resolved types(<=100bases)',useLogX = F,useLogY = T)


plot_count_by_bucket_and_type(cohortSummary(beData, "DBLength<=800,DBLength>=-800,Subclonal=='false',LE=='None',ResolvedType=='COMPLEX'",'DBLength=round(DBLength,-1),ResolvedType'),
                              'DBLength','ResolvedType','DB Length for selected resolved types(<=100bases)',useLogX = F,useLogY = T)




#2. EXACT BREAKS
# RI have a sharp feature for longer INV
plot_count_by_bucket_and_type(cohortSummary(beData,"DBLength>-50,DBLength<=50,Subclonal=='false',ResolvedType=='RECIP_INV',Synthetic=='false'",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for INV in RI by length bucket',useLogX = F,useLogY = F)
# Sharp RI feature exists for Complex short INV also. Excess compared to DEL & DUP and Mostly around 60-500 bases (and not a foldback feature)
plot_count_by_bucket_and_type(cohortSummary(beData,"DBLength>-50,DBLength<=50,Subclonal=='false',ResolvedType=='COMPLEX'",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for short variants in COMPLEX clusters by length bucket',useLogX = F,useLogY = T)
# Structures with 2 local TIs = likely to be sharp break for INV!!!!!
plot_count_by_bucket_and_type(cohortSummary(beData, "DBLength<=50,DBLength>=-50,Subclonal=='false',ResolvedType=='COMPLEX',LocTopTI<6",'DBLength,LocTopTI'),
                              'DBLength','LocTopTI','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY =F)
## RERUN THIS ONE AFTER RECLASSIFICATION OF 2 TI events
plot_count_by_bucket_and_type(cohortSummary(beData, "DBLength<=50,DBLength>=-50,Subclonal=='false',ResolvedType=='COMPLEX',LocTopTI==2",'DBLength,LocTopType'),
                              'DBLength','LocTopType','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY =F)


#3 SYNTH_DUP which are likely LINE
plot_count_by_bucket_and_type(cohortSummary(beData %>% filter(ResolvedType=='DUP'),
                                            "DBLength>-50,DBLength<=50,Synthetic=='true'",'DBLength,IsPolyA'),
                              'DBLength','IsPolyA','SyntheticDelDupLength',useLogX = F,useLogY = F)

#4. SOME SINGLES in COMPLEX 2 clusters and above LOOK LIKE LINE
plot_count_by_bucket_and_type(cohortSummary(beData, "DBLength<=5,DBLength>=-100,Subclonal=='false',LocTopType=='DSB',LE=='None',ClusterCount==2,ResolvedType=='COMPLEX'",'DBLength,ClusterDesc'),
                              'DBLength','ClusterDesc','DBLengthByResolvedType(<=100bases)',useLogX = F,useLogY = F)


#5. LINE Elements - no clear indication of why there are 2 DB peaks.
#The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific.   Also appars to be the same regardless of source element
print(ggplot(data = beData %>% filter(DBLength<=50,ResolvedType %in% c('LINE')) %>% group_by(CancerType,OLPeak=DBLength<(-7)) %>% tally(), aes(x = reorder(CancerType, -n), y = n, fill =OLPeak))
      + geom_bar(stat = "identity", colour = "black") + ylab("DB Count") + xlab("Tumor Type") + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15)) + labs(title = "DB Counts for LINE Elements by cancerType and Overlapping Peak"))
# Also does not seem to depend significantly on refContext.  
View(beData %>% filter(ResolvedType=='LINE',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context,longOverlap=DBLength<(-7)) %>% tally() %>% arrange(-n) %>% spread(longOverlap,n))


################################################
############## FOLDBACKS #######################
################################################

View(svClusters %>% filter(Foldbacks>0) %>% group_by(Subclonal,Foldbacks) %>% tally %>% spread(Subclonal,n))

#1.Foldback length distribution for simple inversions short chained foldbacks
print(ggplot(data = foldbacks %>% filter(Subclonal=='false') %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
             aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=INV, colour='Simple')) + geom_line(aes(y=Combo, colour='Synthetic' )) 
      + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))


# TO DO: WHAT DO WE WANT TO DO WITH OUR DETAILED FB analysis?
# B. When we have 2 FB how often are they facing each other?
# C. How often is the most telomeric FB facing to the centromere
View(merge(foldbacks %>% filter(Subclonal=='false') %>% group_by(SampleId,ChrStart,ArmStart) %>% count %>% group_by(SampleId) %>% count,sampleList,by='SampleId',all = T))
#2. Per sample counts (needs work)
temp = merge(foldbacks %>% filter(Subclonal=='false') %>% group_by(SampleId,ClusterId) %>% count %>% group_by(SampleId) %>% count,sampleList,by='SampleId',all = T)%>% replace_na(list(nn=0))
temp$nn = temp$nn %>% replace_na(0)
ggplot(data=merge(foldbacks %>% filter(Subclonal=='false') %>% group_by(SampleId,ClusterId) %>% count %>% group_by(SampleId) %>% count,sampleList,by='SampleId',all = T)%>% replace_na(list(nn=0)),aes(nn)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # Clusters with Foldbacks  per Sample')
ggplot(data=svClusters %>% group_by(SampleId) %>% summarise(Foldbacks=sum(Foldbacks)),aes(Foldbacks)) + stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # Foldbacks  per Sample')

merge(foldbacks %>% filter(Subclonal=='false') %>% group_by(SampleId,ClusterId) %>% count %>% group_by(SampleId) %>% count,sampleList,by='SampleId',all = T,fill=0)
# LOCATION
telomereOrientedFB=(svData %>% filter(CancerType!='ABreast',Subclonal=='false',IsFoldBack,Type=='INV',PosEnd-PosStart<1e4) %>% mutate(nrows=n(),FBOrient=ifelse(OrientStart*ifelse(ArmStart=='P',1,-1)==1,'TELOMERE','CENTROMERE')) %>% filter(FBOrient=="TELOMERE") %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% summarise(n=n(),avgRep=mean(RepOriginStart)) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
centromereOrientedFB=(svData %>% filter(CancerType!='ABreast',Subclonal=='false',IsFoldBack,Type=='INV',PosEnd-PosStart<1e4) %>% mutate(nrows=n(),FBOrient=ifelse(OrientStart*ifelse(ArmStart=='P',1,-1)==1,'TELOMERE','CENTROMERE')) %>% filter(FBOrient=="CENTROMERE") %>% group_by(chr=ChrStart,pos=round(PosStart,-6),nrows) %>% summarise(n=n(),avgRep=mean(RepOriginStart)) %>% mutate(p=ppois(n,nrows/2.8e9*1e6,FALSE),q=p.adjust(p,"BH",2.8e9/1e6)))
ggplot() + geom_point(data=telomereOrientedFB,aes(pos,n,color="TELOMERE"),shape="x") + geom_point(data=centromereOrientedFB,aes(pos,n,color="CENTROMERE"),shape="x") + facet_wrap(~chr)

#######################################################
################# LINE #############################
#######################################################

# 1. INSERTION MOTIF: ~15% are A-TTTTT with another 30% or so a similar variant
View(beData %>% filter(Subclonal=='false',ResolvedType=='LINE',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
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

View(LineClusters)
View(LineClusters %>% filter(countPOLYA/n<0.3))
scatterPlot(LineClusters%>% filter(countAllSuspected==0),'countPOLYA','countKnown',F,F)
scatterPlot(LineClusters %>% mutate(n=pmin(n,100),countAllSuspected=pmin(100,countAllSuspected+countKnown)),'n','countAllSuspected',F,F)
scatterPlot(LineClusters,'countAllSuspected','countPOLYA',F,F)
#################################################
################# COMPLEX CLUSTER STATS #########
#################################################

#1 LOCAL - REMOTE CONSISTENCY
clusterArmStats = createClusterArmStats(svData)
View(clusterArmStats %>% filter(LocalCount+RemoteCount+SGLCount>10,ResolvedType=='COMPLEX'))

# 2. Per arm by ClusterCount
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='LINE') %>% group_by(SampleId,chr_arm,CC=pmin(ClusterCount,8),id=ClusterId) %>% tally(),
           svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='LINE') %>%group_by(SampleId,chr_arm,CC=pmin(ClusterCount,8),id=ClusterId) %>% tally()) %>% 
       group_by(SampleId,chr_arm,CC,id) %>% tally() %>%
       group_by(SampleId,CC,chr_arm) %>% tally() %>% spread(CC,n) %>%
       filter(chr_arm!='0_P'))

#3. LARGE Clusters - # arms not super correlated with # of SVs
scatterPlot(svClusters %>% filter(ClusterCount>10,ResolvedType!='LINE') %>% mutate(realArmCount=ArmCount-FragmentArms),'realArmCount','ClusterCount',F,F)

View(svClusters %>% filter(ClusterCount>20,ResolvedType!='LINE')  %>% mutate(realArmCount=ArmCount-FragmentArms) %>% filter(realArmCount>10))

#4. CDF of large non-line clusters per sample
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>10,ResolvedType!='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>10 Per Sample')
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>50,ResolvedType!='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>50 Per Sample')

#######################################################
########## Pseudogene insertions ###################
#######################################################

svPseudoGeneLinks = svLinks %>% filter(ExonMatch!='') %>% separate(ExonMatch,c('transcriptId','exonRank','exonLength','other')) #what is the last field in the ExonMatch?

#1. 445 events
View(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd,transcriptId) %>% summarise(count=n(),minExon=min(as.numeric(exonRank)),maxExon=max(as.numeric(exonRank)),range=maxExon-minExon+1))

#2. found in 6% of samples - nearly all with highly activated LINE elements
View(merge(merge(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd) %>% count() %>% group_by(SampleId) %>% count(),sampleList,by='SampleId',all=T),
           svClusters %>% filter(ResolvedType=='LINE') %>% group_by(SampleId) %>% summarise(numLine=n()),by='SampleId',all.x = T) %>% arrange(-nn))

#2. found in 6% of samples - nearly all with highly activated LINE elements
temp = (merge(svLinks %>% filter(ExonMatch!='') %>% group_by(SampleId,ClusterId,GeneStart,GeneEnd) %>% count() %>% group_by(SampleId) %>% summarise(pseudogeneINSCount=n()),svData %>% group_by(SampleId) %>% summarise(count=n(),sumLine=sum(ResolvedType=='Line')),by='SampleId',all=T,fill=0))
View(temp %>% group_by(hasLine=sumLine>0,hasPesudogenINS=!is.na(pseudogeneINSCount)) %>% count %>% spread(hasLine,n))
scatterPlot(temp,'countLine','pseudogeneINSCount',F,F)

##TO DO: Create circos charts with gene transcripts for further analysis

View(svPseudoGeneLinks %>% group_by(SampleId,ClusterId,ClusterCount,ClusterDesc) %>% count())


pgLinks = read.csv('~/Documents/pseudo_gene_links.csv')
pgLinksStart = pgLinks %>% select(SampleId,SvId=Id1,HomOffset=StartHomOffset,PosOffset=StartPosOffset)
pgLinksEnd = pgLinks %>% select(SampleId,SvId=Id2,HomOffset=EndHomOffset,PosOffset=EndPosOffset)
pgLinksCombined = rbind(pgLinksStart,pgLinksEnd)

View(pgLinks %>% mutate(exactMatch=StartPosOffset==0&EndPosOffset==0) %>% filter(exactMatch==F))
View(pgLinks %>% mutate(exactMatch=StartPosOffset==0&EndPosOffset==0,partialExon=(StartPosOffset<0|EndPosOffset<0),partialIntron=(StartPosOffset>0|EndPosOffset>0)) %>% 
       group_by(SampleId,ClusterId,Gene,TransId,ClusterCount) %>% summarise(totalLength=sum(TILength),matchLength=sum(ifelse(exactMatch,TILength,0)),count=n(),countExactMatch=sum(ifelse(exactMatch,1,0)),countPartial=sum(ifelse(partialExon,1,0)),countIntronic=sum(ifelse(partialIntron,1,0)),minExon=min(as.numeric(ExonRank)),maxExon=max(as.numeric(ExonRank)),range=maxExon-minExon+1)
      %>% mutate(diff=ClusterCount-count))

View(pgLinks %>% mutate(exactMatch=StartPosOffset==0&EndPosOffset==0,partialExon=(StartPosOffset<0|EndPosOffset<0),partialIntron=(StartPosOffset>0|EndPosOffset>0)) %>% 
       group_by(SampleId) %>% summarise(totalLength=sum(TILength),matchLength=sum(ifelse(exactMatch,TILength,0)),count=n(),countExactMatch=sum(ifelse(exactMatch,1,0)),countPartial=sum(ifelse(partialExon,1,0)),countIntronic=sum(ifelse(partialIntron,1,0)),minExon=min(as.numeric(ExonRank)),maxExon=max(as.numeric(ExonRank)),range=maxExon-minExon+1)
     )

homOffsetMatches = pgLinksCombined %>% group_by(SampleId,SvId) %>% 
  summarise(HomOffsetUsed=sum(PosOffset==0),HomOffset1=first(HomOffset),HomOffset2=last(HomOffset))

homOffsetMatches = homOffsetMatches %>% filter(HomOffsetUsed==2) %>% mutate(HomOffMatched=HomOffset1==HomOffset2)
View(homOffsetMatches)

View(homOffsetMatches %>% group_by(HomOffMatched) %>% count())

##############################################################
########## FRAGILE SITE ######################################
##############################################################

#1. Fragile Site Length Distribution
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType %in% c('SIMPLE','RECIP_INV'),Type!='BND',Type!='SGL'","LengthBucket,IsFS"),'LengthBucket','IsFS','Length Distribution of SimpleSV & Reciprocal Inversion by IsFS',useLogY =T)

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
scatterPlot(svData %>% filter(ResolvedType=='DEL_EXT_TI',Type=='BND') %>% mutate(diff=SyntheticLen-SyntheticTILen),'SyntheticLen','diff',F,F) 
scatterPlot(svData %>% filter(ResolvedType=='DEL_EXT_TI',Type!='BND'),'SyntheticLen','SyntheticTILen',F,F) 
scatterPlot(svData %>% filter(ResolvedType=='DUP_EXT_TI',Type=='BND'),'SyntheticLen','SyntheticTILen',F,F) 

plot_count_by_bucket_and_type(cohortSummary(svData,"Subclonal=='false',((ResolvedType=='DEL'|ResolvedType=='DUP')&ClusterCount==1)",'IsFS,LengthBucket'),'LengthBucket','IsFS','SimpleSV by Is Fragile Site (log scale)',useLogY =F )


filter = cohortSummary(svData,"ResolvedType=='DUP'",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData %>% filter(SampleId %in% filter),"Subclonal=='false',((ResolvedType=='DEL'|ResolvedType=='DUP')&ClusterCount==1)|ResolvedType=='RECIP_INV'",'SampleId,CancerType,LengthBucket')
                              %>%  mutate(ID = paste(CancerType,SampleId)),'LengthBucket','ID','SimpleSV by Is Genic (log scale)',useLogY =F )

View(svData %>% filter(Type=='SGL',CNChgStart/CNStart>0.1,CNStart>20) %>% group_by(SampleId) %>% count)
View(svData %>% filter(SampleId=='CPCT02330132T',Type=='SGL'))
View(svData %>% filter(Type=='SGL',CNChgStart/CNStart<0.1) %>% group_by(SampleId,CN=pmin(50,round(CNStart,-1))) %>% count() %>% spread(CN,n))

View(svData %>% filter(Type=='SGL',CNChgStart/CNStart<0.1) %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n))
View(svData %>% filter(ResolvedType!='LINE',Type=='SGL',CNChgStart/CNStart<=0.1,CNChgStart<0.5) %>% group_by(IsPolyA) %>% count())
View(merge(merge(svData %>% filter(ResolvedType!='LINE',Type=='SGL',CNChgStart/CNStart>=0.1) %>% group_by(SampleId) %>% count(),
     svData %>% filter(ResolvedType=='LINE',CNChgStart/CNStart>0.1) %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T),
     svData %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T))


View(merge(svClusters %>% filter(ClusterDesc=='BND=3',AcTIOnly==0,TotalTIs==0),svLinks %>% distinct(SampleId,ClusterId,cc=ClusterCount),by=c('SampleId','ClusterId'),all.x=T))# %>% filter(is.na(cc)))
View(svData %>% filter(SampleId=='CPCT02010255TII',ClusterId=='9'))

View(svData %>% group_by(Type,Subclonal) %>% count() %>% spread(Type,n))

View(svData %>% group_by(ResolvedType,Subclonal) %>% count() %>% spread(Subclonal,n))


############# SUBCLONAL INVESTIGATION   #############

#1.Foldbacks < 100 bases look like artefacts
print(ggplot(data = foldbacks %>% filter(FoldbackType=='INV') %>% group_by(FoldbackLenBucket,Subclonal) %>% summarise(Count=n()) %>% spread(Subclonal,Count), 
             aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=true, colour='Subclonal')) + geom_line(aes(y=false, colour='Clonal' )) 
      + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))

# concentrated in a handful of samples

View(foldbacks %>% filter(FoldbackType=='INV',Length<100,NearestLen>-1) %>% group_by(Subclonal,SampleId) %>% count() %>% spread(Subclonal,n))  # LOOK at other damage in 2 samples
View(foldbacks %>% filter(FoldbackType=='INV',NearestLen>-1) %>% group_by(Subclonal,CNChgStart<0.5&CNChgEnd<0.5,Length<100) %>% count() %>% spread(Subclonal,n))  # LOOK at other damage in 2 samples
View(foldbacks %>% filter(FoldbackType=='INV',SampleId!='CPCT02050114T'|SampleId!='CPCT02010304T') %>% group_by(Subclonal,Length<100,roundedVaf=round(Ploidy*2/pmax(CNStart,CNEnd),1)/2) %>% count() %>% spread(Subclonal,n))

## ==> Proposed RULE - exlcude FB < 100 bases IF (CNChngStart<0.5&CNChgEnd<0.5) or VAF < 0.2 

#2. Unbalanced translocation
View(svData %>% filter(Type=='BND') %>% group_by(ResolvedType=='UNBAL_TRANS',Subclonal,roundedVaf=round(Ploidy*2/pmax(CNStart,CNEnd),1)/2) %>% count() %>% spread(Subclonal,n))

#Frequentl artefacts at long homology for isolated or subclonal BND
View(svData  %>% filter(ResolvedType=='UNBAL_TRANS',Type=='BND',grepl('ACACAC',Homology)|grepl('GTGTGT',Homology)|grepl('TTTTTT',Homology)|grepl('AAAAAA',Homology)|grepl('ATATAT',Homology)) %>% group_by(SampleId,Subclonal) %>% count %>% spread(Subclonal,n))

# handful of very high damage samples, particularly CPCT02050114T
View(svData %>% filter(ResolvedType=='LINE'|ResolvedType=='UNBAL_TRANS') %>% unite(res_sub,ResolvedType,Subclonal) %>% group_by(SampleId,res_sub) %>% count() %>% spread(res_sub,n,fill=0) %>% arrange(-UNBAL_TRANS_true))

# Around half of  subclonal unbalanced translocation are likely to be LINE insertions.  Rest are likely FP
View(svData %>% filter(ResolvedType=='LINE'|(ResolvedType=='UNBAL_TRANS')) %>% unite(res_sub,ResolvedType,Subclonal) %>% group_by(SampleId,res_sub) %>% count() %>% spread(res_sub,n,fill=0) %>% group_by(LINE_false>10) %>% summarise(sum(UNBAL_TRANS_true),sum(UNBAL_TRANS_false)) )

View(svData  %>% group_by(SampleId,ResolvedType,Subclonal) %>% count() %>% spread(ResolvedType,n) %>% select(SampleId,UNBAL_TRANS,LINE,everything()) %>% arrange(-UNBAL_TRANS))

## ==> Proposed RULE - Filter BND if nearest variant at both breakends > 5000 bases and (CNChngStart<0.5&CNChgEnd<0.5)    

#3. INF
View(svData %>% filter(grepl('INF',ClusterDesc),ResolvedType=='PAIR_OTHER',grepl('Prox',ClusterReason)) %>% group_by(Type,CNChgStart<0.25|CNChgEnd<0.25) %>% count %>% spread(Type,n))

View(svData %>% filter(grepl('INF',ClusterDesc),ResolvedType=='PAIR_OTHER',grepl('Prox',ClusterReason),CNChgStart<0.25|CNChgEnd<0.25) %>% select(SampleId,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,NearestLen,CNStart,CNEnd,CNChgStart,CNChgEnd,Ploidy,everything()))

View(svData %>% group_by(Type,SampleId) %>% count %>% spread(Type,n))

View(svData %>% filter(ResolvedType=='PAIR_OTHER') %>% group_by(ClusterDesc,grepl('Prox',ClusterReason),Subclonal) %>% count() %>% spread(Subclonal,n))

View(svData %>% filter(ResolvedType=='PAIR_OTHER',ClusterDesc=='DUP=1_INF=1',!grepl('Prox',ClusterReason)))

View(svData %>% filter(Type=='SGL',grepl('LOH_',ClusterReason)) %>% group_by(low=round(ifelse(OrientStart==1,MinorAPStartPost,MinorAPStartPrev),1),high=round(ifelse(OrientStart==-1,MinorAPStartPost,MinorAPStartPrev),1)) %>% count() %>% spread(high,n))

View(svData %>% filter(Type=='',grepl('LOH_',ClusterReason)) %>% group_by(SampleId) %>% count())

View(svData %>% filter(Type=='INF',grepl('LOH_',ClusterReason)) %>% filter(MinorAPStartPost<0.5,MinorAPStartPrev<0.5) )


View(svData %>% filter(grepl('RECIP_TRANS',ResolvedType)) %>% group_by(ClusterDesc) %>% count())
View(svData %>% filter(grepl('SGL_PAIR',ResolvedType)) %>% filter(ClusterDesc=='SGL=2'))
###########################################
############# DRIVERS PER CLUSTER ##############
###########################################

View(svDrivers %>% filter(ClusterId>=0) %>% group_by(Gene,ClusterId,SampleId) %>% count() %>% group_by(ClusterId,SampleId) %>% summarise(count=n(),genes = paste0(Gene, collapse = "",sep=",")) %>% group_by(count) %>% count())

View(merge(svClusters %>% filter(ResolvedType=="COMPLEX") %>% mutate(ClusterCountBucket=2**(pmin(12,pmax(-3,round(log(ClusterCount,2),0))))),svDrivers %>% distinct(SampleId,ClusterId,Gene), by=c('SampleId','ClusterId'),all.x=T) %>% group_by(ClusterCountBucket,hasGene=!is.na(Gene)) %>% count() %>% spread(hasGene,n))

###########################################
############# HIGH LEVEL AMP ##############
###########################################

### CHARACTERISE ALL HIGH LEVEL AMPLIFICATIONS(>8)
highAmplifications = svData %>% group_by(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,CancerType,cancerSubtype) %>% summarise(
  FoldbackCount=sum(ifelse(IsFoldBack,1,0)),
  sumFoldbackPloidy=sum(ifelse(IsFoldBack,Ploidy,0)),
  maxFoldbackPloidy=max(ifelse(IsFoldBack,Ploidy,0)),
  INFCount=sum(ifelse(Type=="INF",1,0)),
  sumINFPloidy=sum(ifelse(Type=="INF",Ploidy,0)),
  maxINFPloidy=max(ifelse(Type=="INF",Ploidy,0)),
  SGLCount=sum(ifelse(Type=="SGL",1,0)),
  sumSGLPloidy=sum(ifelse(Type=="SGL",Ploidy,0)),
  maxSGLPloidy=max(ifelse(Type=="SGL",Ploidy,0)),
  MaxHighAmpPloidyMult=round(max(ifelse(Ploidy>7,pmax(Ploidy/(CNStart-CNChgStart-ifelse(OrientStart==1,MinorAPStartPost,MinorAPStartPrev)),
                  ifelse(Type=='SGL'|Type=='INF',-1,Ploidy/(CNEnd-CNChgEnd-ifelse(OrientEnd==1,MinorAPEndPost,MinorAPEndPrev)))),-1)),1),
  MaxHighAmpPloidyCount=sum(ifelse(ifelse(Ploidy>7,pmax(Ploidy/(CNStart-CNChgStart-ifelse(OrientStart==1,MinorAPStartPost,MinorAPStartPrev)),
                  ifelse(Type=='SGL'|Type=='INF',-1,Ploidy/(CNEnd-CNChgEnd-ifelse(OrientEnd==1,MinorAPEndPost,MinorAPEndPrev)))),-1)>2.3,1,0)),
  MaxCN=max(pmax(CNStart,CNEnd)),
  MaxPloidy=max(Ploidy)
    ) %>% filter(MaxPloidy > 8)

highAmplifications=merge(highAmplifications,beData %>% distinct(SampleId,ClusterId,LocTopId,LocTopType) %>% group_by(SampleId,ClusterId,LocTopType) %>% count() %>% spread(LocTopType,n), by=c('SampleId','ClusterId'))
highAmplifications=merge(highAmplifications,beData %>% distinct(SampleId,ClusterId,LocTopId,LocTopTI) %>% group_by(SampleId,ClusterId) %>%  summarise(shards=sum(LocTopTI)), by=c('SampleId','ClusterId'))
highAmplificationsUnmerged=(merge(highAmplifications,svDrivers %>% filter(grepl('GAIN',EventType)) %>% distinct(SampleId,ClusterId,Gene) %>% group_by(SampleId,ClusterId,Gene) %>% 
                                    summarise(countGenes=n()), by=c('SampleId','ClusterId'),all.x=T))
highAmplifications=(merge(highAmplifications,svDrivers %>% filter(grepl('GAIN',EventType)) %>% distinct(SampleId,ClusterId,Gene) %>% group_by(SampleId,ClusterId) %>% 
                           summarise(countGenes=n(),genes = paste0(Gene, collapse = "",sep=",")), by=c('SampleId','ClusterId'),all.x=T))



# Overall driver proportion
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(hasDriver = !is.na(genes)) %>% count)

# overall double minute proportion
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(hasDM=!is.na(DMSvTypes)) %>% count)

# frequency of DM by gene and cancer type
# generally ~15% DM across board
# Exceptions with higher DM rates:   CNS in general and particularly EGFR,  Colorectal particularly ERBB2,  MDM, MYC
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(CancerType,genes,hasDM=!is.na(DMSvTypes)) %>% count %>% spread(hasDM,n))
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(genes,hasDM=!is.na(DMSvTypes)) %>% count %>% spread(hasDM,n))
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(CancerType,hasDM=!is.na(DMSvTypes)) %>% count %>% spread(hasDM,n))

View(merge(highAmplificationsUnmerged,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(CancerType,Gene,hasDM=!is.na(DMSvTypes)) %>% count %>% spread(hasDM,n))

# Complex event very likely to have gene associations.   Simple DEL and DUP not.
View(merge(highAmplifications,svDoubleMinute, by=c('SampleId','ClusterId'),all.x=T) %>% group_by(ResolvedType.x,hasDriver = !is.na(genes)) %>% count %>% spread(hasDriver,n) )


###########################################
############# DEL & LOH ##############
###########################################

write.csv(svDrivers %>% filter(EventType=='LOH'|EventType=='DEL'),'~/temp.csv')

# HOM DELS are mostly simple DELs, but LOH are more likely complex.
View(svDrivers %>% filter(EventType=='LOH'|EventType=='DEL',Gene=='CDKN2A')  %>% group_by(EventType,ResolvedType) %>% count() %>% spread(EventType,n) %>% arrange(-LOH))

# TP53 & CDH1 have fewest focal LOH events
View(svDrivers %>% filter(grepl('LOH',EventType)|EventType=='DEL')  %>% group_by(EventType,Gene,Chromosome) %>% count() %>% spread(EventType,n,fill=0) %>% 
       mutate(LOHProp=round(LOH/(LOH+LOH_ARM+LOH_CHR),2)) %>% mutate(LOH_ARMProp=round(LOH_ARM/(LOH+LOH_ARM+LOH_CHR),2)) %>% 
       mutate(LOH_CHRProp=round(LOH_CHR/(LOH+LOH_ARM+LOH_CHR),2)))

# Prostate cancer stands out as having the most focal LOH events.   Ovarian the most whole chromosome LOH contributing to drivers
View(svDrivers %>% filter(grepl('LOH',EventType)|EventType=='DEL')  %>% group_by(EventType,CancerType) %>% count() %>% spread(EventType,n,fill=0) %>% 
       mutate(LOHProp=round(LOH/(LOH+LOH_ARM+LOH_CHR),2)) %>% mutate(LOH_ARMProp=round(LOH_ARM/(LOH+LOH_ARM+LOH_CHR),2)) %>% 
       mutate(LOH_CHRProp=round(LOH_CHR/(LOH+LOH_ARM+LOH_CHR),2)))

# simple DELs much more frequently involved in HOM_DELs than in HET_DELs compared to complex / arm level events
View(svDrivers %>% filter(grepl('LOH',EventType)|EventType=='DEL') %>% group_by(Gene,SampleId) %>% mutate(count=n()) %>% ungroup()  %>% group_by(count,EventType,ResolvedType) %>% count %>% spread(ResolvedType,n))

# Prostate cancer has relatively more 'COMPLEX' events involved in TSG drivers
View(svDrivers %>% filter(EventType=='DEL'|EventType=='LOH')  %>% group_by(CancerType,ResolvedType) %>% count %>% spread(ResolvedType,n) %>% arrange(-COMPLEX))

# Complex clusters with HOM_DELETIONS normally contain drivers
View(merge(svData %>% filter(CNStart-CNChgStart<0.5|(CNEnd-CNChgEnd<0.5&ChrEnd!=0)) %>% group_by(SampleId,ClusterId,ResolvedType) %>% count,
           svDrivers %>% distinct(SampleId,ClusterId,Gene),by=c('SampleId','ClusterId'),all.x=T) %>% group_by(ResolvedType,hasDriver=(!is.na(Gene))) %>% count %>% spread(hasDriver,nn))

View(svData %>% filter(SampleId=='COLO829T'))
 
###########################################
############# ViralInsertions #############
###########################################

View(svData %>% filter(VirusName !=''))# %>% group_by(SampleId,ChrStart,ClusterId) %>% count() %>% group_by(SampleId,ChrStart) %>% count())

View(svData %>% filter(VirusName !='') %>% select(VirusName,CancerType,everything()))

###########################################
############# SGLByRepeatClass #############
###########################################

View(svData %>% filter(Type=="SGL") %>% group_by(SampleId,ChrStart,RepeatClass,RepeatType,ResolvedType) %>% count() %>% spread(ResolvedType,n))

View(svData %>% filter(Type=="SGL",RepeatType %in% c('ALR/Alpha'),ResolvedType=='COMPLEX') %>% group_by(SampleId,ChrStart,ClusterId) %>% count() %>% group_by(SampleId,ChrStart) %>% count())

View(svData %>% filter(Type=="SGL",RepeatType %in% c('ALR/Alpha','HSATII')) %>% group_by(SampleId,ChrStart,ClusterCount,RepeatClass,RepeatType) %>% count() %>% spread(RepeatType,n))# %>% spread(ChrStart,n))
View(svData %>% filter(Type=="SGL",RepeatType %in% c('(GAATG)n','(CATTC)n')) %>% group_by(SampleId,ChrStart,ClusterCount,RepeatClass,RepeatType) %>% count() %>% spread(RepeatType,n))# %>% spread(ChrStart,n))

View(svData %>% filter(Type=="SGL",RepeatType %in% c('(GAATG)n','(CATTC)n')) %>% group_by(SampleId,ChrStart,ClusterCount,RepeatClass,RepeatType) %>% count() %>% spread(RepeatType,n))# %>% spread(ChrStart,n))

View(svData %>% filter(Type=="SGL",grepl('SAR',RepeatType)) %>% group_by(SampleId,ChrStart,RepeatClass,RepeatType) %>% count())# %>% spread(RepeatType,n))# %>% spread(ChrStart,n))

View(svData %>% filter(Type=="SGL",RepeatType=='L1PA2') %>% group_by(SampleId,ChrStart,ResolvedType,RepeatClass,RepeatType) %>% count() %>% spread(RepeatType,n))# %>% spread(ChrStart,n))


View(svData %>% filter(Type=="SGL") %>% group_by(RepeatClass,RepeatType,ResolvedType) %>% count()  %>% filter(n>500) %>% spread(ResolvedType,n))

View(svDrivers)

###################################
####### Complex Clusters ##########
###################################

View(svClusters %>% filter(ResolvedType=='COMPLEX') %>% group_by(Foldbacks,CC=pmin(20,ClusterCount),ArmCount) %>% count() %>% spread(CC,n))
merge(svData, sampleCancerTypes, by='SampleId', all.x=T)

View(merge(svClusters,sampleCancerTypes, by='SampleId', all.x=T) %>% filter(ResolvedType=='COMPLEX') %>% filter(ClusterCount<30,ArmCount-FragmentArms>=2))

View(beData %>% filter(ResolvedType=='COMPLEX',Chr!=0,ClusterCount<30) %>% group_by(SampleId,ClusterId,LocTopId,LocTopType,LocTopTI,Chr,Arm) %>% summarise(count=n(),minPos=min(Pos),maxPos=max(Pos),Range=maxPos-minPos) %>% mutate(nonTICount=count-LocTopTI*2) %>% arrange(-nonTICount))


View(beData %>% filter(ResolvedType=='COMPLEX',Chr!=0) %>% group_by(SampleId,ClusterId) %>% count())#  %>% group_by(SampleId,ClusterId,LocTopId,LocTopType,LocTopTI,Chr,Arm) %>% count %>% group_by(SampleId,ClusterId) %>% summarise(nn=n(),BE=sum(n)))# %>% group_by(nn) %>% summarise(n=n(),AvgBE=sum(BE)/n) )

View(beData %>% filter(ResolvedType=='COMPLEX',Chr!=0)  %>% group_by(SampleId,ClusterId,LocTopId,LocTopType,LocTopTI,Chr,Arm) %>% count %>% group_by(SampleId,ClusterId,Chr,Arm) %>% count() %>% group_by(SampleId,ClusterId) %>% count() %>% group_by(n) %>% count)
                       
head(beData)


#############################################
########## CLUSTER HISTORY ###################
#############################################
neoEpitopes = read.csv('~/hmf/analyses/SVAnalysis/SVA_NEO_EPITOPES.csv')
nrow(neoEpitopes)
View(neoEpitopes)

neoEpitopes = neoEpitopes %>% mutate(NovelAALength=stri_length(NovelAminoAcid))

View(neoEpitopes %>% filter(PhaseMatched=='true'&PhaseUp==PhaseDown))
View(neoEpitopes %>% filter(PhaseMatched=='false'))
View(neoEpitopes %>% filter(PhaseMatched=='true'&PhaseUp==PhaseDown&PhaseUp==0))
View(neoEpitopes %>% filter(grepl('_',UpstreamAminoAcids)))
View(neoEpitopes %>% filter(grepl('_',NovelAminoAcid)))
View(neoEpitopes %>% filter(grepl('_',DownstreamAminoAcids)))
View(neoEpitopes %>% group_by(PhaseMatched,PhaseUp,PhaseDown) %>% count())

View(neoEpitopes %>% separate(Fusion,c('UpGene','DownGene'),sep = "_") %>% mutate(sameGene=UpGene==DownGene,PhaseMatched=PhaseUp==PhaseDown)%>% group_by(sameGene,PhaseMatched) %>% count %>% spread(sameGene,n))

View(neoEpitopes %>% group_by(SampleId,SvIdUp,SvIdDown) %>% count() %>% filter(n>1) %>% arrange(-n))
View(neoEpitopes %>% group_by(SampleId,SvIdUp,SvIdDown) %>% count() %>% group_by(SampleId) %>% count)
#############################################
########## CLUSTER HISTORY ###################
#############################################

svClusterHistory = (read.csv(paste(PATH,'SVA_CLUSTERING_HISTORY.csv',sep=''), header = T, stringsAsFactors = F))
svClusterHistory = merge(svClusterHistory,svData %>% mutate(Ploidy=(PloidyMax+PloidyMin)/2) %>% select(SampleId,ChrStart,ChrEnd,PosStart,PosEnd,ClusterId,ClusterCount,ResolvedType,Ploidy,PloidyMin,PloidyMax,CNChgStart,CNChgEnd,Id),by.x=c('SampleId','SvId1'),by.y=c('SampleId','Id'),all.x=T)
svClusterHistory = merge(svClusterHistory,svData %>% mutate(Ploidy=(PloidyMax+PloidyMin)/2) %>% select(SampleId,ChrStart,ChrEnd,PosStart,PosEnd,Ploidy,PloidyMin,PloidyMax,CNChgStart,CNChgEnd,Id),by.x=c('SampleId','SvId2'),by.y=c('SampleId','Id'),all.x=T)
svClusterHistory=svClusterHistory %>% mutate(minCC=pmin(ClusterCount1,ClusterCount2))
head(svClusterHistory)

View(svClusterHistory %>% filter((PloidyMax.x<0.8&PloidyMax.x<PloidyMin.y&Ploidy.x+0.5<Ploidy.y)|
                                   (PloidyMax.y<0.8&PloidyMax.y<PloidyMin.x&Ploidy.y+0.5<Ploidy.x)) %>% 
       mutate(ploidyMax=round(pmin(PloidyMax.x,PloidyMax.y),1),CNMax=round(pmax(pmin(CNChgEnd.x,CNChgEnd.y),pmin(CNChgStart.x,CNChgStart.y)),1)) %>%
       group_by(Reason,ClonalDiscrepancy,ploidyMax,CNMax) %>% count %>% spread(CNMax,n))

View(svClusterHistory %>% filter(ClonalDiscrepancy=='true') %>% group_by(Reason) %>% count)

View(svClusterHistory %>% filter(SampleId == 'CPCT02230088T',ClusterId=='267'))
ggplot(data=svClusterHistory,aes(minCC)) + stat_ecdf(geom = "step", pad = FALSE) + ylim(0.8,1) + labs(title = 'CDF') + facet_wrap(~Reason) +scale_x_log10()
ggplot(data=svClusterHistory,aes(MinDistance)) + stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1)  + labs(title = 'CDF') + facet_wrap(~Reason) +scale_x_log10()
#############################################
########## INDIVIDUAL SAMPLE ANALYSIS #######
#############################################

chartSample = 'CPCT02080153T'; chartChr='7'; 
View(svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>% 
                                                    select(Length,ResolvedType,IsFoldBack,Ploidy,CNStart,CNChgStart,CNEnd,CNChgEnd,PloidyMin,PloidyMax,Type,ClusterCount,everything()) %>% arrange(-Ploidy))
View(svData %>% filter(Type=='SGL',ResolvedType=='COMPLEX') %>% group_by(substring(ClusterReason,1,4)) %>% count)
chartSample = 'CPCT02380011T'; View(svData %>% filter(SampleId==chartSample,ClusterCount>20,) %>% 
                                      select(ResolvedType,IsFoldBack,ChainId,ChainIndex,ChainCount,PloidyMin,PloidyMax,Ploidy,CNStart,CNChgStart,CNEnd,CNChgEnd,Type,ClusterCount,everything()) %>% arrange(-Ploidy))
View(svLinks %>% filter(SampleId=='CPCT02010963TII',ClusterId==524) %>% select(PosStart,PosEnd,everything()))
View(svData %>% filter(SampleId=='CPCT02010898T',ClusterId==87) %>% select(IsFoldBack,Ploidy,PloidyMin,PloidyMax,everything()))
View(beData %>% filter(SampleId=='CPCT02100090T',ClusterId==77))
View(svDrivers %>% filter(SampleId=='CPCT02010569T'))


1#1. Whole Chr View
View(svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>%select(ResolvedType,IsFoldBack,Ploidy,PloidyMin,PloidyMax,CNChgStart,CNChgEnd,Type,ClusterCount,ClusterReason,everything()) %>% arrange(-PloidyMin))
print(ggplot(data = svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>% mutate(
  modEnd=ifelse(ChrEnd==chartChr,ifelse(ChrStart==chartChr,PosEnd,as.numeric(ChrStart)*-1e7),as.numeric(ChrEnd)*-1e7),modStart=ifelse(ChrStart==chartChr,PosStart,PosEnd)),aes(modStart,modEnd)) + 
    geom_point(aes(size = Ploidy,colour=ifelse(ResolvedType == 'LINE'|ClusterCount==2,ResolvedType,ClusterDesc))) + theme_bw())

#1a. all chromosomes
print(ggplot(data = rbind(svData %>% filter(Subclonal=='false',SampleId==chartSample) %>% mutate(modStart=PosStart*1e7,modEnd=ifelse(ChrEnd==ChrStart,PosEnd,as.numeric(ChrEnd)*-1e7)) %>% select(ChrStart,Ploidy,ResolvedType,ClusterCount,ClusterDesc,modStart,modEnd),
                          svData %>% filter(Subclonal=='false',ChrEnd!=ChrStart,SampleId==chartSample,ChrEnd!=0) %>% mutate(modStart=PosEnd*1e7,modEnd=as.numeric(ChrStart)*-1e7) %>% select(ChrStart,Ploidy,ResolvedType,ClusterCount,ClusterDesc,modStart,modEnd)),
             aes(modStart,modEnd)) + geom_point(aes(size = Ploidy,colour=ifelse(ResolvedType == 'LINE'|ClusterCount==2,ResolvedType,ClusterDesc))) + theme_bw() +facet_wrap(~ChrStart))




#2.By CHR for large clusters
View(svData %>% filter(SampleId == chartSample,ClusterCount>2,ResolvedType!='LINE') %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% tally() %>% spread(chr_arm_end,n,fill=""))

View(svData %>% filter(SampleId == chartSample,IsFoldBack) %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% tally() %>% spread(chr_arm_end,n,fill=""))

# 3. Single Cluster View
View(svData %>% filter(SampleId==chartSample,ClusterId==60)%>% arrange(ChrStart,PosStart) 
     %>% select(ResolvedType,LengthBucket,Ploidy,PloidyMin,PloidyMax,IsFoldBack,Type,AsmbStart,AsmbEnd,ClusterReason,everything()))

# 4. Large cluster View
View(svData %>% filter(!IsLowQual,SampleId==chartSample,ClusterCount>=1,ClusterCount>0,ResolvedType!='ASimpleSV')%>% 
       select(ResolvedType,Length,Ploidy,Type,AsmbStart,AsmbEnd,everything()) %>% arrange(-Length))#%>% group_by(ClusterCountBucket,CnChStartBucket) %>% tally() %>% spread(CnChStartBucket,n))


View(svData %>% filter(SampleId==chartSample,ResolvedType!='LowQual') %>% 
       select(IsFoldBack,ResolvedType,Length,ClusterId,Ploidy,Type,ClusterCount,ChrStart,ChrEnd,PosStart,PosEnd,ClusterReason,everything()) %>% arrange(-Ploidy))



###### FILE GENERATION
svDrivers %>% filter(Gene=='MAP2K4',ClusterId!='-1') %>% mutate(temp=paste('/do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" "))# %>% select(temp)


View(svClusterHistory %>% filter(ClonalDiscrepancy=='true'))
write.csv(svClusterHistory %>% filter(ClonalDiscrepancy=='true',Reason=='LOH',ClusterCount<20) %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/clonalDiscrepancy.txt',row.names = FALSE, quote=FALSE)

write.csv(svClusters %>% filter(ResolvedType=='LINE',ClusterCount>30,ClusterCount<80) %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/complexLINE.txt',row.names = FALSE, quote=FALSE)

View(beData %>% filter(ResolvedType=='LINE',ClusterCount>20,ClusterCount>40,Chr!=0) %>% group_by(Chr,SampleId,ClusterId) %>% count() %>% filter(n>2) %>% group_by(SampleId,ClusterId) %>% count())
View(svData %>% filter(ClusterId=='286',SampleId=='CPCT02020501T'))
write.csv(svDrivers %>% filter(Gene=='SMAD4',ClusterId!='-1') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/SMAD4.txt',row.names = FALSE, quote=FALSE)
write.csv(svDrivers %>% filter(Gene=='PTEN',CancerType=='Prostate',ClusterId!='-1') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/PTEN.txt',row.names = FALSE, quote=FALSE)

write.csv(svDrivers %>% filter(Gene=='AR',ClusterId!='-1') %>% distinct(SampleId,ClusterId)  %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/AR.txt',row.names = FALSE, quote=FALSE)
write.csv(svDrivers %>% filter(Gene=='EGFR',ClusterId!='-1') %>% distinct(SampleId,ClusterId)  %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/EGFR.txt',row.names = FALSE, quote=FALSE)

write.csv(svDrivers %>% filter(Gene=='PTEN',CancerType=='Prostate',ClusterId!='-1') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/PTEN.txt',row.names = FALSE, quote=FALSE)


write.csv(svDrivers %>% filter(Gene=='APC',ClusterId!='-1',DriverType!='MUTATION') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/APC.txt',row.names = FALSE, quote=FALSE)

write.csv(svDrivers %>% filter(Gene=='KRAS',ClusterId!='-1',DriverType!='MUTATION') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/KRAS.txt',row.names = FALSE, quote=FALSE)


write.csv(svDrivers %>% filter(Gene=='MDM2',ClusterId!='-1') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/MDM2.txt',row.names = FALSE, quote=FALSE)

write.csv(svFusions %>% filter(KnownType=='Known')  %>% filter(RegionTypeDown=='Exonic'|BreakendExonDown!=FusedExonDown|BreakendExonUp!=FusedExonUp|ExonsSkippedDown>0|ExonsSkippedUp>0)  %>% 
            distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/UnusualFusions.txt',row.names = FALSE, quote=FALSE)

write.csv(svFusions %>% filter(KnownType=='Known') %>% 
            distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")) %>% select(temp),'~/KnownFusions.txt',row.names = FALSE, quote=FALSE)

View()
write.csv(head(svDrivers %>% filter(grepl('RECIP',ResolvedType),ClusterId!='-1') %>% distinct(SampleId,ClusterId) %>% mutate(temp=paste('./do_run_linx_vis ',SampleId,"-clusterId",ClusterId,sep=" ")),100) %>% select(temp),'~/RECIP.txt',row.names = FALSE, quote=FALSE)


View(svDrivers %>% filter(Gene=='SMAD4',ClusterId!='-1') )

##########################

View(svClusters %>% group_by(SampleId,ResolvedType) %>% summarise(clusters=n()) %>% spread(ResolvedType,clusters))
View(beData  %>% group_by(SampleId,ResolvedType,ClusterId,Chr) %>% count %>% group_by(SampleId,ResolvedType,Chr) %>% summarise(BE=sum(n),clusters=n()) %>% spread(chr,clusters))

View(beData  %>% filter(Chr>0) %>% group_by(SampleId,ResolvedType,ClusterId,Chr) %>% count %>% group_by(SampleId,ResolvedType,Chr) %>% summarise(BE=sum(n),clusters=n()))

View(svClusters %>% filter(SampleId=='DRUP01050035T') %>% group_by(SampleId,ResolvedType,ClusterDesc) %>% summarise(variants=sum(ClusterCount),clusters=n()))

View(svData %>% group_by(NearestLenBucket))


#############################

###### DEALING WITH INF ############

#1. Inconsistent CNChg
# Find all non SGL, non INF variants where ploidy and CN change on one end are consistent, but the CNChg on the other end is significantly lower ( perhaps use higher thresholds? - 0.6, 20%)
# TEST IF A NEIGHBOURING INF can explain the CN discrepancy.   If so:
  # Filter the INF
  # Mark the CNChg on that end as unreliable (ie. don't use in ploidy & ploidy uncertainty calc and don't use in 'available chromatid ploidy' calculations)
# IF NOT NEIGBOURING INF & Breakend is isolated & NOT a foldback  then:
  # Create an offsetting INF with opposite orientation ploidy = CN discrepancy,ploid uncertainty = ..
  # Mark the CNChg on that end as unreliable  (ie. don't use in ploidy & ploidy uncertainty calc and don't use in 'available chromatid ploidy' calculations)
# Rough query for CNStart
View(svData %>% filter(Type!='INF',Type!='SGL',abs(CNChgStart-Ploidy)<0.5|abs(CNChgStart-Ploidy)/CNStart<0.15,!(CNChgStart>0.5&Ploidy>0.5),(pmin(Ploidy,CNChgStart)-CNChgEnd)/CNEnd>0.2,pmin(Ploidy,CNChgStart)-CNChgEnd>0.6) %>% 
       filter(IsFoldBack,NearestLen!=0,NearestLen>5000|(ClusterCount==2&grepl('INF',ClusterDesc))) %>%
       select(Type,ResolvedType,CNChgStart,CNChgEnd,Ploidy,PloidyMin,PloidyMax,CNStart,CNEnd,NearestLen,NearestType,everything())
     )

#2. 

View(svData %>% filter(ResolvedType=='PAIR_OTHER',ClusterDesc=='DEL=1_INF=1',NearestLen<5000) %>%  select(Type,ResolvedType,CNChgStart,CNChgEnd,Ploidy,PloidyMin,PloidyMax,CNStart,CNEnd,NearestLen,NearestType,everything()))

View(svData %>% filter(grepl('SGL',ClusterDesc),ResolvedType!='COMPLEX',ResolvedType!='LINE')%>% group_by(ResolvedType,ClusterDesc) %>% count)

View(svData %>% filter(Type=='INF')%>% group_by(ResolvedType,VAF=pmin(floor(CNChgStart/CNStart*20)/20,1)) %>% count %>% spread(VAF,n))



write.csv(svData %>% filter(ResolvedType=='LINE') %>% group_by(SampleId) %>% count(),'~/temp.csv')

View(svDrivers %>% filter(EventType=='GAIN',ClusterId==-1))



#1.Foldback length distribution for simple inversions short chained foldbacks
print(ggplot(data = foldbacks %>% filter(Subclonal=='false',ClusterDesc=='INV=2',ResolvedType=='RECIP_INV_DEL_DUP') %>% group_by(LengthBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
             aes(x=LengthBucket, y=Count))
      + geom_line(aes(y=INV, colour='Simple'))  #+ geom_line(aes(y=Combo, colour='Synthetic' )) 
      + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))

print(ggplot(data = foldbacks %>% filter(Subclonal=='false',ClusterDesc=='INV=2',ResolvedType=='RECIP_INV_DEL_DUP') %>% group_by(LengthBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
             aes(x=LengthBucket, y=Count))
      + geom_line(aes(y=INV, colour='Simple'))  #+ geom_line(aes(y=Combo, colour='Synthetic' )) 
      + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))


View(svData %>% filter(Subclonal=='false',ClusterDesc=='INV=2',ResolvedType=='RECIP_INV_DEL_DUP') %>% select(,everything()))


View(svData %>% filter(ClusterDesc=='INV=2',PosEnd-PosStart>1e7))



######################################################
############# PAIR length characteristics  ############
######################################################
pairClusters = svClusters %>% filter(ClusterCount==2,(ResolvedType=='DEL'&Synthetic=='true')|(ResolvedType=='DUP'&Synthetic=='true')|
                                   ResolvedType=='DUP_TI'|ResolvedType=='DEL_TI'|ResolvedType=='RECIP_TRANS'|ResolvedType=='RECIP_INV'
                                   |ResolvedType=='RECIP_TRANS_DUPS'|ResolvedType=='RECIP_TRANS_DEL_DUP'
                                   |ResolvedType=='RECIP_INV_DUPS'|ResolvedType=='RECIP_INV_DEL_DUP')

pairClusters = pairClusters %>% separate('Annotations',c('SynLength1','SynLength2','SynGapLength'),sep=';') %>% 
  mutate(SynGapLength=as.numeric(as.character(SynGapLength)),SynLength1=as.numeric(as.character(SynLength1)),SynLength2=as.numeric(as.character(SynLength2)))

ggplot(data=pairClusters %>% filter(Subclonal=='false'),aes(SynLength1)) + stat_ecdf(geom = "step", pad = FALSE) + labs(title = 'CDF') + facet_wrap(~ResolvedType) + scale_x_log10()
ggplot(data=pairClusters %>% filter(Subclonal=='false'),aes(SynLength2)) + stat_ecdf(geom = "step", pad = FALSE) + labs(title = 'CDF') + facet_wrap(~ResolvedType) + scale_x_log10()
ggplot(data=pairClusters %>% filter(Subclonal=='false'),aes(SynGapLength)) + stat_ecdf(geom = "step", pad = FALSE) + labs(title = 'CDF') + facet_wrap(~ResolvedType) + scale_x_log10()

ggplot(data=pairClusters %>% filter(Subclonal=='false'),aes(SynLength1,SynLength2)) + geom_point() + facet_wrap(~ResolvedType) + scale_x_log10() + scale_y_log10()

ggplot(data=pairClusters %>% filter(Subclonal=='false') %>% filter(Synthetic!='Afalse') %>% unite(temp,ResolvedType,ClusterDesc),aes(SynLength1,SynLength2)) + geom_hex(bins=30) + facet_wrap(~temp) + scale_x_log10() + scale_y_log10()+
  scale_fill_gradient(low = 'light blue', high = 'dark blue')
#############################


View(svFusions %>% filter(KnownType=='Known') %>% select(CodingTypeUp,CodingTypeDown,ExonsSkippedUp,ExonsSkippedDown,FusedExonUp,FusedExonDown,everything()))


View(svClusters %>% filter(ResolvedType=='COMPLEX'))

View(svData %>% filter(ResolvedType!='LINE',ClusterCount<=3) %>% group_by(IsPolyA,ResolvedType,ClusterDesc) %>% count() %>% spread(IsPolyA,n))

head(svDrivers)
View(svDrivers %>% group_by(EventType) %>% count)

View(svDoubleMinute %>% filter(SampleId=='CPCT02040228T'))


View(svData %>% filter(ResolvedType=='LINE',ResolvedType!='DUP_BE',ClusterCount<4) %>% 
       group_by(ClusterDesc,IsPolyA) %>% count() %>% spread(IsPolyA,n))

View(svData %>% filter(ResolvedType=='SGL_PAIR_INS'|ResolvedType=='LINE') %>% group_by(ResolvedType,SampleId) %>% count %>% spread(ResolvedType,n))

View(svData %>% filter(ClusterDesc=='BND=1_INF=1',ClusterCount) %>% 
       group_by(ResolvedType,IsPolyA) %>% count() %>% spread(IsPolyA,n))

View(svFusions %>% filter(SampleId=='WIDE01010089T'))

View(svData %>% filter(Subclonal=='true') %>% group_by(SampleId,ClusterId) %>% mutate(maxClusterPloidy=max(Ploidy)) %>% ungroup() %>% filter(maxClusterPloidy>2) %>% 
       group_by(SampleId,ClusterId,ResolvedType,PloidyBucket) %>% count() %>% spread(PloidyBucket,n))

View(svData %>% group_by(SampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n))

View(svData %>% filter(ClusterCount<5) %>% group_by(ResolvedType,ClusterDesc,IsPolyA) %>% count() %>% spread(IsPolyA,n))

View(svData %>% filter(ResolvedType!='LINE',ClusterCount<50,IsPolyA) %>% group_by(ResolvedType,ClusterDesc,SampleId,ClusterId) %>% count() %>% filter(n>1))
BND=1_INF=1_SGL=1

View(svData %>% filter(ClusterDesc=='BND=1_SGL=1',ResolvedType=='LINE') %>% select(IsPolyA,LEStart,LEEnd,everything()))

View(svData %>% filter(Subclonal=='false',ResolvedType=='UNBAL_TRANS_TI') %>% group_by(ClusterCount) %>% tally() )
View(svData %>% filter(Subclonal=='false',ResolvedType=='DUP_TI') %>% group_by(ClusterCount) %>% tally() )

write.csv(svDrivers %>% filter(Gene=='TP53'),'~/temp.csv')

View(svDrivers)# %>% filter(Gene=='TP53') %>% group_by(SampleId) %>% count)
View(svData %>% group_by(SampleId) %>% count)

View(svDrivers)
View(svClusters %>% filter(ClusterCount>10,ClusterCount<60,AssemblyTIs>ClusterCount-5,BndCount>ClusterCount/2))

View(svLinks %>% filter(LinkReason=="ASSEMBLY") %>% group_by(ResolvedType,ClusterDesc,ClusterCount,Id2,PosEnd) %>% count() %>% filter(n>1))
