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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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

scatterPlot<-function(data,xVar,YVar,useLogX = TRUE,useLogY = TRUE) {
  plot<-ggplot(data=data,aes_string(xVar,YVar))+geom_point() 
  if (useLogX == TRUE) {
    plot<-plot+scale_x_log10()
  }
  if (useLogY == TRUE) {
    plot<-plot+scale_y_log10()
  }
  print(plot)
}

sv_set_common_fields<-function(cluster){ cluster %>% mutate( 
    IsLINE = ifelse(LEStart!='None'|LEEnd!='None',T,F),
    IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F),
    IsGenicStart = ifelse(GeneStart!='',T,F),
    IsGenicEnd = ifelse(GeneEnd!='',T,F),
    Length = ifelse(as.character(ChrStart)!=as.character(ChrEnd)|Type=='INS'|ArmEnd!=ArmStart, -1, PosEnd-PosStart),
    DoubleDupBE = ifelse(DupBEStart=='true'&DupBEEnd=='true',T,F),
    SingleDupBE = ifelse(DoubleDupBE==0&(DupBEStart=='true'|DupBEEnd=='true'),T,F),
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
    IsLowQual = ResolvedType=='LowQual'
    
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
    LengthBucket=ifelse(Type=='BND'|Type=='INS'|PosEnd-PosStart==0|ArmEnd!=ArmStart,0,2**(pmin(25,pmax(5,round(log(PosEnd-PosStart,2),0))))),
    HomLenBucket=2**(round(log(nchar(Homology),2),0)),
    NearestLenBucket=2**(round(log((NearestLen),2),0)),
    SynDelDupLenBucket=4**(round(log((SynDelDupLen),4),0)),
    SynDelDupTILenBucket=4**(round(log((SynDelDupTILen),4),0)),
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
              countRecipTrans=sum(ResolvedType=='RecipTrans'),
              #countShortTI=sum(ShortTICount),
              countFS=sum(countFS))
}

#############################################
################ LOADING ####################
#############################################
PATH='~/Dropbox/HMF Australia team folder/Structural Variant Analysis/'
svData = read.csv(paste(PATH,'CLUSTER.csv',sep=''), header = T, stringsAsFactors = F)
svClusters = (read.csv(paste(PATH,'SVA_CLUSTERS.csv',sep=''), header = T, stringsAsFactors = F))
svLinks = (read.csv(paste(PATH,'SVA_LINKS.csv',sep=''), header = T, stringsAsFactors = F))
sampleCancerTypes= (read.csv(paste(PATH,'sample_cancer_types.csv',sep=''), header = T, stringsAsFactors = F))
svChordStatus = read.csv(paste(PATH,'SVChordStatus.csv',sep=''), header = T, stringsAsFactors = F)
svDrivers = read.csv(paste(PATH,'SVDrivers.csv',sep=''), header = T, stringsAsFactors = F)
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
sampleList = svData %>% distinct(SampleId,CancerType)
#Annotations
svData=createBuckets(svData)
svData=sv_set_common_fields(svData)
svSampleSummary=createSampleSummary(svData)
foldbacks=createFoldbacks(svData)
dbData = rbind(svData %>% filter(DBLenStart>-31) %>% mutate(DBLength = DBLenStart,Assembled = ifelse(AsmbMatchStart=="MATCH","Assembled","NotAssembled"),LE=LEStart,RefContext=RefContextStart,Orient=OrientStart),
               svData %>% filter(DBLenEnd>-31) %>% mutate(DBLength = DBLenEnd, Assembled = ifelse(AsmbMatchEnd=="MATCH","Assembled","NotAssembled"),LE=LEEnd,RefContext=RefContextEnd,Orient=OrientEnd)) %>%
               mutate(DBLenBucket = ifelse(DBLength==0,0,ifelse(DBLength<0,-(2**round(log(-DBLength,2))),2**round(log(DBLength,2)))))


#############################################
########## OVERVIEW OF COHORT ###############
#############################################

# 1. By RESOLOVED Type
View(svData %>% filter(!IsLowQual) %>% group_by(ResolvedType,Type) %>% count() %>% spread(Type,n))

# 2. By CLUSTERCOUNTBUCKET
View(svData %>% filter(!IsLowQual) %>% group_by(ClusterCountBucket,Type) %>% count() %>% spread(Type,n))


##############################################################
########## Short Deletion Bridge Analysis ####################
#############################################################

#1.DBLength by ResolvedType => less than 50 bases (NB - DBLength of 1 means exact break - shoould correct this)
# Note the LINE double peak and the sharp feature for Reciprocal Inversion, Reciprocal Translocations and None
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual",'DBLength,ResolvedType'),'DBLength','ResolvedType','DBLengthByResolvedType(<50bases)',useLogX = F,useLogY = F)

plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual",'DBLength,IsLine=ResolvedType=="Line"'),'DBLength','IsLine','DBLengthByIsLine(<50bases)',useLogX = F,useLogY = F)

#2. DB Length - Exact break short INV feature also exists in chains but only when INV length is <1e3   
# only exists for Short INV
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,ResolvedType %in% c('SimpleChain','ComplexChain')",'DBLength,LengthBucket=Length<0|Length>1e3'),'DBLength','LengthBucket','DBLength for Complex Chains short or long variant',useLogX = F,useLogY = F)
# Occurs prinicipally in the most complex clusters.  Short Inversions are much more likely than other variant types than other short variants.
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length<1e3,ResolvedType %in% c('SimpleChain','ComplexChain')",'DBLength,ClusterCountBucket'),'DBLength','ClusterCountBucket','DBLength for short variants in chains by cluster count',useLogX = F,useLogY = F)
# Is not primarily a foldback feature
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length<1e3,ResolvedType %in% c('SimpleChain','ComplexChain')",'DBLength,IsFoldBack'),'DBLength','IsFoldBack','DBLength for short variants in chains by isFoldbacl',useLogX = F,useLogY = F)
# Mostly around 60-500 bases.
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,Length<1e3,ResolvedType %in% c('SimpleChain','ComplexChain')",'DBLength,LengthBucket'),'DBLength','LengthBucket','DBLength for short variants in chains by length bucket',useLogX = F,useLogY = F)

#3. 'NONE' clusters - appear to have characteristics of short INVs and likely LINE elements for BND + SGL
plot_count_by_bucket_and_type(cohortSummary(dbData,"DBLength<=50,!IsLowQual,ResolvedType %in% c('None')",'DBLength,LengthBucket=Length<0|Length<1e3'),'DBLength','LengthBucket','DBLength for Complex Chains ByLengthBucket>1e3(<50bases)',useLogX = F,useLogY = F)

#4. LINE Elements - no clear indication of why there are 2 peaks.
#The 2 observed DB length peaks for LINE elements are NOT sample or cancer type specific.   Also appars to be the same regardless of source element
print(ggplot(data = dbData %>% filter(DBLength<=50,ResolvedType %in% c('Line','SglPair_INS')) %>% group_by(CancerType,OLPeak=DBLength<(-7)) %>% count(), 
             aes(x = reorder(CancerType, -n), y = n, fill =OLPeak))
      + geom_bar(stat = "identity", colour = "black") + ylab("DB Count") + xlab("Tumor Type") + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15)) + labs(title = "DB Counts for LINE Elements by cancerType and Overlapping Peak"))
# Also does not seem to depend significantly on refContext.  
View(dbData %>% filter(ResolvedType=='Line',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context),longOverlap=DBLength<(-7) %>% count() %>% arrange(-n) %>% spread(longOverlap,n))

#######################################################
########## LINE specific properties ###################
#######################################################

# 1. INSERTION MOTIF: ~15% are A-TTTTT with another 30% or so a similar variant
View(dbData %>% filter(ResolvedType=='Line',LE=='None',IsPolyA,RefContext!='') %>% mutate(context=stri_reverse(ifelse(Orient==-1,base_complements(substring(RefContext,10,15)),(substring(RefContext,7,12))))) %>% 
       group_by(context) %>% count() %>% arrange(-n))

# 2. CLUSTER COUNT vs TOTAL COUNT Count of clusters per sample is very related to the total number of variants (averaging at 3-4 per cluster)
scatterPlot(svClusters %>% filter(ResolvedType=='Line') %>% group_by(SampleId) %>% summarise(countLine=n(),sumLine=sum(ClusterCount)),'countLine','sumLine',F,F)

# 3. LINE Clusters affect a predictable number of chromosome arms that grows with the svCount in the cluster
scatterPlot(svClusters %>% filter(ResolvedType=='Line'),'ClusterCount','ArmCount',F,F)

# 4. LINE elements mmuch more common in Esophagus and Stomach cancers
ggplot(data=merge(sampleList,svClusters %>% filter(ResolvedType=='Line',!(ClusterDesc %in% c('SGL','SGL=2','INS'))) %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of LINE Clusters per Sample') + facet_wrap(~CancerType)

# 5. Known and suspected LINE element source locations => many potentially novel source locations that feature
View(rbind(svData %>% filter(ResolvedType=='Line') %>% filter(LEStart!='None') %>% group_by(chr=ChrStart,pos=round(PosStart,-4),SampleId) %>% summarise(n=n(),countSuspected=sum(LEStart=='Suspect'),
           countKnown=sum(LEStart=='Known'),countPOLYA=sum(IsPolyA&Type=='BND'),countSuspectKnown=sum(LEStart=='Known;Suspect')),
           svData %>% filter(ResolvedType=='Line') %>% filter(LEEnd!='None') %>% group_by(chr=ChrEnd,pos=round(PosEnd,-4),SampleId) %>% summarise(n=n(),countSuspected=sum(LEEnd=='Suspect'),
           countKnown=sum(LEEnd=='Known'),countPOLYA=sum(IsPolyA&Type=='BND'),countSuspectKnown=sum(LEEnd=='Known;Suspect'))) %>%
       group_by(chr,pos,SampleId) %>% summarise(n=sum(n),countSuspected=sum(countSuspected),countSuspectKnown=sum(countSuspectKnown),countKnown=sum(countKnown),countPOLYA=sum(countPOLYA)) %>% 
       group_by(chr,pos) %>% summarise(n=sum(n),countSuspected=sum(countSuspected),countSuspectKnown=sum(countSuspectKnown),countKnown=sum(countKnown),countPOLYA=sum(countPOLYA),numSamples=n()) )

################################################
############## FOLDBACKS #######################
################################################
foldbackClusters=(foldbacks %>% filter(!IsLowQual) %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(FBcount=n()/2,
      ClusterCount=first(ClusterCount),CNMax=max(pmax(AdjCNStart,AdjCNEnd)),FBPloidyMax=max(pmax(AdjCNChgStart,AdjCNChgEnd,Ploidy))))

# TO DO: FOLDBACKS - combos > 5kb length should be just excluded by charles

#1.Foldback length distribution for simple inversions short chained foldbacks
print(ggplot(data = foldbacks %>% filter(FoldbackType=="INV"|FoldbackChainLengthBucket<5000) %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=INV, colour='INV')) + geom_line(aes(y=Combo, colour='Combo ( FB Chain Lenght < 5k)' )) + theme_bw()
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + labs(title = "Foldback Length Distribution"))

#2. Longer chains are more prone to mis-chaining and do not have a similar distribution.
print(ggplot(data = foldbacks %>% filter(FoldbackType=="Combo") %>% group_by(FoldbackLenBucket,LongChain=FoldbackChainLengthBucket>5000) %>% summarise(Count=n()) %>% spread(LongChain,Count), 
      aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=`TRUE`, colour='Chain Length > 5k')) + geom_line(aes(y=`FALSE`, colour='Chain Length < 5k'))
      + scale_x_log10() + theme_bw() + labs(title = "Combo Foldback Length Distribution by Chain Length"))      

#3. Count of Foldbacks per cluster
# TO Do - Why do so many samples have 5+ separate clusters with foldbacks???
ggplot(data=merge(sampleList,foldbackClusters %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of Foldback Clusters per Sample')
View(merge(sampleLikst,foldbackClusters %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T) %>% replace_na(list(n=0)))
View(foldbackClusters %>% group_by(SampleId,CC=pmin(5,ClusterCount)) %>% count() %>% arrange(-n) %>% spread(CC,n)) 

#4. Combo foldbacks counts increase with simple foldback counts
scatterPlot(foldbacks%>% group_by(SampleId) %>% summarise(simpleFoldbackCount=sum(FoldbackType=='INV'),comboFoldbackCount=sum(FoldbackType!='INV'&FoldbackChainLengthBucket<5000)),
            'simpleFoldbackCount','comboFoldbackCount',F,F)

############################################################
########### Reciprocal Inversions ##########################
############################################################

# #1. TO DO:  show the length distribution and come up with a standardised definition
# eg.
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType=='DEL_Int_TI',Type=='INV',!IsLowQual,SynDelDupTILen>=0.9*SynDelDupLen",'CancerType,LengthBucket') 
                              ,'LengthBucket','CancerType','Reciprical INV',useLogY =F )

#2. Reciprocal INV NOT enriched in typical DEL lengths OR FS!
plot_count_by_bucket_and_type(cohortSummary(svData,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsFS",'CancerType,LengthBucket'),'LengthBucket','CancerType','AllDel&Dup',useLogY =T )


############################################################
########### Reciprocal Translocations ######################
############################################################

# TO DO:  show the length distribution and come up with a standardised definition

############################################################
########### Short TIs #####################################
############################################################

#1. Short TIs present in nearly all samples at a rate of ....
# TO DO: show this

#2. Short TIs create 'synthetic' events that match the profile
scatterPlot(svSampleSummary,'countSimpleDEL','countDEL_Ext_TI',F,F) # DEL correlated with synthetic DEL
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_Ext_TI',F,F)  # DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUP','countDUP_Int_TI',F,F)  # Int DUP correlated with synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPShort','countDUP_Ext_TIShort',F,F)  # Short DUP correlated with short synthetic DUP
scatterPlot(svSampleSummary,'countSimpleDUPLong','countDUP_Ext_TILong',F,F)  # Long DUP correlated with long synthetic DUP

# TO DO: Show the same for DEL_Int_TI after removing reciprocal inversions
scatterPlot(svSampleSummary,'countSimpleDEL','countDEL_Int_TI') # this NEEDS MORE THOUGHT

##############################################################
########## FRAGILE SITE ######################################
##############################################################

#1. Fragile Site Length Distribution
plot_count_by_bucket_and_type(cohortSummary(svData,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV')","LengthBucket,IsFS"),'LengthBucket','IsFS','Length Distribution by isFS',useLogY =F)

#2. FS like length dels outside FS are correlated, but not all 20k-500k are FS Dels
scatterPlot(svData %>% filter(Type=='DEL'&ResolvedType=='SimpleSV') %>% group_by(SampleId) %>% summarise(countFS=sum(IsFS),countFSLike=sum(!IsFS&PosEnd-PosStart>2e4&PosEnd-PosStart<5e5)),'countFS','countFSLike',F,F) 

#3. Extreme tail possible for FS length DELS => Unrelated process
ggplot(data=svData %>% filter(ResolvedType=='SimpleSV',PosEnd-PosStart>20000,PosEnd-PosStart<500000,Type=='DEL') %>% group_by(SampleId,IsFS) %>% count() %>% spread(IsFS,n,fill=0)) + 
  stat_ecdf(aes(`TRUE`,color='isFS'),geom = "step", pad = FALSE) + stat_ecdf(aes(`FALSE`,color='FSLength non-FS'),geom = "step", pad = FALSE) + 
  ylim(0,1) + labs(title = 'CDF # of FS length DELs')

#4. Location of FS + potentially novel FS
View(svData %>% filter(ResolvedType=='SimpleSV',Type=='DEL'|Type=='DUP',PosEnd-PosStart>2e4,PosEnd-PosStart<5e5) %>% group_by(chr=ChrStart,pos=round(PosStart,-6), IsFS) %>% summarise(n=n()) %>% spread(IsFS,n,fill=0) %>% mutate(Total=`TRUE`+`FALSE`))

#5. Types of events in FS
# TO DO: split out DEL_INT_TI into RecipInv=T/F
View(svData %>% group_by(IsFS,ResolvedType,Type) %>% count() %>% spread(IsFS,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))

#6. TO DO: CHECK GENIC vs non-GENIC

#################################################
################# COMPLEX CLUSTER STATS #########
#################################################

#1 LOCAL - REMOTE CONSISTENCY
clusterArmStats = createClusterArmStats(svData)
View(clusterArmStats %>% filter(LocalCount+RemoteCount+SGLCount>10))

# 2. Per arm by ClusterCount
View(rbind(svData %>% unite(chr_arm,ChrStart,ArmStart) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='Line') %>% group_by(SampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count(),
           svData %>%  unite(chr_arm,ChrEnd,ArmEnd) %>% filter(!IsLowQual,Type!='NONE',ResolvedType!='Line') %>%group_by(SampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count()) %>% 
       group_by(SampleId,chr_arm,CC,id) %>% count() %>%
       group_by(SampleId,CC,chr_arm) %>% count() %>% spread(CC,n) %>%
       filter(chr_arm!='0_P'))

#3. LARGE Clusters - # arms not super correlated with # of SVs
scatterPlot(svClusters %>% filter(ClusterCount>5,ResolvedType!='LINE') %>% mutate(realArmCount=ArmCount-FragmentArms),'realArmCount','ClusterCount',F,F)

#4. CDF of large non-line clusters per sample
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>10,ResolvedType!='LINE') %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>10 Per Sample')
ggplot(data=merge(sampleList,svClusters %>% filter(ClusterCount>50,ResolvedType!='LINE') %>% group_by(SampleId) %>% count(),by='SampleId',all.x=T) %>% replace_na(list(n=0)),aes(n)) + 
  stat_ecdf(geom = "step", pad = FALSE) + ylim(0,1) + labs(title = 'CDF # of NON LINE Clusters with ClusterCount>50 Per Sample')

##############################################################
########## SIMPLE SV SIGNATURES ##############################
##############################################################

#1. top50 simple
filter = cohortSummary(svData,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV')",'SampleId') %>% top_n(50,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','AllDel&Dup',useLogY =F)

#2. SHORT DELS: <100 bases length is related 1k to 10k bucket
scatterPlot(svData %>% filter(Type=='DEL',ResolvedType=='SimpleSV') %>% group_by(SampleId) %>% summarise(DElLT100bases=sum(PosEnd-PosStart<100),DELLT1kto10kbases=sum(PosEnd-PosStart>1e3&PosEnd-PosStart<1e4)),
            'DElLT100bases','DELLT1kto10kbases',F,F) 

#3. ULTRA SHORT DUP: Top 20:
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length<100",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Ultra Short Dup',useLogY =F)

#4.ULTRA SHORT DEL: Top 20:
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length<100",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Ultra Short Del',useLogY =F)

#5. SHORT DEL: Top 30:
filter = cohortSummary(svData,"Type=='DEL',ResolvedType=='SimpleSV',Length>1000,Length<5000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','Ultra Short Del',useLogY =F)

#5. BRCA Length DUP (5k-100k): Top 30:
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>5000,Length<100000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','BRCA LENGTH DUP',useLogY =F)

#5. CDK12 Length DUP (100k+): Top 30:
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>100000",'SampleId') %>% top_n(30,count) %>% .$SampleId
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId %in% filter),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              %>%  mutate(ID = paste(CancerType,SampleId)) ,'LengthBucket','ID','CDK12 LENGTH DUP',useLogY =F)


#CDK12 AND CCNE1 associated with 2 different Long DUP types
filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length>100000",'SampleId,CancerType') %>% top_n(1600,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDrivers %>% filter(gene %in% c('CDK12','CCNE1')),by.x='SampleId',by.y='sampleId',all.x=T))# %>% group_by(gene) %>% count())

filter = cohortSummary(svData,"Type=='DUP',ResolvedType=='SimpleSV',Length<100",'SampleId,CancerType') %>% top_n(20,count) %>% select(SampleId,CancerType,count)
View(merge((filter),svDrivers %>% filter(),by.x='SampleId',by.y='sampleId',all.x=T))# %>% group_by(gene) %>% count())

View(merge((filter),svDrivers %>% filter(),by.x='SampleId',by.y='sampleId',all.x=T))# %>% group_by(gene) %>% count())

head(svData %>% filter(ClusterId==70), table(LEStart, LEEnd))

###################################################
########### Driver Variants #######################
###################################################

# 1. TO DO: REVISIT when Charles repopulates this field
#driverVariants=svData %>% filter(DriverStart!=''|DriverEnd!='') %>% mutate(DriverInfo = pmax(DriverStart,DriverEnd)) %>% separate(DriverInfo,c('DriverType','DriverGene','DriverSubType'),sep = ';') %>% unite(DriverCombined,DriverType,DriverSubType)
#svData = svData %>% group_by(SampleId,ClusterId) %>% mutate(isDriver=sum(ifelse(DriverStart!=''|DriverEnd!='',1,0))>1) %>% ungroup()
#View(svData %>% group_by(ResolvedType,ClusterCountBucket,isDriver,ClusterId,SampleId) %>% count() %>% group_by(ResolvedType,ClusterCountBucket,isDriver) %>% count() %>% spread(isDriver,nn))
#View(svData %>% group_by(ResolvedType,ClusterCountBucket,isDriver)  %>% count() %>% spread(isDriver,n))


########################################
########## GENIC #######################
########################################

# TO DO: REVIVE THIS ANALYSIS AFTER CHARLES UPDATES WITH GENE ANNOTATIONS
#View(svData %>% group_by(ResolvedType,IsGenicStart) %>% count())# %>% spread(IsGenicStart,n) %>% mutate(proportion=`TRUE`/(`TRUE`+`FALSE`)))
#plot_count_by_bucket_and_type(cohortSummary(svData ,"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsFS","LengthBucket,IsGenic=IsGenicStart"),'LengthBucket','IsGenic','Length Distribution by isGenicEnd',useLogY =F)

##############################################################
############# INSERT SEQUENCE ################################
##############################################################

#1. INS SEQUENCE DISTRIBTUTION - 
# TO DO: Need to check after Daniel fix
plot_count_by_bucket_and_type(cohortSummary(svData,"ResolvedType!='Line',nchar(InsertSeq)<10,nchar(Homology)==0","Insert,ResolvedType"),'InsLen','ResolvedType','',F,T)
plot_count_by_bucket_and_type(cohortSummary(svData,"nchar(Homology)==0","InsLenBucket,ResolvedType"),'InsLenBucket','ResolvedType','',T,F)

#2. Frequency of short distribution is relatively random - rate of As and Ts is proportional to ref genome rate.
# TO DO: Need to check after Daniel fix
View(svData %>% filter(!IsLowQual,ResolvedType!='Line',nchar(InsertSeq)<=3) %>% group_by((substring(InsertSeq,1,1))) %>% count())

#3. Source of long insert sequences
# TO DO: How much is local vs remote?   How much is human vs non-human ref genome

#4. VIRAL INSERSTIONS
# TO DO: Count frequency per viral genome cancer type once Charles has loaded BEALN


##############################################################
############# SUPSICIOUS Findings ############################
##############################################################

#1. SUSPECT SGL EVENTS
# TO DO: Check if insert sequences match viruseses.   If not will probably filter in PUPRPLE where CN change < 0.15 and SVType = SGL and no immediate deletion bidge
View(svData %>% filter(AdjCNStart>0,!IsLowQual) %>% group_by(SampleId,ChrStart) %>% summarise(count=n(),countSGL=sum(ifelse(Type=='SGL',1,0)),
             countLT0.2VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.2,1,0)),countLT0.15VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.15,1,0)),
             countLT0.1VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.1,1,0)),maxCN=max(AdjCNStart)) %>% arrange(-countLT0.15VafSGL))
View(svData %>% filter(AdjCNStart>0,!IsLowQual) %>% group_by(ResolvedType,PolyA=grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)) %>% summarise(count=n(),countSGL=sum(ifelse(Type=='SGL',1,0)),
           countLT0.2VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.2,1,0)),countLT0.15VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.15,1,0)),
           countLT0.1VafSGL=sum(ifelse(Type=='SGL'&AdjCNChgStart/AdjCNStart<0.1,1,0)),maxCN=max(AdjCNStart)) %>% arrange(-countLT0.15VafSGL))

#2. What to do about POLYG INS SEQ?
# Mostly SGL in high copy number regions
# TO DO: Decide if additional filtering rule needed apart from #1 above
View(svData %>% filter(grepl('GGGGGGGGGGGG',InsertSeq)|grepl('CCCCCCCCCCCC',InsertSeq),!IsLowQual) %>% group_by(round(AdjCNStart,-1),Type) %>% count() %>% spread(Type,n))


#3. POLYA in non LINE ELEMENTS
#TO DO: why so many POLY A and T in unexpected cluster types
View( dbData %>% group_by(IsLINE,ResolvedType,IsPolyA) %>% count() %>% spread(IsPolyA,n,fill=0))

#4. TO DO: Investigate events in general with Qual score < 450

#5. Recurrent Variants
# TO DO: check after reload of Daniel PON update
View(svData %>% filter(PosEnd != -1,SampleId %in% highestPurityCohortSummary$sampleId) %>% group_by(ChrStart,ChrEnd,PosStart,PosEnd,Type,Length=PosEnd-PosStart) %>% 
       summarise(n=n(),min(Ploidy),max(Ploidy),mean(Ploidy),min(AdjCNChgStart),max(AdjCNChgStart),mean(AdjCNChgStart),min(AdjCNChgEnd),max(AdjCNChgEnd),mean(AdjCNChgEnd)) %>% filter(n>1))

#6. Overlapping DB are associated with much shorter INV foldback events - these may be some other type of event????
print(ggplot(data = foldbacks %>% filter(FoldbackType=="INV") %>% group_by(FoldbackLenBucket,DBOverlap=((DBLenStart>-30&DBLenStart==0)|(DBLenEnd>-30&DBLenEnd<=0))) %>% summarise(Count=n()) %>% spread(DBOverlap,Count), 
             aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=`TRUE`, colour='Short Overlapping DB'))
      + geom_line(aes(y=`FALSE`, colour='Other' ))
      + theme(panel.grid.major = element_line(colour="grey", size=0.5),panel.grid.minor.x = element_line(colour="grey", size=0.5))
      + scale_x_log10() + scale_y_log10()
      + labs(title = "Foldback Length Distribution"))

#7. LONE INV - vast majority are Foldbacks
View(svData %>% filter(ClusterCount==1,Type =='INV',!IsLowQual,ArmStart==ArmEnd) %>% arrange(SampleId,ChrStart,PosStart) %>%
       group_by(IsFoldBack,LengthBucket) %>% count() %>% spread(IsFoldBack,n))
#select(ResolvedType,len,Ploidy,Type,ClusterCount,ClusterReason,everything()))

#8. LONE BND - unclear what these are
View(svData %>% filter(ClusterCount==1,Type =='BND',!IsLowQual)%>%  mutate(len=PosEnd-PosStart) %>% arrange(SampleId,ChrStart,PosStart) %>%
       group_by(ResolvedType,PloidyBucket) %>% count() %>% spread(ResolvedType,n))

#9. Should we still have imprecise?
View(svData %>% group_by(Imprecise) %>% count())

#############################################
########## INDIVIDUAL SAMPLE ANALYSIS #######
#############################################

chartSample = 'CPCT02010794T'; chartChr = '6';

#0. Length distribution 
plot_count_by_bucket_and_type(cohortSummary(svData%>% filter(SampleId ==chartSample),"(ResolvedType=='SimpleSV')|(ResolvedType=='DEL_Int_TI'&Type=='INV'),!IsLowQual",'SampleId,CancerType,LengthBucket') 
                              ,'LengthBucket','SampleId','',useLogY =F)

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
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% count() %>% spread(chr_arm_end,n,fill=""))

View(svData %>% filter(SampleId == chartSample,IsFoldBack) %>% unite(chr_arm_start,ChrStart,ArmStart) %>% unite(chr_arm_end,ChrEnd,ArmEnd) %>%
       group_by(ClusterId,ResolvedType,SampleId,chr_arm_start,chr_arm_end) %>% count() %>% spread(chr_arm_end,n,fill=""))

# 3. Single Cluster View
View(svData %>% filter(SampleId==chartSample,ClusterId==385)%>% arrange(ChrStart,PosStart) 
     %>% select(ResolvedType,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,Type,ClusterCount,ClusterReason,everything()))

# 4. Large cluster View
View(svData %>% filter(!IsLowQual,SampleId==chartSample,ClusterCount>=1,ClusterCount>0,ResolvedType!='ASimpleSV')%>% 
       select(ResolvedType,Length,Ploidy,Type,ClusterCount,ClusterReason,everything()) %>% arrange(-Length))#%>% group_by(ClusterCountBucket,CnChStartBucket) %>% count() %>% spread(CnChStartBucket,n))







############################################### EXPERIMENTAL BELOW HERE ##############################################



##################################################
############### SV CLUSTERS ######################
##################################################

View(svClusters %>% filter(ClusterCount==3) %>% group_by(AcComplexFoldback,AcSimpleFoldback,AcComplexOther,AcMultipleDsb,AcDsb,AcRemoteTI,AcSoloSv) %>% count())
scatterPlot(svClusters %>% filter(ClusterCount>5),'AcComplexFoldback','AcSimpleFoldback')
scatterPlot(svClusters %>% filter(ClusterCount>5),'AcMultipleDsb','AcDsb')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcDsb')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcRemoteTI')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','AcSoloSv')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','ShortTIRemotes')
scatterPlot(svClusters %>% filter(ClusterCount>5),'ClusterCount','UnlinkedRemotes')

View(svClusters %>% filter(ClusterCount>0,ResolvedType!='Line',AcSimpleFoldback>0) %>% group_by(SampleId) %>% count())

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



############################################################
########### Basic Signature Charting #######################
############################################################

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

