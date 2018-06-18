detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)

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

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_arrange = function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_mutate = function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}

#' Modified version of dplyr's summarise that uses string arguments
#' @export
s_summarise = function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
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

breakEnd_Analysis_By_Filter_And_Group<-function(filterString = "",groupByString = "")
{

  filtered = cluster %>% s_filter(filterString)
  filtered1 = filtered[,c("ChrStart","PosStart","OrientStart","ArmStart","Ploidy","Type","SampleId","FSStart","LEStart","ArmCountStart","ArmExpStart")]
  filtered2 = filtered[,c("ChrEnd","PosEnd","OrientEnd","ArmEnd","Ploidy","Type","SampleId","FSEnd","LEEnd","ArmCountEnd","ArmExpEnd")]
  colnames(filtered2)<-colnames(filtered1)
  breakend=rbind(filtered1,filtered2)
  summary = breakend %>% s_group_by(groupByString) %>%
    summarise(count=n(),countSample=n_distinct(SampleId),
              countLE=sum(LEStart=='true'),
              countFS=sum(FSStart=='true'),
              countBND=sum(Type=='BND'),
              countINV=sum(Type=='INV'),
              countDEL=sum(Type=='DEL'),
              countDUP=sum(Type=='DUP')) %>%
    arrange(-count)

}

cohortSummary<-function(cluster,filterString = "",groupByString = "")
{
  summary = (cluster %>% s_filter(filterString) %>% s_group_by(groupByString)
             %>% summarise(count=n(),
                           countInSingleC=sum(ClusterCount<=1),
                           countIn2To5C=sum(ClusterCount>1 & ClusterCount<=5),
                           countIn5To20C=sum(ClusterCount>5 & ClusterCount<=20),
                           countInGT20C=sum(ClusterCount>20),
                           countLE=sum(LEStart!='false'|LEEnd!='false'),
                           countFS=sum(FSStart!='false'|FSEnd!='false'),
                           countDupBE=sum(DupBEStart=='true'|DupBEEnd=='true'),
                           countBND=sum(Type=='BND'),
                           countINV=sum(Type=='INV'),
                           countDEL=sum(Type=='DEL'),
                           countDUP=sum(Type=='DUP'))
             %>% arrange(-count) %>% as.data.frame)
}

cdfCounts<-function(countsData) {
  plots <- list()
  i = 1
  for (col in colnames(countsData)) {
    if (substr(col,1,5)=='count') {
      plots[[i]]<-ggplot(aes_string(col),data=countsData) + stat_ecdf(geom = "step", pad = FALSE)
      i=i+1
    }
  }
  multiplot(plotlist=plots,cols=3)
}


scatterCounts<-function(countsData,bucket) {
  plots <- list()
  i = 1
  for (col in colnames(countsData)) {
    if (substr(col,1,5)=='count') {
      plots[[i]]<-ggplot(data=countsData,aes_string(bucket,col))+geom_line()+scale_x_log10()
      i=i+1
    }
  }
  multiplot(plotlist=plots,cols=3)
}

plot_count_by_bucket_and_type<-function(countsData,bucket,facetWrap,titleString ="",useLogX = TRUE) {
  plot <- ggplot(data=countsData,aes_string(x=bucket))+geom_line(aes(y=countDEL,colour='DEL'))+
    geom_line(aes(y=countDUP,colour='DUP'))+geom_line(aes(y=countINV,colour='INV'))+geom_line(aes(y=countBND,colour='BND'))+facet_wrap(as.formula(paste("~", facetWrap)))+labs(title = titleString)
  if (useLogX == TRUE) {
      plot<-plot+scale_x_log10()
  }
  print(plot)
}

signature_by_sample_and_type<-function(cluster,filter,signatureName,bucket='LengthBucket'){
  summary = cohortSummary(cluster,'',paste('SampleId,primaryTumorLocation,',bucket,sep=""))
  summary = summary[summary$SampleId %in% filter,]
  summary$ID = paste(summary$cancerType,summary$SampleId)
  plot_count_by_bucket_and_type(summary,bucket,'ID',paste('SIGNATURE:',signatureName))
}


################################################################
##################   CODE STARTS HERE   ########################
################################################################
### 0. LOAD Data

#LOAD and ADD Buckets
cluster = read.csv('~/hmf/analyses/cluster/CLUSTER_V24.csv')
#cluster2 = cluster %>% separate(ChrArmStats,c('ArmStartBECount','ArmEndBECount','ArmMedianBECount'),sep=":")
cluster$ClusterCountBucket=2**(round(log(cluster$ClusterCount,2),0))
cluster$ArmCountBucket=2**(round(log(cluster$ArmCountStart,2),0))
cluster$HomLenBucket=2**(round(log(stri_length(cluster$Homology),2),0))
cluster$stressedArm=(cluster$ArmCountStart>1.3*cluster$ArmExpStart+6|cluster$ArmCountEnd>1.3*cluster$ArmExpEnd+6)
cluster$CNStartChBucket=2**(pmin(7,pmax(-3,round(log(cluster$AdjCNChgStart,2),0))))
cluster$CNEndChBucket=2**(pmin(7,pmax(-3,round(log(cluster$AdjCNChgEnd,2),0))))
cluster$PloidyBucket=2**(pmin(7,pmax(-3,round(log(cluster$Ploidy,2),0))))
#cluster$NearestLengthBucket=2**(pmin(20,pmax(0,round(log(cluster$NearestLen,2),0))))
#cluster$NearestTILengthBucket=2**(pmin(25,pmax(0,round(log(cluster$NearestTILen,2),0))))
#cluster$NearestDBLengthBucket=2**(pmin(25,pmax(0,round(log(cluster$NearestDBLen,2),0))))

cluster$LengthBucket=ifelse(cluster$Type=='BND'|cluster$Type=='INS'|cluster$PosEnd-cluster$PosStart==0|cluster$ArmEnd!=cluster$ArmStart,
                            0,2**(round(log(cluster$PosEnd-cluster$PosStart,2),0)))

# Enrich with Tumor type
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
clinical = purple::query_clinical_data(dbProd)[,c('sampleId','primaryTumorLocation')]
cluster = (merge(cluster,clinical,by.x="SampleId",by.y="sampleId",all.x=TRUE))


################################################################
### 1. PON ANALYSIS
######## MAKE THIS CONSISTENT #############

# The PON mainly removes short DELS
View(cluster %>% group_by(LengthBucket,Type,inPONRegion=PONRegionCount>1) %>% summarise(count=n()) %>% unite(PONType,Type,inPONRegion) %>% spread(PONType,count))

# 95% of the PON variants are in single clusters
View(cluster %>% group_by(ClusterCount=pmin(10,ClusterCount),Type,inPONRegion=PONRegionCount>1) %>% summarise(count=n()) %>% unite(PONType,Type,inPONRegion) %>% spread(PONType,count))

# The ploidies of the PON filtered variants are not necessarily low
View(cluster %>% group_by(PloidyBucket,Type,inPONRegion=PONRegionCount>1)  %>% summarise(count=n()) %>% unite(PONType,Type,inPONRegion) %>% spread(PONType,count))

# We still have a number of recurrentVariants (need to remove multiple biopsies from analysis)
recurrentVariants = (cluster %>% filter(PONCount<2) %>% group_by(ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,Type,len=ifelse(Type=='BND',0,PosEnd-PosStart))
                     %>% summarise(count=n(),mean(ClusterCount))
                     %>% filter(count>2))
View(recurrentVariants)

# FILTER FOR PONCount <2 for all subsequent analyses
cluster = cluster %>% filter(PONCount<2)

################################################################
### 2. Per Sample / Feature Counts

##### All Samples by sampleID
summary = cohortSummary(cluster,'','SampleId')
cdfCounts(summary)

##### All Samples by Length
summary = cohortSummary(cluster,'','LengthBucket')
scatterCounts(summary,'LengthBucket')


################################################################
### 3. By Cancer Type

##### Lengths by cancer Type
summary = cohortSummary(cluster,'','LengthBucket,primaryTumorLocation')
plot_count_by_bucket_and_type(summary,'LengthBucket','primaryTumorLocation',useLogX = T)

##### DEL LENGTHS for FS by cancer Type
summary = cohortSummary(cluster,"FSStart!='false'|FSEnd!='false'",'LengthBucket,cancerType')
plot_count_by_bucket_and_type(summary,'LengthBucket','cancerType')

summary = cohortSummary(cluster,'cancerType=="Colorectal"',"LengthBucket,IsFS=FSStart=='true'|FSEnd=='true'")
plot_count_by_bucket_and_type(summary,'LengthBucket','IsFS')
head(cluster)

################################################################
### 4. SIGNATURES

##### LONG DUPS ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>1e5&LengthBucket<5e6','SampleId,primaryTumorLocation')
filter$excessDUP = filter$countDUP - 0.5 * filter$countINV -0.5 * filter$countBND
filter = filter %>% arrange (-excessDUP) %>% filter(row_number() <= 30) %>% .$SampleId
#### TOPN(30)
signature_by_sample_and_type(cluster,filter,"Long Dup")



##### Short DELS ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket<1e4','SampleId,cancerType')
filter = filter %>% arrange (-countDEL) %>% filter(row_number()<=30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"Short DEL")

##### BRCA RS3 ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>4e3,LengthBucket<5e4','SampleId,cancerType')
filter = filter %>% arrange (-countDUP) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"BRCA")

##### High mid length DUP + DEL ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>4e3,LengthBucket<5e5','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>0.5) %>% arrange (-countDUP) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(filter,"Mid Length Dup + DEL")

##### Ovarian short-mid DEL ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>1e4,LengthBucket<5e4','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>3) %>% arrange (-countDEL) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(filter,"shortish DEL")

##### ChromosomalStress ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>1e6','SampleId,cancerType')
filter = filter %>% filter(countINV/(countDUP+countDEL)>0.9) %>% arrange (-countINV) %>% filter(row_number() <=30) %>% .$SampleId
signature_by_sample_and_type(filter,"ChromosomalStress")

##### FS like DEL ############
filter = cohortSummary(cluster,'Type=="BND"|LengthBucket>5e4,LengthBucket<5e5','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>3)  %>% arrange (-countDEL) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"FS like DEL")

##### Esophagus ############
filter = cohortSummary(cluster,'cancerType=="Ovary"','SampleId,cancerType')
filter = filter %>% arrange (-count) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"Ovary")

##### Overall ############
filter = cohortSummary(cluster,'','SampleId,cancerType')
filter = filter %>% arrange (-count) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"ALL")

##### Colorectal ############
filter = cohortSummary(cluster,'cancerType=="Colorectal"','SampleId,cancerType')
View(filter)
filter = filter %>% arrange (-count) %>% filter(row_number() <=30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"Colorectal")

################################################################
### 4. BE Analysis

# Overall Counts of LE, FS, Duplicate BE
View(cluster %>% group_by(hasLE=LEStart=='true'|LEEnd=='true',potentialhasLE=LEStart=='ident'|LEEnd=='ident',hasFS=FSStart=='true'|FSEnd=='true',Type) %>%
       summarise(count=n()) %>%
       spread (Type,count,fill=0) %>%
       as.data.frame)

# Sites Enriched with Long DUP
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',LengthBucket>1e5&LengthBucket<5e6",
                                                       "ChrStart,position=round(PosStart,-6)"))
breakendSummary$excessDup = breakendSummary$countDUP-breakendSummary$countINV/2
View(breakendSummary %>% arrange(-excessDup))

# HotSpots for BFB (shows Key Oncogene Amplification sites)
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',ChrStart==5",
                                           "ChrStart,position=round(PosStart,-6)"))
breakendSummary$excessDel = breakendSummary$countDEL-breakendSummary$countINV/2
breakendSummary$excessDup = breakendSummary$countDUP-breakendSummary$countINV/2
View(breakendSummary %>% arrange(-countINV))



# Long INV,DUP,DEL by Chr Arm
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',LengthBucket>1e6",
                                                       "ChrStart,ArmStart,SampleId"))
print(ggplot(data=breakendSummary,aes(countINV,countDEL))+geom_point())
print(ggplot(data=breakendSummary,aes(countINV,countDUP))+geom_point())


################################################################
### 5. Cluster Analysis

### Single clusters BY TYPE AND Annotation
temp =cluster %>%
  group_by(ClusterCountBucket,Annotations,isDup=DupBEStart=='true'|DupBEEnd=='true') %>%
  summarise(count=n())  %>%
  spread(isDup,count)
View(temp)

################################################################
### 6. Neighbour Proximity

temp = cluster %>% group_by(SampleId,ClusterId) %>% mutate(ClusterChrCount=n_distinct(ChrStart,ChrEnd)) %>% ungroup() %>% as.data.frame

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLen>0,ClusterCount==1",'NearestLengthBucket,ClusterCountBucket')
plot_count_by_bucket_and_type(summary,'NearestLengthBucket','ClusterCountBucket')

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount==3,stressedArm==FALSE",'NearestLengthBucket,Desc')
plot_count_by_bucket_and_type(summary,'NearestLengthBucket','Desc')

summary = cohortSummary(temp,"FSStart=='false',FSEnd=='false',LengthBucket>=0,NearestLen>=0,
                        ClusterCount==2,NearestTILengthBucket<NearestDBLengthBucket",'NearestTILengthBucket,ClusterChrCount')
plot_count_by_bucket_and_type(summary,'NearestTILengthBucket','ClusterChrCount',,TRUE)

summary = cohortSummary(temp,"FSStart=='false',FSEnd=='false',stressedArm==T,LengthBucket>=0,NearestLen>=0,
                        ClusterCount<=2,NearestDBLengthBucket<=5e6",'NearestDBLengthBucket,ClusterChrCount')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','ClusterChrCount',,TRUE)

summary = cohortSummary(temp,"FSStart=='false',FSEnd=='false',LengthBucket>=0,NearestLen>=0,
                        ClusterCount<=200,NearestDBLengthBucket<=3e6,ArmCountStart<200,ArmCountEnd<200",'NearestDBLengthBucket,stressedArm')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','stressedArm',,TRUE)

summary = cohortSummary(temp,"FSStart=='false',FSEnd=='false',LengthBucket>=0,NearestLen==0|NearestLen>30,
                        ClusterCount<=200,NearestTILengthBucket<=3e6,ArmCountStart<20,ArmCountEnd<20",'NearestTILengthBucket,stressedArm')
plot_count_by_bucket_and_type(summary,'NearestTILengthBucket','stressedArm',,TRUE)

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',stressedArm==F,LengthBucket>=0,NearestLen>0,ClusterCount==2,NearestDBLengthBucket>NearestTILengthBucket",'NearestTILengthBucket,ClusterDesc')
plot_count_by_bucket_and_type(summary,'NearestTILengthBucket','ClusterDesc',,TRUE)

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',stressedArm==F,LengthBucket>=0,NearestLen>0,ClusterCount==3,NearestDBLengthBucket<NearestTILengthBucket",'NearestDBLengthBucket,ClusterDesc')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','ClusterDesc',,TRUE)

 summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount==2,stressedArm==FALSE",'NearestDBLengthBucket,Desc')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','Desc')


summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount==3,stressedArm==FALSE",'NearestLengthBucket,NearestType')
plot_count_by_bucket_and_type(summary,'NearestLengthBucket','Desc',,TRUE)

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount>=2,ClusterCount<=2",'NearestTILengthBucket,ArmCountBucket')
plot_count_by_bucket_and_type(summary,'NearestTILengthBucket','ArmCountBucket',,TRUE)

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount>=2,ClusterCount<=10,stressedArm==FALSE",'NearestDBLengthBucket,ArmCountBucket')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','ArmCountBucket',,TRUE)

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount>=2,ClusterCount<=2,stressedArm==FALSE,NearestDBLengthBucket<NearestTILengthBucket",'NearestDBLengthBucket,ArmCountBucket')
plot_count_by_bucket_and_type(summary,'NearestDBLengthBucket','ArmCountBucket',,TRUE)

##############

# SINGLE ARM INVERSIONS
breakendSummary=breakEnd_Analysis_By_Filter_And_Group("ArmCountStart==2","SampleId,ChrStart,ArmStart")
View(breakendSummary %>%filter(countINV==count) %>% group_by(SampleId) %>% summarise(count=n()) %>% arrange(-count))
View(breakendSummary %>%filter(countINV==count) %>% group_by(SampleId) %>% summarise(count=n()) %>% arrange(-count))

summary = cohortSummary(cluster,"ArmCountStart==4,ClusterCount>=1,Desc!='CRS'",'PloidyBucket,Type')
plot_count_by_bucket_and_type(summary,'PloidyBucket','Type',,TRUE)

summary = cohortSummary(cluster,"ArmCountStart>=2,ArmCountStart<=10,ClusterCount>1,Desc!='CRS',Type=='DUP'",'CNStartChBucket,LengthBucket')
plot_count_by_bucket_and_type(summary,'CNStartChBucket','LengthBucket',,TRUE)

summary = cohortSummary(cluster,"ArmCountStart>=2,ArmCountStart<=5000,ClusterCount>1,Desc!='CRS',Type=='INV'",'CNStartChBucket,LengthBucket')
plot_count_by_bucket_and_type(summary,'CNStartChBucket','LengthBucket',,TRUE)



#############

filter = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount<=10,stressedArm==FALSE,NearestDBLengthBucket<NearestTILengthBucket,NearestDBLengthBucket<1e5",'SampleId,cancerType')
filter = filter %>% arrange (-count) %>% filter(row_number() <= 30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"ALL",'NearestDBLengthBucket')


summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',NearestLength>0,ClusterCount==3",'ArmCountBucket,Desc')
plot_count_by_bucket_and_type(summary,'ArmCountBucket','Desc',,TRUE)



################
###### 7. Short INV

summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',PloidyBucket>8","LengthBucket,stressedArm")
plot_count_by_bucket_and_type(summary,'LengthBucket','stressedArm',FALSE)


summary = cohortSummary(cluster,"FSStart=='false',FSEnd=='false',LengthBucket<3e3,PloidyBucket<2","CNStartChBucket,stressedArm=ArmCountStart>1.3*ArmExpStart+6|ArmCountEnd>1.3*ArmExpEnd+6")
plot_count_by_bucket_and_type(summary,'CNStartChBucket','stressedArm',FALSE)



##RANDOM
cluster %>% filter(SampleId =='CPCT02070300T')

summary = cohortSummary(cluster,"SampleId=='CPCT02070300T'","LengthBucket")
plot_count_by_bucket_and_type(summary,'LengthBucket','stressedArm',FALSE)

summary = cohortSummary(cluster,'','ArmStartBECount')
View(summary)
head(cluster[order(cluster$ArmStartBECount),])

temp = cohortSummary(cluster,'','SampleId,')
cdfCounts(temp)

temp = (cluster %>% filter(ChrEnd==5))
head(temp)
ggplot(aes(PosEnd),data=cluster) + stat_ecdf(geom = "step", pad = FALSE)+facet_wrap(~ChrEnd)

##### Analysis of Nearest Length
head(cluster)
ggplot(aes(NearestLength),data=cluster) + stat_ecdf(geom = "step", pad = FALSE) +scale_x_log10() + facet_wrap(~NearestType)
ggplot(aes(NearestLength),data=cluster) + stat_ecdf(geom = "step", pad = FALSE) +scale_x_log10() + facet_wrap(~ClusterCountBucket)


##### Cluster filtering by

cluster2 = cluster %>% filter(ArmCountStart<1.5*ArmExpStart+10&ArmCountEnd<1.5*ArmExpEnd+10)

##### NON Chromosomal Stress with long variants ############
filter = cohortSummary(cluster,'ArmCountStart<2*ArmExpStart+5,ArmCountEnd<2*ArmExpEnd+5,Type=="BND"|LengthBucket>1e6','SampleId,cancerType')
filter = filter %>% filter(countINV/(countDUP+countDEL)>0.9) %>% arrange (-countINV) %>% filter(row_number() <=30) %>% .$SampleId
signature_by_sample_and_type(cluster,filter,"Long non chromosomal stress")


#######################

breakEnd_Analysis_By_Filter_And_Group<-function(filterString = "",groupByString = "")
{

  filtered = cluster %>% s_filter(filterString)
  filtered1 = filtered[,c("ChrStart","PosStart","OrientStart","ArmStart","Ploidy","Type","SampleId","FSStart","LEStart","ArmCountStart","ArmExpStart")]
  filtered2 = filtered[,c("ChrEnd","PosEnd","OrientEnd","ArmEnd","Ploidy","Type","SampleId","FSEnd","LEEnd","ArmCountEnd","ArmExpEnd")]
  colnames(filtered2)<-colnames(filtered1)
  breakend=rbind(filtered1,filtered2)
  summary = breakend %>% s_group_by(groupByString) %>%
    summarise(count=n(),countSample=n_distinct(SampleId),
              countLE=sum(LEStart=='true'),
              countFS=sum(FSStart=='true'),
              countBND=sum(Type=='BND'),
              countINV=sum(Type=='INV'),
              countDEL=sum(Type=='DEL'),
              countDUP=sum(Type=='DUP')) %>%
    arrange(-count)

}

temp = cluster %>% filter(SampleId %in% filter,Type=='DUP')

temp = breakEnd_Analysis_By_Filter_And_Group("SampleId %in% filter,Type=='DUP'","round(PosStart,-6),ChrStart")
View(temp)
nrow(cluster2)

ggplot(aes(NearestLength),data=cluster[cluster$arm]) + stat_ecdf(geom = "step", pad = FALSE) +scale_x_log10() + facet_wrap(~ClusterCountBucket)


View(cluster %>% filter(SampleId=='CPCT02160019T') %>% group_by(ChrStart,ArmStart) %>% summarise(mean(ArmCountStart),mean(ArmExpStart)))

head(temp)
ggplot(data=temp,aes=(temp$DB,temp$TI))+geom_point()
ggplot(temp, aes(DB,TI)) + geom_point()

View(cluster %>% filter(SampleId=='CPCT02030256T',stressedArm==T) %>% group_by(ClusterId) %>% count())

summary = cohortSummary(cluster,"SampleId=='CPCT02330059T'",'LengthBucket,stressedArm')
plot_count_by_bucket_and_type(summary,'LengthBucket','stressedArm')

    cluster$DBorTI=ifelse(cluster$ClusterCount==1,'None',ifelse(cluster$ClusterCount>1&cluster$NearestDBLen>cluster$NearestTILen&cluster$NearestTILen>30,ifelse(cluster$NearestTILen<300,'shortTI','longTI'),'DB'))
View(cluster %>% filter(stressedArm==T,ArmStart==ArmEnd,ArmCountStart>20,ChrStart==ChrEnd)   %>% group_by(SampleId,ChrEnd,ArmEnd,DBorTI) %>% summarise(count=n()) %>% spread(DBorTI,count))
