detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

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
  filtered1 = filtered[,c("ChrStart","PosStart","OrientStart","ArmStart","Ploidy","Type","SampleId","FSStart","LEStart")]
  filtered2 = filtered[,c("ChrEnd","PosEnd","OrientEnd","ArmEnd","Ploidy","Type","SampleId","FSEnd","LEEnd")]
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

cohortSummary<-function(filterString = "",groupByString = "")
{
  summary = (cluster %>% s_filter(filterString) %>% s_group_by(groupByString)
             %>% summarise(count=n(),
                           countInSingleC=sum(ClusterCount<=1),
                           countIn2To4C=sum(ClusterCount>=5),
                           countInGTE5C=sum(ClusterCount>=5),
                           countLE=sum(LEStart=='true'|LEEnd=='true'),
                           countFS=sum(FSStart=='true'|FSEnd=='true'),
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
      plots[[i]]<-ggplot(aes_string(col),data=countsData) + stat_ecdf(geom = "step", pad = FALSE)+scale_x_log10()
      i=i+1
    }
  }
  multiplot(plotlist=plots,cols=3)
}

perCancerTypePDF<-function(variantType) {
  plots <- list()
  i = 1
  for (ct in unique(cluster$cancerType)) {
    summary = (cluster %>% filter(PONCount<2,cancerType==ct,Type==variantType) %>% group_by(LengthBucket) %>% summarise(count=n()))
    plots[[i]]<-ggplot(data=summary,aes(LengthBucket,count))+geom_point()+scale_x_log10()
    i=i+1
  }
  multiplot(plotlist=plots,cols=5)
}


scatterCounts<-function(countsData,bucket) {
  plots <- list()
  i = 1
  for (col in colnames(countsData)) {
    if (substr(col,1,5)=='count') {
      plots[[i]]<-ggplot(data=countsData,aes_string(bucket,col))+geom_point()+scale_x_log10()+
      i=i+1
    }
  }
  multiplot(plotlist=plots,cols=3)
}


################################################################
##################   CODE STARTS HERE   ########################
################################################################
### 0. LOAD Data



#LOAD and ADD Buckets
cluster = read.csv('~/hmf/analyses/cluster/CLUSTER_V6.csv')
cluster$ClusterCountBucket=2**(round(log(cluster$ClusterCount,2),0))
cluster$PloidyBucket=2**(pmin(7,pmax(-3,round(log(cluster$Ploidy,2),0))))
cluster$LengthBucket=ifelse(cluster$Type=='BND'|cluster$Type=='INS'|cluster$PosEnd-cluster$PosStart==0,
                            0,2**(round(log(cluster$PosEnd-cluster$PosStart,2),0)))

# Enrich with Tumor type
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
clinical = purple::query_clinical_data(dbProd)[,c('sampleId','cancerType')]
cluster = (merge(cluster,clinical,by.x="SampleId",by.y="sampleId",all.x=TRUE))


################################################################
### 1. PON ANALYSIS

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
summary = cohortSummary('','SampleId')
cdfCounts(summary)

##### All Samples by Length
summary = cohortSummary('','LengthBucket')
scatterCounts(summary,'LengthBucket')


################################################################
### 3. By Cancer Type

##### Lengths by cancer Type
for (variantType in c('INV','DUP','DEL') ) {
  summary = (cluster %>% filter(Type==variantType,FSStart=='false'|FSEnd=='false') %>% group_by(LengthBucket,cancerType) %>% summarise(count=n()))
  print(ggplot(data=summary,aes(LengthBucket,count))+geom_point()+scale_x_log10()+facet_wrap(~cancerType)+labs(title=variantType))

}

##### DEL LENGTHS for FS by cancer Type
summary = (cluster %>% filter(Type=='DEL',FSStart=='true'|FSEnd=='true') %>% group_by(LengthBucket,cancerType) %>% summarise(count=n()))
print(ggplot(data=summary,aes(LengthBucket,count))+geom_point()+scale_x_log10()+facet_wrap(~cancerType)+labs(title=variantType))


################################################################
### 4. SIGNATURES

signature_by_sample_and_type<-function(filter,signatureName){
  summary = cohortSummary('','SampleId,cancerType,LengthBucket')
  summary = summary[summary$SampleId %in% filter,]
  summary$ID = paste(summary$cancerType,summary$SampleId)
  for (variantType in c('INV','DUP','DEL') ) {
    print(ggplot(data=summary,aes_string('LengthBucket',paste('count',variantType,sep='')))+
            geom_point()+scale_x_log10()+facet_wrap(~ID)+labs(title=paste(signatureName,"signature for",variantType)))
  }
}



##### LONG DUPS ############
filter = cohortSummary('Type=="BND"|LengthBucket>1e5&LengthBucket<5e6','SampleId,cancerType')
filter$excessDUP = filter$countDUP - 0.5 * filter$countINV -0.5 * filter$countBND
filter = filter %>% arrange (-excessDUP) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"Long Dup")

##### Short DELS ############
filter = cohortSummary('Type=="BND"|LengthBucket<1e4','SampleId,cancerType')
filter = filter %>% arrange (-countDEL) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"Short DEL")

##### BRCA RS3 ############
filter = cohortSummary('Type=="BND"|LengthBucket>4e3,LengthBucket<5e4','SampleId,cancerType')
filter = filter %>% arrange (-countDUP) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"BRCA")

##### High mid length DUP + DEL ############
filter = cohortSummary('Type=="BND"|LengthBucket>4e3,LengthBucket<5e5','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>0.5) %>% arrange (-countDUP) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"Mid Length Dup + DEL")

##### Ovarian short-mid DEL ############
filter = cohortSummary('Type=="BND"|LengthBucket>1e4,LengthBucket<5e4','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>3) %>% arrange (-countDEL) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"shortish DEL")

##### ChromosomalStress ############
filter = cohortSummary('Type=="BND"|LengthBucket>1e6','SampleId,cancerType')
filter = filter %>% filter(countINV/(countDUP+countDEL)>0.9) %>% arrange (-countINV) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"ChromosomalStress")

##### FS like DEL ############
filter = cohortSummary('Type=="BND"|LengthBucket>5e4,LengthBucket<5e5','SampleId,cancerType')
filter = cohortSummary('Type=="BND"|LengthBucket>5e4,LengthBucket<5e5','SampleId,cancerType')
filter = filter %>% filter(countDEL/countDUP>3)  %>% arrange (-countDEL) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"FS like DEL")

##### Esophagus ############
filter = cohortSummary('cancerType=="Esophagus"','SampleId,cancerType')
filter = filter %>% arrange (-count) %>% filter(row_number() <= 20) %>% .$SampleId
signature_by_sample_and_type(filter,"Esophagus")


################################################################
### 4. BE Analysis

# Overall Counts of LE, FS, Duplicate BE
View(cluster %>% group_by(hasLE=LEStart=='true'|LEEnd=='true',hasDupBE=DupBEStart=='true'|DupBEEnd=='true',hasFS=FSStart=='true'|FSEnd=='true',Type)
     %>% summarise(count=n())
     %>% spread (Type,count,fill=0)
     %>% as.data.frame)

# Sites Enriched with Long DUP
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',LengthBucket>1e5&LengthBucket<5e6",
                                                       "ChrStart,position=round(PosStart,-6)"))
breakendSummary$excessDup = breakendSummary$countDUP-breakendSummary$countINV/2
View(breakendSummary %>% arrange(-excessDup))

# HotSpots for BFB (shows Key Oncogene Amplification sites)
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',LengthBucket>1e6",
                                           "ChrStart,position=round(PosStart,-6)"))
breakendSummary$excessDel = breakendSummary$countDEL-breakendSummary$countINV/2
breakendSummary$excessDup = breakendSummary$countDUP-breakendSummary$countINV/2
View(breakendSummary %>% arrange(countINV))


# Long INV,DUP,DEL by Chr Arm
breakendSummary=(breakEnd_Analysis_By_Filter_And_Group("LEStart=='false',LEEnd=='false',FSStart=='false',FSEnd=='false',LengthBucket>1e6",
                                                       "ChrStart,ArmStart,SampleId"))
print(ggplot(data=breakendSummary,aes(countINV,countDEL))+geom_point())
print(ggplot(data=breakendSummary,aes(countINV,countDUP))+geom_point())


################################################################
### 4. Cluster Analysis

### Single clusters BY TYPE AND Annotation
temp =cluster %>%
  group_by(ClusterCountBucket,Annotations,isDup=DupBEStart=='true'|DupBEEnd=='true') %>%
  summarise(count=n())  %>%
  spread(isDup,count)
View(temp)