library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

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



#LOAD and ADD Buckets
cluster = read.csv('~/hmf/analyses/cluster/CLUSTER_V6.csv')
cluster$ClusterCountBucket=2**(round(log(cluster$ClusterCount,2),0))
cluster$PloidyBucket=2**(pmin(7,pmax(-3,round(log(cluster$Ploidy,2),0))))

########### CLUSTER COUNT ANALYSIS #################


countSummary = (cluster %>% filter(PONCount==0) 
                %>% group_by(ClusterCount,Annotations,inPONRegion=PONRegionCount>0) 
                %>% summarise(count=n())
                %>% unite(PONAnn,Annotations,inPONRegion)
                %>% spread(PONAnn,count))
View(countSummary)


########## INDIVIDUAL SAMPLE ANALYSIS   ################

# By Cluster
sample = 'DRUP01050018T'
sampleSummary = (cluster %>% filter(PONCount==0,SampleId==sample) 
                 %>% group_by(ClusterId)#,len=ifelse(Type=='BND',0,round(PosEnd-PosStart,-6))) 
                 %>% summarise(count=n(),
                               countLE=sum(LEStart=='true'|LEEnd=='true'),
                               countFS=sum(FSStart=='true'|FSEnd=='true'),
                               countDupBE=sum(DupBEStart=='true'|DupBEEnd=='true'),
                               countArm=mean(ArmCount),
                               countBND=sum(Type=='BND'),
                               countINV=sum(Type=='INV'),
                               countDEL=sum(Type=='DEL'),
                               countDUP=sum(Type=='DUP'),
                               countRegionPon=sum(PONRegionCount>0)) 
                 %>% arrange(-count) %>% as.data.frame)
View(sampleSummary)

# By Cluster Count
countSummary = (cluster %>% filter(PONCount==0,SampleId==sample) 
                %>% group_by(ClusterCount,Desc,Annotations) 
                %>% summarise(count=n())
                %>% spread(Annotations,count))
View(countSummary)

# Detailed Cluster Analysis
clusterNum = '5'
clusterSummary = (cluster %>% filter(SampleId==sample,ClusterId==clusterNum) 
                 %>% arrange(ChrStart,PosStart) %>% as.data.frame)
View(clusterSummary)


############# Sample Overview #############
cohortSummary = (cluster %>% filter(PONCount==0) 
                 %>% group_by(SampleId)#,len=ifelse(Type=='BND',0,round(PosEnd-PosStart,-6))) 
     %>% summarise(count=n(),
                   countC=n_distinct(ClusterId),
                   countInSingleC=sum(ClusterCount<=1),
                   countInGT5C=sum(ClusterCount>=5),
                   countLE=sum(LEStart=='true'|LEEnd=='true'),
                   countFS=sum(FSStart=='true'|FSEnd=='true'),
                   countDupBE=sum(DupBEStart=='true'|DupBEEnd=='true'),
                   countBND=sum(Type=='BND'),
                   countINV=sum(Type=='INV'),
                   countDEL=sum(Type=='DEL'),
                   countDUP=sum(Type=='DUP'),
                   countRegionPon=sum(PONRegionCount>0)) 
     %>% arrange(SampleId) %>% as.data.frame)
View(cohortSummary)
#ggplot(aes(countCluster),data=cohortSummary) + stat_ecdf(geom = "step", pad = FALSE) +   facet_wrap( ~len )


###### HOTSPOTS:  FRAGILE SITES / LINE ELEMENTS / DUPLICATE ###########

# Overall Counts
View(cluster %>% group_by(hasLE=LEStart=='true'|LEEnd=='true',hasDupBE=DupBEStart=='true'|DupBEEnd=='true',hasFS=FSStart=='true'|FSEnd=='true',Type) 
     %>% summarise(count=n()) 
     %>% spread (Type,count,fill=0)
     %>% as.data.frame)

# HotSpotAnlalysis
breakEnd_Analysis_By_Filter_And_Group<-function(filterString = "",groupByString = "")
{
  
  filtered = cluster %>% s_filter(filterString)
  filtered1 = filtered[,c("ChrStart","PosStart","OrientStart","ArmStart","Ploidy","Type","SampleId")]
  filtered2 = filtered[,c("ChrEnd","PosEnd","OrientEnd","ArmEnd","Ploidy","Type","SampleId")]
  colnames(filtered2)<-colnames(filtered1)
  breakend=rbind(filtered1,filtered2)
  summary = spread(breakend %>% s_group_by(groupByString) %>% summarise(count=n()) ,Type,count,fill=0)
  summary$Total = rowSums(summary[,-(1:2)])
  summary %>% arrange(-Total)
  
}

# LINE AND DupBE 1k buckets  
View(breakEnd_Analysis_By_Filter_And_Group("LEStart=='true'|LEEnd=='true',PONCount<2","ChrStart,Type,position=round(PosStart,-3)"))
View(breakEnd_Analysis_By_Filter_And_Group("DupBEStart=='true'|DupBEEnd=='true',PONCount<2","ChrStart,Type,position=round(PosStart,-3)"))
# Fragile Site
View(breakEnd_Analysis_By_Filter_And_Group("FSStart=='true'|FSEnd=='true',PONCount<2","ChrStart,Type,position=round(PosStart,-6)"))

#NONE OF THE ABOVE
View(breakEnd_Analysis_By_Filter_And_Group("FSStart=='false',FSEnd=='false',LEStart=='false',LEEnd=='false',DupBEStart=='false',DupBEEnd=='false',Type=='BND'|PosEnd-PosStart<1e9,PONCount<2",
                                           "ChrStart,Type,position=round(PosStart,-6)"))



variants_by_bucket<-function(filterString = "",groupByString = "")
{
  
  cluster %>% s_filter(filterString)%>% 
    s_group_by(groupByString) %>% 
    summarise(count=n(),unique=n_distinct(SampleId,ClusterId)) %>% 
    as.data.frame
  
}

# Duplicate Sites
variants_by_bucket("PONCount<2,DupBEStart=='true'|DupBEEnd=='true'","ClusterCountBucket") 
variants_by_bucket("PONCount<2,DupBEStart=='true'|DupBEEnd=='true'","PloidyBucket") 

# Fragile Sites
variants_by_bucket("PONCount<2,FSStart=='true'|FSEnd=='true'","ClusterCountBucket") 
variants_by_bucket("PONCount<2,FSStart=='true'|FSEnd=='true'","PloidyBucket") 
variants_by_bucket("PONCount<2,FSStart=='true'|FSEnd=='true'","ChrStart") %>% arrange(-count)






######################### OLD - TO BE FIXED ##################

########### Proximity analysis ########################


# BY CLUSTER COUNT
count10k = cluster10k %>% group_by(SampleId) %>% summarise(count=n_distinct(ClusterId),countClust=sum(ClusterCount>2),countLE=sum(LEStart=='True'|LEEnd=='true')) %>% arrange(SampleId) %>% as.data.frame
count100k = cluster100k %>% group_by(SampleId) %>% summarise(count=n_distinct(ClusterId),countClust=sum(ClusterCount>2),countLE=sum(LEStart=='True'|LEEnd=='true')) %>% arrange(SampleId) %>% as.data.frame
comparison = (merge(count10k,count100k,by="SampleId",suffixes = c("10k","100k"),all=TRUE,fill=0))
comparison$diff= comparison$count10k-comparison$count100k
View(comparison)

View(cluster %>% filter(SampleId=='CPCT02010304T') %>% group_by(ClusterId) %>% summarise(count=n()))
View(cluster %>% filter(SampleId=='CPCT02010304T'))#,ClusterId==43))

# BY ARM COUNT
field = "ArmCount"
count10k = cluster10k %>% group_by(ArmCount) %>% summarise(count=n()) %>% arrange(ArmCount) %>% as.data.frame
count100k = cluster100k %>% group_by(ArmCount) %>% summarise(count=n()) %>% arrange(ArmCount) %>% as.data.frame
comparison = (merge(count10k,count100k,by="ArmCount",suffixes = c("10k","100k"),all=TRUE,fill=0))
comparison$diff= comparison$count10k-comparison$count100k
View(comparison)

### Single clusters BY TYPE AND Annotation
temp = cluster10k %>% 
  filter(ClusterCount==2) %>%
  group_by(Type,Annotations) %>%
  summarise(count=n())
temp = dcast(temp,Annotations ~ Type,value="count")
temp[is.na(temp)] <- 0
head(temp,20)


########### PON VALIDATION ANALYSIS ########################

# RECURRENT VARIANTS not IN PON
recurrentVariants = (cluster %>% filter(PONCount<2) 
                     %>% group_by(ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,Type,len=ifelse(Type=='BND',0,PosEnd-PosStart))
                     %>% summarise(count=n(),mean(ClusterCount))
                     %>% filter(count>2))
%>% arrange(-count)
View(recurrentVariants)

# PON FILTERING IN VARIANTS IN CLUSTERS >= 2
PONClusteredSummary = (cluster %>% filter(ClusterCount>1) 
                       %>% group_by(SampleId,ClusterId,ClusterCount,Desc)  #,len=ifelse(Type=='BND',0,round(PosEnd-PosStart,-6))) 
                       %>% summarise(count=n(),countInPON=sum(PONCount>1))
                       %>% filter(countInPON>0)
                       %>% arrange(ClusterCount,countInPON) %>% as.data.frame)
View(PONClusteredSummary)

### By TYPE & PON_COUNT
View(cluster %>%  filter(PONCount>=1) %>% group_by(PONCount,Type) %>% summarise(count=n()) %>% spread(Type,count))

#BY Ploidy, CopyNumberChange Bucket 
PONClusteredSummary = (cluster %>% filter(PONCount>1) 
                       %>% group_by(Type,cnEnd=round(pmin(2,pmax(AdjCNChgEnd,0)),1),cnStart=round(pmin(2,pmax(AdjCNChgStart,0)),1))  #,len=ifelse(Type=='BND',0,round(PosEnd-PosStart,-6))) 
                       %>% summarise(count=n()) %>% spread(cnEnd,count) %>% as.data.frame)
View(PONClusteredSummary)


########## CLUSTER STATISTICS    ##################

temp = filtered %>% filter(ClusterCount==2,TICount == 1)
head(temp[order_by(temp$TILens),])
temp$TILens<-as.numeric(temp$TILens)
#temp = temp[,c("Desc","TILens")] %>% arrange(TILens)
##############
temp =cluster %>% 
  filter(PONCount==0) %>%
  group_by(isDup=DupBEStart=='true'|DupBEEnd=='true',ClusterCount) %>%
  summarise(count=n())
temp = dcast(temp,ClusterCount ~ isDup,value="count")
temp[is.na(temp)] <- 0
View(temp)

View(cluster %>% filter(ClusterCount==108))

### Single clusters BY TYPE AND Annotation
temp =cluster %>% 
  filter(PONCount==0,ClusterCount==3) %>%
  group_by(isDup=DupBEStart=='true'|DupBEEnd=='true',Annotations) %>%
  summarise(count=n())
temp = dcast(temp,Annotations ~ isDup,value="count")
temp[is.na(temp)] <- 0
head(temp,20)

cluster %>% filter(DupBEStart=='true'|DupBEEnd=='true',ClusterCount==1)

