detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
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

get_chromosome_copy_number<-function(dbConnect)
{
  query = paste("select purity.sampleId,cancerType,purity.gender,chromosome, ploidy,purity, ",
    "round(sum(bafCount*actualBaf*copyNumber)/sum(bafCount),1) as lwMajorAlleleAvg, ",
    "round(sum((end-start)*copyNumber)/sum(end-start),1) as lwCopyNumberAvg, ",
    "round(sum(if(copyNumber>=2.5,end-start,0))/sum(end-start),2) as amplifiedProportion, ",
    "round(sum(if(copyNumber*actualBaf>=1.5,end-start,0))/sum(end-start),2) as majorAlleleAmplifiedProportion, ",
    "round(sum(if(copyNumber<=1.5,end-start,0))/sum(end-start),2) as deletedProportion, ",
    "round(sum(if(segmentStartSupport='TELOMERE',copyNumber,0)),1) as pT, ",
    "round(sum(if(segmentEndSupport='CENTROMERE',copyNumber,0)),1) as pC, ",
    "round(sum(if(segmentStartSupport='CENTROMERE',copyNumber,0)),1) as qC, ",
    "round(sum(if(segmentEndSupport='TELOMERE',copyNumber,0)),1) as qT ",
    "from copyNumber inner join purity on copyNumber.sampleId = purity.sampleId ",
    "inner join clinical on copyNumber.sampleId = clinical.sampleId ",
    "where status <> 'NO_TUMOR' and qcStatus = 'PASS' ",
    "group by 1,2,3,4,5,6;",sep='')
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
}

bucketSummary<-function(cohort,bucketString = "")
{
  summary = (cohort %>% s_group_by(paste(bucketString,'chromosome',sep=','))
             %>% summarise(count=n())
             %>% spread(bucketString,count,fill=0)
             %>% arrange(chromosome) %>% as.data.frame)
}

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

chrCopyNumber<-get_chromosome_copy_number(dbProd)
chrCopyNumber$centromereChange = chrCopyNumber$qC-chrCopyNumber$pC
chrCopyNumber$telomereChange = chrCopyNumber$qT-chrCopyNumber$pT
chrCopyNumber$centromereChangeBucket=ifelse(abs(chrCopyNumber$centromereChange)<0.5,0,sign(chrCopyNumber$centromereChange)*2**(pmin(3,pmax(0,round(log(abs(chrCopyNumber$centromereChange),2),0)))))
chrCopyNumber$telomereChangeBucket=ifelse(abs(chrCopyNumber$centromereChange)<0.5,0,sign(chrCopyNumber$telomereChange)*2**(pmin(3,pmax(0,round(log(abs(chrCopyNumber$telomereChange),2),0)))))
chrCopyNumber$pTBucket=ifelse(abs(chrCopyNumber$pT)<0.5,0,2**(pmin(3,pmax(0,round(log(chrCopyNumber$pT,2),0)))))
chrCopyNumber$qTBucket=ifelse(abs(chrCopyNumber$qT)<0.5,0,2**(pmin(3,pmax(0,round(log(chrCopyNumber$qT,2),0)))))
chrCopyNumber$qCBucket=ifelse(abs(chrCopyNumber$qC)<0.5,0,2**(pmin(3,pmax(0,round(log(chrCopyNumber$qC,2),0)))))
chrCopyNumber$pCBucket=ifelse(abs(chrCopyNumber$pC)<0.5,0,2**(pmin(3,pmax(0,round(log(chrCopyNumber$pC,2),0)))))
chrCopyNumber$pTRelBucket=ifelse(abs(chrCopyNumber$pT/chrCopyNumber$ploidy)<0.2,0,2**(pmin(3,pmax(-2,round(log(chrCopyNumber$pT/chrCopyNumber$ploidy,2),0)))))
chrCopyNumber$qTRelBucket=ifelse(abs(chrCopyNumber$qT/chrCopyNumber$ploidy)<0.2,0,2**(pmin(3,pmax(-2,round(log(chrCopyNumber$qT/chrCopyNumber$ploidy,2),0)))))
chrCopyNumber$pCRelBucket=ifelse(abs(chrCopyNumber$pC/chrCopyNumber$ploidy)<0.2,0,2**(pmin(3,pmax(-2,round(log(chrCopyNumber$pC/chrCopyNumber$ploidy,2),0)))))
chrCopyNumber$qCRelBucket=ifelse(abs(chrCopyNumber$qC/chrCopyNumber$ploidy)<0.2,0,2**(pmin(3,pmax(-2,round(log(chrCopyNumber$qC/chrCopyNumber$ploidy,2),0)))))

View(bucketSummary(chrCopyNumber,'pCRelBucket'))
View(bucketSummary(chrCopyNumber,'qCRelBucket'))
View(bucketSummary(chrCopyNumber,'pTRelBucket'))
View(bucketSummary(chrCopyNumber,'qTRelBucket'))

View(bucketSummary(chrCopyNumber,'telomereChangeBucket'))
View(bucketSummary(chrCopyNumber,'centromereChangeBucket'))
View(bucketSummary(chrCopyNumber,'qTBucket'))

summary = (chrCopyNumber %>% group_by(telomereChangeBucket,chromosome)
           %>% summarise(count=n())
           %>% spread(telomereChangeBucket,count,fill=0)
           %>% arrange(chromosome) %>% as.data.frame)
View(summary)

head(chrCopyNumber)


print(ggplot(data=countsData,aes_string(x=bucket))+geom_line(aes(y=countDEL,colour='DEL'))
      +geom_line(aes(y=countDUP,colour='DUP'))+geom_line(aes(y=countINV,colour='INV'))+
        scale_x_log10()+facet_wrap(as.formula(paste("~", facetWrap)))+labs(title = titleString))