library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
library(stringr)

# install.packages("devtools") 
# library(BiocManager)
# install.packages('argparser')
# install.packages("data.table")
#install.packages("hexbin")
library(hexbin)

#install.packages("ggplot2") 

sv_set_common_fields<-function(svData)
{
  svData = svData %>% mutate(IsLINE = ifelse(LEStart!='None'|LEEnd!='None',T,F),
                             IsFS = ifelse(FSStart!='false'|FSEnd!='false',T,F),
                             Length = ifelse(as.character(ChrStart)!=as.character(ChrEnd)|Type=='INS'|ArmEnd!=ArmStart, -1, PosEnd-PosStart),
                             TICount = ifelse(LnkLenStart>0,0.5,0)+ifelse(LnkLenEnd>0,0.5,0),
                             DBCount = ifelse(DBLenStart>=0,0.5,0)+ifelse(DBLenEnd>=0,0.5,0),
                             AsmbTICount = ifelse(grepl('asm',AsmbStart),0.5,0)+ifelse(grepl('asm',AsmbEnd),0.5,0),
                             InferTICount = TICount - AsmbTICount,
                             ClusterSize = ifelse(ClusterCount==1,'Single',ifelse(ClusterCount<=3,'Small','Large')),
                             IsChained = (ChainCount>=1),
                             FoldbackCount = ifelse(FoldbackLenStart>=0,0.5,0)+ifelse(FoldbackLenEnd>=0,0.5,0))
  return (svData)
}

sv_load_and_prepare<-function(filename, filterLowQual=T)
{
  svData = read.csv(filename)
  sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
  svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
  
  if(filterLowQual)
    svData = svData %>% filter(ResolvedType!='DUP_BE'&ResolvedType!='POLY_G_C')
  
  svData = sv_set_common_fields(svData)
  return (svData)  
}

clusters_load<-function(filename)
{
  clusters = read.csv(filename)
  sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
  clusters = merge(clusters, sampleCancerTypes, by='SampleId', all.x=T)
  return (clusters)  
}

plot_cross_chr_sample_clusters<-function(svData, chromosomes, specificSample, minClusterCount = 2, sizeByCount=T)
{
  relevantClusters = svData %>% filter(SampleId==specificSample) %>% filter(Type=='BND' & ChrStart %in% chromosomes & ChrEnd %in% chromosomes) %>% select(ClusterId)
  ccDataStart = svData %>% filter(SampleId==specificSample & ClusterId %in% relevantClusters$ClusterId & ChrStart %in% chromosomes)
  ccDataStart$Chr = ccDataStart$ChrStart
  ccDataStart$Pos = ccDataStart$PosStart
  ccDataEnd = svData %>% filter(SampleId==specificSample & ClusterId %in% relevantClusters$ClusterId & ChrEnd %in% chromosomes)
  ccDataEnd$Chr = ccDataEnd$ChrEnd
  ccDataEnd$Pos = ccDataEnd$PosEnd
  ccData = rbind(ccDataStart,ccDataEnd)
  ccData$PosBucket = round(ccData$Pos/1e6)
  ccData$Cluster = paste(ccData$ClusterId, ccData$ClusterDesc, sep='_')
  
  ccDataSum = ccData %>% filter(ClusterCount>=minClusterCount&grepl('BND',ClusterDesc))%>% group_by(Cluster,Chr,PosBucket) %>% summarise(Count=n())
  ccDataSum$Point = ifelse(sizeByCount,ccDataSum$Count,1)

  # clusterColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown", "darkgreen", "deepskyblue", "tan")
  
  print(ggplot(data = ccDataSum, aes(x=PosBucket, y=Chr, colour=Cluster)) 
        + geom_point(aes(size = Point, colour=factor(Cluster))))
}

plot_cross_chr_sample_clusters(svData, c(1,3,5,17), 'CPCT02050182T', 1, T)


view_cluster_sv<-function(sampleId, clusterId)
{
  View(svData %>% filter(SampleId==sampleId&ClusterId==clusterId))  
}

view_clusters_sv<-function(sampleId, clusterIds)
{
  View(svData %>% filter(SampleId==sampleId & ClusterId %in% clusterIds))  
}

view_sv_id<-function(svId)
{
  View(svData %>% filter(Id==svId))  
}

view_sv_ids<-function(svIdList)
{
  View(svData %>% filter(Id %in% svIdList))  
}

view_chromosome_sv<-function(sampleId, chrStart, chrEnd, requireBoth=F)
{
  data = svData %>% filter(SampleId==sampleId)
  
  if(requireBoth && chrStart != '' && chrEnd != '')
  {
    data = data %>% filter(ChrStart==chrStart&ChrEnd==chrEnd) %>% arrange(PosStart)
  }
  else if(chrStart != '' && chrEnd != '')
  {
    data = data %>% filter(ChrStart==chrStart | ChrEnd==chrEnd) %>% arrange(PosStart)
  }
  else if(chrStart != '')
  {
    data = data %>% filter(ChrStart==chrStart) %>% arrange(PosStart)
  }
  else if(chrEnd != '')
  {
    data = data %>% filter(ChrEnd==chrEnd) %>% arrange(PosEnd)
  }
  
  View(data)
}

view_chromosome_sv<-function(svData, sampleId, chrStart, chrEnd, requireBoth=F)
{
  data = svData %>% filter(SampleId==sampleId)
  
  if(requireBoth && chrStart != '' && chrEnd != '')
  {
    data = data %>% filter(ChrStart==chrStart&ChrEnd==chrEnd) %>% arrange(PosStart)
  }
  else if(chrStart != '' && chrEnd != '')
  {
    data = data %>% filter(ChrStart==chrStart | ChrEnd==chrEnd) %>% arrange(PosStart)
  }
  else if(chrStart != '')
  {
    data = data %>% filter(ChrStart==chrStart) %>% arrange(PosStart)
  }
  else if(chrEnd != '')
  {
    data = data %>% filter(ChrEnd==chrEnd) %>% arrange(PosEnd)
  }
  
  View(data)
}

view_cluster_arm_breakends<-function(svData, sampleId, clusterId, chromosome, arm)
{
  specCluster = svData %>% filter(SampleId==sampleId)
  
  if(clusterId!="")
    specCluster = svData %>% filter(ClusterId==clusterId)
  
  startData = specCluster %>% filter(ChrStart==chromosome)
              
  if(arm != '') 
    startData = startData %>% filter(ArmStart==arm)

  startData = startData %>% 
    mutate(Chr=ChrStart,Pos=PosStart,Orient=OrientStart,OtherPos=PosEnd,IsStart=T,CNChg=CNChgStart)
  
  endData = specCluster %>% filter(ChrEnd==chromosome)
  
  if(arm != '')
    endData = endData %>% filter(ArmEnd==arm)

  endData = endData %>% 
    mutate(Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,OtherPos=PosStart,IsStart=F,CNChg=CNChgEnd)
  
  specClusterChr = rbind(startData,endData) %>% 
    mutate(PloidyCalc=(PloidyMin+PloidyMax)*0.5) %>%
    select(ClusterId,Id,Chr,IsStart,Pos,OtherPos,Orient,Type,PloidyCalc,PloidyMin,PloidyMax,CNChg,CNStart,CNEnd,
           FoldbackLnkStart,FoldbackLnkEnd,CNChgStart,CNChgEnd,Ploidy,everything()) %>% arrange(Pos)
  
  return (specClusterChr)
}

view_cluster_sv_links<-function(sampleId, clusterId)
{
  View(svData %>% filter(SampleId==sampleId&ClusterId==clusterId)
       %>% select(SampleId,ClusterId,ClusterDesc,ResolvedType,Id,Type,ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd,LnkLenStart,LnkLenEnd,AsmbMatchStart,AsmbMatchEnd,
                  DBLenStart,DBLenEnd,ChainIndex,Consistency,FoldbackLnkStart,FoldbackLnkEnd,AdjCNStart,AdjCNChgStart,AdjCNEnd,AdjCNChgEnd,Ploidy))  
}

view_clusters_sv_links<-function(sampleId, clusterIds)
{
  View(svData %>% filter(SampleId==sampleId&ClusterId %in% clusterIds)
       %>% select(SampleId,ClusterId,ClusterDesc,ResolvedType,Id,Type,ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd,LnkLenStart,LnkLenEnd,AsmbMatchStart,AsmbMatchEnd,
                  DBLenStart,DBLenEnd,ChainIndex,Consistency,FoldbackLnkStart,FoldbackLnkEnd,AdjCNStart,AdjCNChgStart,AdjCNEnd,AdjCNChgEnd,Ploidy))  
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
# generic plotting functions
plot_length_facetted<-function(data, filterStr, groupByStr, bucketStr, facetStr, titleStr, logScaleX=T, logScaleY=F)
{
  if(filterStr != '')
    data = data %>% s_filter(filterStr)
  
  dataSummary = data %>% s_group_by(groupByStr) %>% summarise(Count=n()) %>% as.data.frame
  
  plot = (ggplot(data = dataSummary, aes_string(x=bucketStr, y='Count'))
          + geom_line()
          + labs(title = titleStr))
  
  if(facetStr != "")
    plot = plot + facet_wrap(as.formula(paste("~", facetStr)))
  
  if(logScaleX)
    plot = plot + scale_x_log10()
  
  if(logScaleY)
    plot = plot + scale_y_log10()

    print(plot)
}

classify_del_dup_length<-function(type,length)
{
  shortDupThreshold=500
  medDupThreshold=8e4
  longDupThreshold=1e6
  shortDelThreshold=500
  medDelThreshold=1e4
  longDelThreshold=5e5
  
  classifyStr = ifelse(type=='DEL',
                       ifelse(length<=shortDelThreshold,'DEL_SHORT',
                              ifelse(length>shortDelThreshold&length<medDelThreshold,'DEL_MEDIUM',
                                     ifelse(length>medDelThreshold&length<longDelThreshold,'DEL_LONG','DEL_V_LONG'))),
                       ifelse(type=='DUP',ifelse(length<=shortDupThreshold,'DUP_SHORT',
                                                 ifelse(length>shortDupThreshold&length<medDupThreshold,'DUP_MEDIUM',
                                                        ifelse(length>medDupThreshold&length<longDupThreshold,'DUP_LONG','DUP_V_LONG'))),'NONE'))

    return (classifyStr)
}

# print(classify_del_dup_length('DEL',500))
# print(classify_del_dup_length('DEL',1e6))
# print(classify_del_dup_length('DUP',1e5))


calc_fisher_et<-function(allSamples,g1Samples,g2Samples,g1Name,g2Name,log=T,returnDF=F,testLabel='')
{
  scAll = nrow(allSamples)
  scWithG1 = nrow(g1Samples)
  scWithG2 = nrow(g2Samples)
  scWithG1WithG2 = nrow(g1Samples %>% filter(SampleId %in% g2Samples$SampleId))
  scNoG1NoG2 = nrow(allSamples %>% filter(!(SampleId %in% g1Samples$SampleId)&!(SampleId %in% g2Samples$SampleId)))
  
  scWithG1NoG2 = nrow(g1Samples %>% filter(!(SampleId %in% g2Samples$SampleId)))
  scNoG1WithG2 = nrow(g2Samples %>% filter(!(SampleId %in% g1Samples$SampleId)))

  fishMatrix = rbind(c(scWithG1WithG2,scNoG1WithG2), c(scWithG1NoG2,scNoG1NoG2))
  
  # expected count of this Cancer within enriched samples
  expected = round(scWithG1 * scWithG2 / scAll, 2)
  fetProb = fisher.test(fishMatrix, alternative="greater")$p.value

  if(log)
  {
    print(sprintf("all=%d with%s=%d with%s=%d both=%d neither=%d with%sNo%s=%d no%sWith%s=%d", 
                  scAll,g1Name,scWithG1,g2Name,scWithG2,scWithG1WithG2,scNoG1NoG2,
                  g1Name,g2Name,scWithG1NoG2,g1Name,g2Name,scNoG1WithG2))

    print(sprintf("expected=%.2f prob=%f",expected, fetProb))
  }
  
  if(returnDF)
  {
    result = data.frame(matrix(ncol = 10, nrow = 0))
    colnames(result) = c('Test','All',g1Name,g2Name,'Both','Neither',
                         sprintf('With%sNo%s',g1Name,g2Name),sprintf('No%sWith%s',g1Name,g2Name),'Expected','Prob')
    result[1,1] = testLabel
    result[1,2] = scAll
    result[1,3] = scWithG1
    result[1,4] = scWithG2
    result[1,5] = scWithG1WithG2
    result[1,6] = scNoG1NoG2
    result[1,7] = scWithG1NoG2
    result[1,8] = scNoG1WithG2
    result[1,9] = expected
    result[1,10] = fetProb
    
    return (result)
  }

  return (fetProb)
}


load_cancer_types<-function(sourceFile,keepSubtypes=T,minSampleCount=10)
{
  sampleCancerTypes = read.csv(sourceFile)
  
  if(!keepSubtypes)
  {
    sampleCancerTypes = sampleCancerTypes %>% select(SampleId,CancerType)
  }
  
  if(minSampleCount>0)
  {
    minorTypes = sampleCancerTypes %>% group_by(CancerType) %>% count %>% filter(n<minSampleCount)
    sampleCancerTypes = sampleCancerTypes %>% mutate(CancerType=ifelse(CancerType %in% minorTypes$CancerType,'Other',as.character(CancerType)))
  }
  
  return (sampleCancerTypes)
}


##################
## MISC Analysis

rm(clusters)
rm(svData)

rm(svData)


svData = read.csv('~/logs/SVA_SVS.csv')
View(svData %>% filter(Type=='INF'))

svData = sv_load_and_prepare('~/logs/LNX_SVS.csv',F)
svData = sv_load_and_prepare('~/data/sv/LNX_SVS.csv')
svData = sv_load_and_prepare('~/data/sv/drivers/LNX_SVS.csv')
clusters = read.csv('~/data/sv/LNX_CLUSTERS.csv')
clusters = read.csv('~/data/sv/drivers/LNX_CLUSTERS.csv')

svData = read.csv('~/data/sv/drivers/LNX_SVS.csv')
clusters = read.csv('~/logs/LNX_CLUSTERS.csv')

svData = read.csv('~/data/sv/drivers/LNX_SVS.csv')
clusters = read.csv('~/logs/tmp/LNX_CLUSTERS.csv')
dmData = read.csv('~/logs/tmp/LNX_DOUBLE_MINUTES.csv')
View(dmData)

View(clusters %>% filter(ClusterId %in% c(61,25)))


ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)

ensemblGeneData38 = read.csv('~/data/sv/ensembl_hg38_gene_data.csv')
View(ensemblGeneData38)

ensemblTransExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'))


View(svData %>% filter(Type=='SGL'&ResolvedType=='DUP'))

View(svData %>% filter(SampleId=='WIDE01010667T') %>% filter(Type=='BND'|Type=='INF'))
View(svData %>% filter(SampleId=='CPCT02020691T'))
View(svData %>% filter(SampleId=='CPCT02020691T'&ClusterId==100) %>% 
       select(SampleId,Id,ClusterId,ClusterCount,ClusterDesc,LEStart,LEEnd,InsertSeq,RepeatClass,Annotations,everything()))

View(ensemblGeneData %>% filter(GeneId %in% c('ENSG00000265150','ENSG00000258486','ENSG00000202198','ENSG00000266037','ENSG00000263740','ENSG00000265735')))
View(ensemblGeneData %>% filter(Chromosome=='14'&GeneStart>50.05e6&GeneEnd<50.4e6))

lineElements = read.csv('~/data/line_elements.csv')
write.csv(lineElements %>% select(Chromosome,PosStart,PosEnd), '~/data/sv/line_elements.bed',row.names = F, quote = F)
fragileSites = read.csv('~/data/fragile_sites_hmf.csv')
write.csv(fragileSites %>% select(Chromosome,PosStart,PosEnd), '~/data/sv/fragile_sites.bed',row.names = F, quote = F)


View(head(sampleCancerTypes,100))


apcRna = read.csv('~/data/sv/fusion_like/apc_rna_read_counts.csv')
View(apcRna %>% group_by(Junction) %>% summarise(Samples=n(),
                                                 MinReads=min(ReadCount),
                                                 MedianReads=median(ReadCount),
                                                 MaxReads=max(ReadCount)))

bachRecords = read.csv('~/data/runs/WIDE01010111T/WIDE01010111T.bachelor.germline_variant.tsv',sep='\t')
View(bachRecords)
View(bachRecords %>% group_by(chromosome,position) %>% count)
rm(bachRecords)



View(fragileSites)
fragileSites = merge(fragileSites,ensemblGeneData %>% select(Gene=GeneName,GeneStart,GeneEnd),by='Gene',all.x=T)
ensemblGeneData38 = read.csv('~/data/sv/ensembl_hg38_gene_data.csv')
fragileSites = merge(fragileSites,ensemblGeneData38 %>% select(Gene=GeneName,GeneStartHg38=GeneStart,GeneEndHg38=GeneEnd),by='Gene',all.x=T)
write.csv(fragileSites %>% mutate(Width=GeneEndHg38-GeneStartHg38) %>% select(Chromosome,GeneStartHg38,GeneEndHg38,Width,Gene),
          '~/data/sv/fragile_sites_hmf_hg38.csv', quote = F, row.names = F)

# foldbacks by arm
fbStart = svData %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
fbEnd = svData %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
foldbacks = rbind(fbStart,fbEnd)
View(foldbacks)
foldbacks = foldbacks %>% mutate(IsChained=(OtherId!=Id),
                                 SingleBreakend=(OtherId==Id&FoldbackLength==0),
                                 FoldbackId=ifelse(Id<OtherId,Id,OtherId))


#########



# clustering reasons for long DEL-DUP vs long INV overlapping clusters
tmpSvData = read.csv('~/logs/LNX_SVS.csv')
clusters = read.csv('~/logs/LNX_CLUSTERS.csv')
nrow(clusters)

clusterHistory = read.csv('~/logs/LNX_CLUSTERING_HISTORY.csv')
View(clusterHistory)

View(clusters)
View(clusters %>% filter(grepl('LongDelDup',ClusterReasons)|grepl('LongInv',ClusterReasons)) %>% group_by(LongDelDup=grepl('LongDelDup',ClusterReasons),
                                                                                                          LongInv=grepl('LongInv',ClusterReasons)) %>% count)

View(clusters %>% filter(grepl('LongDelDup',ClusterReasons)))

View(tmpSvData %>% filter(grepl('LongDelDup',ClusterReason)&!grepl('LongInv',ClusterReason)) %>% mutate(Length=PosEnd-PosStart) %>%
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,Length,PosStart,PosEnd,ChrStart,Ploidy,ChainIndex,everything()))


View(ensemblGeneData %>% filter(GeneName %in% c('IGH','IGL','IGK')))
View(ensemblGeneData %>% filter(Chromosome==14&GeneStart>=106e6&GeneEnd<=106.4e6))

# DEBUG

View(ctCleanClusters %>% filter(SampleId=='CPCT02090058T'&ClusterId==156) %>% 
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,DeletedPerc,TotalDeleted,DeleteRatio,TotalDBLength,TotalRange,RangeRatio,ChainedLength,
              OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

View(svData %>% filter(SampleId=='CPCT02330080T',ChrStart==4))
View(svData %>% filter(SampleId=='CPCT02090058T',ChrStart==4))
View(view_cluster_arm_breakends(svData,'CPCT02090058T','',5,'Q'))



View(driverGenes)

# unique foldback data
View(foldbacks %>% group_by(SampleId,ClusterId,FoldbackId) %>% summarise(Chr=first(Chr),Arm=first(Arm),FoldbackLength=first(FoldbackLength)))


View(svData %>% filter(FoldbackLenStart==0|FoldbackLenEnd==0))
View(svData %>% filter(FoldbackLenStart==0) %>% group_by(SameSV=(Id==FoldbackLnkStart)) %>% count())
View(svData %>% filter(FoldbackLenEnd==0) %>% group_by(SameSV=(Id==FoldbackLnkEnd)) %>% count())

fbArmData = foldbacks %>% group_by(SampleId,Chr,Arm) %>% summarise(FoldbackCount=sum(ifelse(IsChained,0.5,1)))
fbArmData = fbArmData %>% ungroup() %>% mutate(FoldbackCount=ifelse(FoldbackCount<1,1,round(FoldbackCount)))



foldbackSummary = foldbacks %>% group_by(SampleId,Chr,Arm) %>% summarise(FoldbackCount=n()/2)



## ASSEMBLY anchor distances - should use min or max?

links = read.csv('~/data/sv/drivers/LNX_LINKS.csv')
# shortLinks = links %>% filter(IsAssembled=='true')
View(shortLinks)

newLinks = read.csv('~/data/sv/drivers/LNX_LINKS.csv')
shortLinks = newLinks %>% filter(TILength<=1000)

linkSvData = svData %>% filter(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd)|LnkLenStart<1000|LnkLenEnd<1000)
nrow(linkSvData)
View(linkSvData)
linkSvData = linkSvData %>% mutate(SvLength=ifelse(Type=='DUP'|Type=='DEL'|Type=='INV',PosEnd-PosStart,-1))
View(linkSvData)
View(linkSvData %>% select(Id,Type,SvLength))

shortLinks = merge(shortLinks,linkSvData %>% select(SampleId,Id,PosStart,AnchorSv1Start=AnchorStart,HomSv1Start=HomologyStart,Sv1Length1=SvLength),
                  by.x=c('SampleId','Id1','PosStart'),by.y=c('SampleId','Id','PosStart'),all.x=T)

shortLinks = merge(shortLinks,linkSvData %>% select(SampleId,Id,PosEnd,AnchorSv1End=AnchorEnd,HomSv1End=HomologyEnd,Sv1Length2=SvLength),
                  by.x=c('SampleId','Id1','PosStart'),by.y=c('SampleId','Id','PosEnd'),all.x=T)

shortLinks = merge(shortLinks,linkSvData %>% select(SampleId,Id,PosStart,AnchorSv2Start=AnchorStart,HomSv2Start=HomologyStart,Sv2Length1=SvLength),
                  by.x=c('SampleId','Id2','PosEnd'),by.y=c('SampleId','Id','PosStart'),all.x=T)

shortLinks = merge(shortLinks,linkSvData %>% select(SampleId,Id,PosEnd,AnchorSv2End=AnchorEnd,HomSv2End=HomologyEnd,Sv2Length2=SvLength),
                  by.x=c('SampleId','Id2','PosEnd'),by.y=c('SampleId','Id','PosEnd'),all.x=T)

nrow(shortLinks %>% filter(!is.na(AnchorSv1Start)))
nrow(shortLinks %>% filter(!is.na(AnchorSv1End)))
nrow(shortLinks %>% filter(!is.na(AnchorSv2Start)))
nrow(shortLinks %>% filter(!is.na(AnchorSv2End)))

View(shortLinks)

shortLinks = shortLinks %>% mutate(AnchorDistSv1=ifelse(!is.na(AnchorSv1Start),AnchorSv1Start,AnchorSv1End),
                                 AnchorDistSv2=ifelse(!is.na(AnchorSv2Start),AnchorSv2Start,AnchorSv2End),
                                 Sv1Length=ifelse(!is.na(Sv1Length1),Sv1Length1,Sv1Length2),
                                 Sv2Length=ifelse(!is.na(Sv2Length1),Sv2Length1,Sv2Length2),
                                 HomLengthSv1=ifelse(!is.na(AnchorSv1Start),stri_length(HomSv1Start),stri_length(HomSv1End)),
                                 HomLengthSv2=ifelse(!is.na(AnchorSv2Start),stri_length(HomSv2Start),stri_length(HomSv2End)),
                                 AdjAnchorDistSv1=AnchorDistSv1+HomLengthSv2,
                                 AdjAnchorDistSv2=AnchorDistSv2+HomLengthSv1,
                                 MinAnchorDist=pmin(AnchorDistSv1,AnchorDistSv2)-HomLengthSv2-HomLengthSv1,
                                 MaxAnchorDist=pmax(AnchorDistSv1,AnchorDistSv2),
                                 DiffMin=TILength-MinAnchorDist,
                                 DiffMax=TILength-MaxAnchorDist,
                                 MinHasShortSv=ifelse(AnchorDistSv1<AnchorDistSv2,Sv2Length>0&Sv2Length<1000,Sv1Length>0&Sv1Length<1000),
                                 MaxHasShortSv=ifelse(AnchorDistSv1>AnchorDistSv2,Sv2Length>0&Sv2Length<1000,Sv1Length>0&Sv1Length<1000))

shortLinks = shortLinks %>% filter(TILength>0)

asmbLinks = shortLinks %>% filter(IsAssembled=='true')

# View(asmbLinks %>% mutate(maxMin=MaxAnchorDist-MinAnchorDist,) %>% group_by(DiffMin>(10),DiffMax<(-10)) %>% count())
View(asmbLinks %>% filter(DiffMin<=0,DiffMin>-60) %>% group_by(DiffMin) %>% count())

View(asmbLinks %>% mutate(Min_LT_TI=DiffMin>(-10),Max_LT_TI=DiffMax>(-10)) %>% group_by(Min_LT_TI,Max_LT_TI) %>% count())
View(shortLinks %>% mutate(Min_LT_TI=DiffMin>(-10),Max_LT_TI=DiffMax>(-10)) %>% group_by(IsAssembled,Min_LT_TI,Max_LT_TI) %>% count())

View(ampDrivers %>% group_by(SampleId,Gene) %>% summarise(Clusters=n(),
                                                          DupCount=sum(ResolvedType=='DUP'),
                                                          SimpleGrpCount=sum(ResolvedType=='SIMPLE_GRP')) %>% filter(DupCount>1|SimpleGrpCount>=1))



View(head(asmbLinks,100))
View(asmbLinks %>% filter(is.na(AnchorDistSv1)|is.na(AnchorDistSv2)))
View(shortLinks %>% filter(TILength<30))

View(shortLinks %>% filter(is.na(DiffMin)))


######
## Double Minutes

dmProd = read.csv('~/data/sv/dm/LNX_DOUBLE_MINUTES_OLD.csv')
View(dmProd)
nrow(dmProd) # 603
nrow(dmProd %>% group_by(SampleId) %>% count) # 537
dmProdSamples = dmProd %>% group_by(SampleId) %>% count
write.csv(dmProd %>% group_by(SampleId) %>% count %>% select(SampleId),'~/logs/dm_sample_ids.csv',row.names = F,quote = F)

dmTest = read.csv('~/logs/LNX_DOUBLE_MINUTES.csv')
View(dmTest)

dmNew = read.csv('~/data/sv/dm/LNX_DOUBLE_MINUTES.csv')
View(dmNew)
nrow(dmNew) # 1032
nrow(dmNew %>% group_by(SampleId) %>% count) # 537

# recommended filters
View(dmNew %>% filter(ClosedSegLength>=5e3&OpenBreakends==0) %>% arrange(-IntExtMaxJcn))


# DMs with no open breakends
View(dmNew %>% filter(OpenBreakends==0) %>% arrange(-IntExtMaxJcn))


dmInDb = read.csv('~/logs/dm_sample_ids_in_db.csv')
View()
View(dmProdSamples %>% filter(!(SampleId %in% dmInDb$SampleId)))
write.csv(dmProdSamples %>% filter(!(SampleId %in% dmInDb$SampleId)))
# dm = merge(dm,dmDrivers %>% select(SampleId,ClusterId,HasDriver='DRIVER'),by=c('SampleId','ClusterId'),all.x=T)

dmClusters = read.csv('~/logs/LNX_CLUSTERS.csv')
View(dmClusters)

dmSvs = read.csv('~/data/sv/dm/LNX_SVS.csv')
View(dmSvs)
nrow(dmSvs)
View(dmSvs %>% filter(grepl('MAJOR',ClusterReason)))
View(dmSvs %>% filter(Annotations!=''))

View(dmSvs %>% filter(grepl('HIGH',ClusterReason)))

highJcnSVs = dmSvs %>% filter(grepl('HIGH',ClusterReason))
nrow(highJcnSVs)

oldDmSvs = read.csv('~/data/sv/dm/LNX_SVS_PROD.csv')

highJcnSVs = merge(highJcnSVs,oldDmSvs %>% select(SampleId,Id,OldClusterCount=ClusterCount,OldClusterReason=ClusterReason),by=c('SampleId','Id'),all.x=T)
View(highJcnSVs %>% mutate(Jcn=(JcnMin+JcnMax)*0.5) %>% select(SampleId,Id,ClusterCount,OldClusterCount,ClusterReason,OldClusterReason,ResolvedType,Jcn,everything()))

View(highJcnSVs %>% mutate(Jcn=(JcnMin+JcnMax)*0.5) %>% 
       filter(ClusterCount!=OldClusterCount) %>%
       select(SampleId,Id,ClusterCount,OldClusterCount,ClusterReason,OldClusterReason,ResolvedType,Jcn,everything()))

View(highJcnSVs %>% filter(ClusterCount!=OldClusterCount) %>% group_by(SampleId,ClusterId,ClusterCount,OldClusterCount) %>% count)

View(dmSvs %>% filter(Type=='INF'&NearestLen>=0) %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,Jcn,NearestLen,NearestType,DBLenStart,LnkSvStart,LnkLenStart,DMSV,everything()))

View(dmSvs %>% filter(Type=='INF') %>% 
       mutate(InDB=(DBLenStart==-2000),InTI=(LnkLenStart>=0),ZeroLengthTI=(LnkLenStart==0)) %>%
       group_by(InDB,InTI,ZeroLengthTI) %>% count)

View(dmSvs %>% filter(Type=='INF'&LnkLenStart>=0) %>% 
       group_by(TILength=ifelse(LnkLenStart==0,0,2**round(log(LnkLenStart,2)))) %>% count)

dmLinks = read.csv('~/data/sv/dm/LNX_LINKS.csv')
View(dmLinks)
zeroLenLinks = dmLinks %>% filter(TILength==0)
zeroLenLinks = merge(zeroLenLinks,dmSvs %>% select(SampleId,Id1=Id,Type1=Type),by=c('SampleId','Id1'),all.x=T)
zeroLenLinks = merge(zeroLenLinks,dmSvs %>% select(SampleId,Id2=Id,Type2=Type),by=c('SampleId','Id2'),all.x=T)

View(zeroLenLinks %>% mutate(OtherType=ifelse(Type=='INF',as.character(Type2),as.character(Type))) %>% group_by(OtherType) %>% count)
View(zeroLenLinks %>% mutate(OtherType=ifelse(Type=='INF',as.character(Type2),as.character(Type)),
                             ClusterSize=2**round(log(ClusterCount,2))) %>% group_by(OtherType,ClusterSize) %>% count %>% spread(OtherType,n,fill=0))


View(zeroLenLinks %>% select(SampleId,Id1,Id2,Type,Type2,ClusterId,ClusterCount,ResolvedType,TILength,PosStart,PosEnd,Jcn,everything()))
tmpLinks = read.csv('~/logs/LNX_LINKS.csv')
View(tmpLinks)
View(tmpLinks %>% mutate(Length=abs(PosStart-PosEnd)))


tmpSVs = read.csv('~/logs/LNX_SVS.csv')
View(tmpSVs)
View(tmpSVs %>% filter(Annotations!=''))

View(tmpSVs %>% filter(Type=='INF') %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,Jcn,NearestLen,NearestType,DBLenStart,LnkSvStart,LnkLenStart,everything()))


View(tmpSVs %>% group_by(SampleId,ChrStart,Pos=PosStart) %>% summarise(SVs=n(),
                                                                   InfCount=sum(Type=='INF'),
                                                                   SvId1=first(Id),
                                                                   SvId2=last(Id),
                                                                   Type1=first(Type),
                                                                   Type2=last(Type),
                                                                   JcnFirst=first((JcnMin+JcnMax)*0.5),
                                                                   JcnLast=last((JcnMin+JcnMax)*0.5),
                                                                   AnchorFirst=first(AnchorStart),
                                                                   AnchorLast=last(AnchorStart),
                                                                   ) %>% filter(SVs>1&InfCount==1))

View(tmpSVs %>% group_by(SampleId,ChrStart,Pos=ifelse(Type=='INF',PosStart,PosEnd)) %>% summarise(SVs=n(),
                                                                   InfCount=sum(Type=='INF'),
                                                                   SvId1=first(Id),
                                                                   SvId2=last(Id),
                                                                   Type1=first(Type),
                                                                   Type2=last(Type),
                                                                   JcnFirst=first((JcnMin+JcnMax)*0.5),
                                                                   JcnLast=last((JcnMin+JcnMax)*0.5),
                                                                   AnchorFirst=first(ifelse(Type=='INF',AnchorStart,AnchorEnd)),
                                                                   AnchorLast=last(ifelse(Type=='INF',AnchorStart,AnchorEnd))) %>% filter(SVs>1&InfCount==1))


gripssData = read.csv('~/logs/LNX_SVS.csv')
View(gripssData)

View(gripssData %>% filter(grepl('trs',AsmbStart)|grepl('trs',AsmbEnd)) %>% select(SampleId,Id,AsmbStart,AsmbEnd,LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd,everything()))



dmDrivers = read.csv('~/logs/LNX_DRIVERS.csv')
dmDrivers = dmDrivers %>% filter(Gene=='CDK4'|Gene=='MDM2')
View(dmDrivers)
write.csv(dmDrivers %>% filter(ClusterId>=0) %>% group_by(SampleId,ClusterId) %>% count %>% select(SampleId,ClusterId),
          '~/logs/cdk4_mdm2_sample_clusters.csv',row.names = F,quote = F)

dmDrivers = merge(dmDrivers,dmClusters %>% select(SampleId,ClusterId,Foldbacks,Annotations),by=c('SampleId','ClusterId'),all.x=T)
dmDrivers = merge(dmDrivers,dm %>% select(SampleId,ClusterId,DMSvCount,DMSvTypes,Chromosomes,MaxCopyNumber,MaxPloidy),by=c('SampleId','ClusterId'),all.x=T)
View(dmDrivers)

write.csv(dmDrivers,'~/logs/cdk4_mdm2_driver_data.csv',row.names = F,quote = F)




# short TIs which could instead be DBs
shortInfs = newLinks %>% filter(TILength>0&TILength<=1000&IsAssembled=='false'&LocTopTypeStart=='TI_ONLY'&LocTopTypeEnd=='TI_ONLY')
View(shortInfs)

shortInfs = merge(shortInfs,linkSvData %>% select(SampleId,Id,PosStart,CNSv1Start=CNStart,CNChgSv1Start=CNChgStart,Sv1Ploidy1=Ploidy),
                   by.x=c('SampleId','Id1','PosStart'),by.y=c('SampleId','Id','PosStart'),all.x=T)

shortInfs = merge(shortInfs,linkSvData %>% select(SampleId,Id,PosEnd,CNSv1End=CNEnd,CNChgSv1End=CNChgEnd,Sv1Ploidy2=Ploidy),
                   by.x=c('SampleId','Id1','PosStart'),by.y=c('SampleId','Id','PosEnd'),all.x=T)

shortInfs = merge(shortInfs,linkSvData %>% select(SampleId,Id,PosStart,CNSv2Start=CNStart,CNChgSv2Start=CNChgStart,Sv2Ploidy1=Ploidy),
                   by.x=c('SampleId','Id2','PosEnd'),by.y=c('SampleId','Id','PosStart'),all.x=T)

shortInfs = merge(shortInfs,linkSvData %>% select(SampleId,Id,PosEnd,CNSv2End=CNEnd,CNChgSv2End=CNChgEnd,Sv2Ploidy2=Ploidy),
                   by.x=c('SampleId','Id2','PosEnd'),by.y=c('SampleId','Id','PosEnd'),all.x=T)

shortInfs = shortInfs %>% mutate(CNSv1=ifelse(!is.na(CNSv1Start),CNSv1Start,CNSv1End),
                                 CNChgSv1=ifelse(!is.na(CNChgSv1Start),CNChgSv1Start,CNChgSv1End),
                                 PloidySv1=ifelse(!is.na(Sv1Ploidy1),Sv1Ploidy1,Sv1Ploidy2),
                                 CNSv2=ifelse(!is.na(CNSv2Start),CNSv2Start,CNSv2End),
                                 CNChgSv2=ifelse(!is.na(CNChgSv2Start),CNChgSv2Start,CNChgSv2End),
                                 PloidySv2=ifelse(!is.na(Sv2Ploidy1),Sv2Ploidy1,Sv2Ploidy2))

View(shortInfs)

shortInfs = shortInfs %>% mutate(CnChgDiff=abs(CNChgSv1-CNChgSv2),
                                 PloidyDiff=abs(PloidySv1-PloidySv2),
                                 CnChgDiffBucket=2**round(log(CnChgDiff,2)),
                                 PloidyBucket=2**round(log(Ploidy,2)),
                                 TILengthBucket=2**round(log(TILength,2)))

View(shortInfs %>% filter(CNSv1>=1&CNSv2>=1) %>% group_by(CnChgDiffBucket) %>% count())
View(shortInfs %>% filter(CNSv1>=1&CNSv2>=1) %>% group_by(CnChgDiffBucket,LinkReason) %>% count() %>% spread(LinkReason,n))
View(shortInfs %>% filter(CNSv1>=1&CNSv2>=1&LinkReason=='ONLY') %>% group_by(CnChgDiffBucket,TILengthBucket) %>% count() %>% spread(TILengthBucket,n))


View(shortInfs %>% group_by(CnChgDiffBucket=2**round(log(CnChgDiff/pmax(CNSv1,CNSv2),2))) %>% count())

print(ggplot(data = asmbLinks, aes(x=TILength))
      + geom_point(aes(y=MaxAnchorDist,colour='MaxAnchorDist'))
      + geom_point(aes(y=MinAnchorDist,colour='MinAnchorDist'))
      + labs(title = "TI length sv Anchor distances"))

print(ggplot(data = shortLinks, aes(x=TILength))
      + geom_point(aes(y=DiffMin,colour='MinAnchorDist'))
      + facet_wrap(~IsAssembled)
      + labs(title = "TI length sv Anchor distances")
      + ylim(-500,500))

print(ggplot(data = asmbLinks, aes(x=TILength))
      + geom_point(aes(y=DiffMax,colour='MinAnchorDist'))
      + labs(title = "TI length sv Anchor distances")
      + facet_wrap(~MaxHasShortSv)
      + ylim(-500,500))

print(ggplot(data = shortLinks, aes(x=TILength))
      + geom_point(aes(y=DiffMax,colour='MaxAnchorDist'))
      + facet_wrap(~IsAssembled)
      + labs(title = "TI length sv Anchor distances")
      + ylim(-500,500))

print(ggplot(data = asmbLinks %>% mutate(maxMin=MaxAnchorDist-MinAnchorDist), aes(x=TILength))
      + geom_point(aes(y=maxMin,colour='maxMin'))
      + labs(title = "TI length sv Anchor distances")
      + ylim(-500,500))


rm(dmSvs)
rm(matchedSVs)
rm(mbClusterData)
rm(mbAmpSVs)
rm(mbSvCombined)
rm(mbSvData)
rm(mbSVs)



# SIMPLE GROUPS

View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% group_by(HasAssembledLinks=AssemblyTIs>0,
                                                                  FullyAssembled=(ClusterCount==AssemblyTIs+1)) %>% count())

View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% group_by(ClusterDesc,HasAssembledLinks=AssemblyTIs>0,
                                                                  FullyAssembled=(ClusterCount==AssemblyTIs+1)) %>% count()
     %>% spread(ClusterDesc,n))


View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% filter(ClusterDesc=='DEL=2'&AssemblyTIs==0))
View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% group_by(ClusterDesc,AssemblyTIs==ClusterCount-1) %>% count())
View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% filter((grep('DEL=2',ClusterDesc)|grep('DEL=3',ClusterDesc))&AssemblyTIs==0))










#######
# Neo-Epitopes

neoEpitopes = read.csv('~/logs/SVA_NEO_EPITOPES.csv')
neoEpitopes = read.csv('~/data/sv/fusions/SVA_NEO_EPITOPES.csv')
nrow(neoEpitopes)
View(neoEpitopes)
View(neoEpitopes %>% filter(PhaseMatched=='true'&PhaseUp==PhaseDown))
View(neoEpitopes %>% filter(PhaseMatched=='false'))
View(neoEpitopes %>% filter(PhaseMatched=='true'&PhaseUp==PhaseDown&PhaseUp==0))
View(neoEpitopes %>% filter(grepl('_',UpstreamAminoAcids)))
View(neoEpitopes %>% filter(grepl('_',NovelAminoAcid)))
View(neoEpitopes %>% filter(grepl('_',DownstreamAminoAcids)))
View(neoEpitopes %>% group_by(SameGene) %>% count())
View(neoEpitopes %>% group_by(RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown) %>% count())
View(neoEpitopes %>% filter(PhaseUp==-1))
View(neoEpitopes %>% filter(PhaseDown==-1))

View(neoEpitopes %>% group_by(SampleId,SvIdUp,SvIdDown) %>% count() %>% filter(n>1) %>% arrange(-n))

neoEpitopes = neoEpitopes %>% mutate(NovelAALength=stri_length(NovelAminoAcid))

View(neoEpitopes %>% group_by(PhaseUp,PhaseDown) %>% count())
View(neoEpitopes %>% group_by(PhaseDown) %>% count())

baseAminoAcids = read.csv('~/data/base_amino_acid.csv')
View(baseAminoAcids)
write.csv(baseAminoAcids %>% select(AminoAcid,Bases), '~/data/base_amino_acid2.csv', row.names=F, quote=F)


######
# Clonal / Subclonal cluster merges
clonalMerges = read.csv('~/data/sv/clonal_discrep_clustering.csv')
View(clonalMerges)
View(clonalMerges %>% group_by(ClusterReason) %>% count())
write.csv(clonalMerges %>% select(-Time),'~/data/sv/subclonal_cluster_merges.csv', row.names = F, quote = F)



######
# Clonal / Subclonal cluster merges
clonalMerges = read.csv('~/data/sv/clonal_discrep_clustering.csv')
View(clonalMerges)
View(clonalMerges %>% group_by(ClusterReason) %>% count())
write.csv(clonalMerges %>% select(-Time),'~/data/sv/subclonal_cluster_merges.csv', row.names = F, quote = F)


clusterHistory = read.csv('~/data/sv/clustering_history.csv')
View(clusterHistory)
View(head(clusterHistory,100))
View(clusterHistory %>% filter(ClonalDiscrepancy=='true'))
View(clusterHistory %>% filter(ClonalDiscrepancy=='true') %>% group_by(Reason) %>% count())
View(clusterHistory %>% group_by(ClusterReason) %>% count())

svData = read.csv('~/logs/SVA_SVS.csv')

View(svData %>% filter(ClusterId==459))

View(svData %>% filter(ClusterId==524) %>% select(SampleId,Id,ChrStart,PosStart,PosEnd,OrientStart,OrientEnd,AnchorStart,AnchorEnd,Homology))


clusters = read.csv('~/data/sv/fusions/LNX_CLUSTERS.csv')
clusters = read.csv('~/logs/LNX_CLUSTERS.csv')

View(clusters %>% filter(ResolvedType=='UNBAL_TRANS_TI'))
View(clusters %>% filter(ResolvedType=='UNBAL_TRANS_TI') %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ArmCount,TotalTIs,ShortTIs,OriginArms,FragmentArms))

View(clusters %>% filter(ResolvedType=='UNBAL_TRANS_TI') %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ArmCount,TotalTIs,ShortTIs,OriginArms,FragmentArms)
     %>% filter(TotalTIs>ShortTIs+1))


########
# Spanning dup breakends

dupBreakends = read.csv('~/logs/dup_breakends.csv')
View(dupBreakends %>% filter(!grepl('eqv',AsmbStart)&!grepl('eqv',AsmbEnd)&stri_length(InsertSeq)>=30&Type!='SGL'))
View(dupBreakends %>% group_by(IsEqv=grepl('eqv',AsmbStart)|grepl('eqv',AsmbEnd),HasInsertSeq=stri_length(InsertSeq)>=30,Type) %>% count())


#######
# Pseudogenes

links = read.csv('~/logs/SVA_LINKS.csv')
pgLinks = links %>% filter(ExonMatch!='')
View(pgLinks)

pgLinks = pgLinks %>% 
  separate(ExonMatch,c('Gene','TransId','ExonRank','ExonLength','StartHomOffset','StartPosOffset','HomMismatchStart','EndHomOffset','EndPosOffset','HomMismatchEnd'),sep=';')

pgLinks = within(pgLinks,rm(HomMismatchStart))
pgLinks = within(pgLinks,rm(HomMismatchEnd))

pgLinks = pgLinks %>% mutate(StartPosOffset=as.numeric(as.character(StartPosOffset)),
                             EndPosOffset=as.numeric(as.character(EndPosOffset)),
                             StartHomOffset=as.numeric(as.character(StartHomOffset)),
                             EndHomOffset=as.numeric(as.character(EndHomOffset)),)

pgLinks = pgLinks %>% mutate(ExonTruncStart=StartPosOffset<0,
                             ExonTruncEnd=EndPosOffset<0,
                             ExonExtendedStart=StartPosOffset>0,
                             ExonExtendedEnd=EndPosOffset>0)

write.csv(pgLinks,'~/logs/pseudo_gene_links.csv', row.names = F, quote = F)

View(pgLinks)

View(pgLinks %>% filter(SampleId=='CPCT02010976T'&ClusterId==242))


View(pgLinks %>% group_by(StartHomOffset,EndHomOffset) %>% count())
View(pgLinks %>% group_by(HomMismatchStart,HomMismatchEnd) %>% count())
View(pgLinks %>% group_by(ExonTruncStart,ExonTruncEnd,ExonExtendedStart,ExonExtendedEnd) %>% count())

pgLinksStart = pgLinks %>% select(SampleId,SvId=Id1,HomOffset=StartHomOffset,PosOffset=StartPosOffset)
pgLinksEnd = pgLinks %>% select(SampleId,SvId=Id2,HomOffset=EndHomOffset,PosOffset=EndPosOffset)
pgLinksCombined = rbind(pgLinksStart,pgLinksEnd)
View(pgLinksStart)

homOffsetMatches = pgLinksCombined %>% group_by(SampleId,SvId) %>% 
  summarise(HomOffsetUsed=sum(PosOffset==0),HomOffset1=first(HomOffset),HomOffset2=last(HomOffset))


homOffsetMatches = homOffsetMatches %>% filter(HomOffsetUsed==2) %>% mutate(HomOffMatched=HomOffset1==HomOffset2)
View(homOffsetMatches)

View(homOffsetMatches %>% group_by(HomOffMatched) %>% count())


##########
# Shards vs recip inversions
tiLinks = read.csv('~/logs/SVA_LINKS.csv')
tiLinks = read.csv('~/data/sv/drivers/SVA_LINKS.csv')

clusters = read.csv('~/data/sv/drivers/SVA_CLUSTERS.csv')
clusters = read.csv('~/logs/SVA_CLUSTERS.csv')

pairClusters = clusters %>% filter(ClusterCount==2&(ResolvedType=='DUP'|ResolvedType=='DEL'|ResolvedType=='RECIP_INV'))
nrow(pairClusters)
View(pairClusters)
pairClusters = pairClusters %>% separate(Annotations,c('SyntheticLength','TiLength','SvGapLength'),sep=';')

View(tiLinks %>% group_by(ResolvedType) %>% count())

pairLinks = tiLinks %>% filter(ClusterCount==2&(ResolvedType=='DUP'|ResolvedType=='DEL'|ResolvedType=='RECIP_INV'))
View(pairLinks)
nrow(pairLinks)

pairLinks = merge(pairLinks,pairClusters %>% select(SampleId,ClusterId,SyntheticLength,TiLength,SvGapLength),
                  by=c('SampleId','ClusterId'),all.x=T)

pairLinks$SyntheticLength = as.numeric(pairLinks$SyntheticLength)

pairLinks = pairLinks %>% mutate(TiLengthBucket=ifelse(TILength>0,2**round(log(TILength,2)),0),
                                 SynLengthBucket=ifelse(SyntheticLength>0,2**round(log(SyntheticLength,2)),0))

pairLinks = pairLinks %>% filter(ResolvedType!='RECIP_INV'|LocationType=='Internal')
View(pairLinks %>% filter(ResolvedType=='RECIP_INV'&LocationType!='Internal'))

View(pairLinks %>% filter(LocationType!='Remote'&(ClusterDesc=='INV=2'|ClusterDesc=='DEL=2')) %>% 
       group_by(ResolvedType,ClusterDesc,LocationType,TiLengthBucket) %>% count())

View(pairLinks %>% filter(LocationType!='Remote') %>% group_by(ResolvedType,ClusterDesc,LocationType) %>% count())

pairLinks = pairLinks %>% mutate(Category=paste(ResolvedType,ClusterDesc,LocationType,sep=' '))

print(plot_length_facetted(pairLinks,"TiLengthBucket<1e4&(ClusterDesc=='INV=2'|ClusterDesc=='DEL=2')",'TiLengthBucket,Category','TiLengthBucket','Category','DELs vs RECIP-INVs'))

print(plot_length_facetted(pairLinks,"SynLengthBucket<1e6&(ClusterDesc=='INV=2'|ClusterDesc=='DEL=2')",'SynLengthBucket,Category','SynLengthBucket','Category','DELs vs RECIP-INVs'))

print(plot_length_facetted(pairLinks,"SynLengthBucket<1e6&ResolvedType='DEL'",'SynLengthBucket,Category','SynLengthBucket','Category','DELs vs RECIP-INVs'))

# plot_length_facetted<-function(data, filterStr, groupByStr, bucketStr, facetStr, titleStr, logScaleX=T, logScaleY=F)
  
print(plot_length_facetted(pairLinks,"SynLengthBucket<1e9&ResolvedType=='DUP'",'SynLengthBucket,Category','SynLengthBucket','Category','Synthetic DUP lengths'))



################
# Pair Resolution

clusters = read.csv('~/logs/SVA_CLUSTERS.csv')
clusters = read.csv('~/data/sv/drivers/SVA_CLUSTERS.csv')
View(clusters %>% group_by(ResolvedType,Synthetic) %>% count() %>% spread(Synthetic,n))
#View(clusters %>% filter(ResolvedType=='DUP_TI'))

View(clusters %>% filter(ClusterCount>1&ResolvedType!='COMPLEX'&ResolvedType!='LINE'))


# length characteristics
pairClusters = clusters %>% filter(Annotations!='') %>%
  filter((ResolvedType=='DEL'&Synthetic=='true')|(ResolvedType=='DUP'&Synthetic=='true')
         |ResolvedType=='DUP_TI'|ResolvedType=='DEL_TI'|ResolvedType=='RECIP_TRANS'|ResolvedType=='RECIP_INV'
         |ResolvedType=='RECIP_TRANS_DUPS'|ResolvedType=='RECIP_TRANS_DEL_DUP'
         |ResolvedType=='RECIP_INV_DUPS'|ResolvedType=='RECIP_INV_DEL_DUP')

View(pairClusters)
View(pairClusters %>% filter(Annotations==''))

pairClusters = pairClusters %>% separate('Annotations',c('SynLength1','SynLength2','SynGapLength'),sep=';')
View(pairClusters)

pairClusters = pairClusters %>% mutate(SynLength1=as.numeric(as.character(SynLength1)),
                                       SynLength2=as.numeric(as.character(SynLength2)),
                                       SynGapLength=as.numeric(as.character(SynGapLength)))

View(pairClusters %>% filter(ResolvedType=='RECIP_INV'))




# AMP gene chaining
View(view_cluster_arm_breakends(svData,'CPCT02220026T','','19','P'))

View(view_cluster_arm_breakends(svData,'CPCT02010623T','','10',''))
View(view_cluster_arm_breakends(svData,'CPCT02140050T','','8','Q'))

clusteringHistory = read.csv('~/data/sv/SVA_CLUSTER_HISTORY.csv')
View(clusteringHistory %>% filter(SampleId=='CPCT02010692TII'))

clusteringHistory = clusteringHistory %>% mutate(MinClusterCount=ifelse(ClusterCount1<ClusterCount2,ClusterCount1,ClusterCount2))
View(clusteringHistory)
View(clusteringHistory %>% filter(SampleId=='CPCT02010933T'))

View(svData %>% filter(SampleId=='CPCT02040160T'&ClusterId==12) %>% select(SampleId,Id,ChrStart,PosStart,PosEnd,Type,Ploidy,AsmbStart,AsmbEnd))

View(subsetSvData %>% filter(SampleId=='CPCT02010509T'&ClusterId==135))


view_cluster_sv('CPCT02120069T',285)
View(subsetSvData %>% filter(SampleId=='CPCT02120069T'&ClusterId==285))



# collapsing of INFs
inferredSegments = read.csv('~/data/sv/inferred_segments.csv')
inferredSegments = inferredSegments %>% select(-Time)
View(inferredSegments)
View(inferredSegments %>% group_by(InfCount) %>% count())



infData = svData %>% filter(Type=='INF')

View(infData %>% group_by(SampleId,ChrStart,ArmStart) %>% count() %>% filter(n>10) %>% arrange(-n))


# chaining investigations
tmpLinks = read.csv('~/logs/SVA_LINKS.csv')
View(tmpLinks)

View(tmpLinks %>% filter(SampleId=='DRUP01070036T'&ClusterId==46) 
     %>% select(Id1,Id2,PosStart,PosEnd,LinkIndex,LinkReason,Ploidy,PloidyUncertainty,ChainId,ChainIndex,ChainCount,everything()))


subsetSvData = read.csv('~/logs/SVA_SVS.csv')
View(subsetSvData %>% filter(ResolvedType=='LOW_VAF'))

View(subsetSvData %>% filter(SampleId=='DRUP01070036T'&ClusterId==46))

clusterArmData = view_cluster_arm_breakends(subsetSvData,'DRUP01070036T',46,'X','Q')
clusterArmData = merge(clusterArmData,tmpLinks %>%)


View(subsetSvData %>% filter(SampleId=='CPCT02220026T'&(ChrStart==5|ChrEnd==5)))



View(svData)


View(driverClusters %>% filter(ArmClusterCount==AcTIOnly))

# Compare DM lists
prodDMs = read.csv('~/data/sv/SVA_DM_PROD.csv')
ppDMs = read.csv('~/data/sv/DM_REF.csv')
nrow(ppDMs)
nrow(prodDMs)

View(prodDMs)
prodDMs = prodDMs %>% mutate(SampleClusterId = paste(SampleId,ClusterId,sep='_'))
ppDMs = ppDMs %>% mutate(SampleClusterId = paste(SampleId,ClusterId,sep='_'))
View(ppDMs %>% filter(!(SampleClusterId %in% prodDMs$SampleClusterId)))


links = read.csv('~/logs/SVA_LINKS.csv')
subsetClusters = read.csv('~/logs/SVA_CLUSTERS.csv')
View(subsetClusters)
View(links %>% filter(OverlapCount>10))

View(links %>% filter(SampleId=='CPCT02130013T'&ClusterId==56))
View(links %>% filter(SampleId=='CPCT02011051T'&ClusterId==235&(Id1==202|Id2==202)))

View(links %>% filter(SampleId=='CPCT02011051T'&ClusterId==235) %>% 
       filter(ChrArm=='3_Q'&TILength>1e3&PosStart>=136729497&PosEnd<=139079969))


View(svData %>% filter((Type=='SGL'|Type=='INF')&(RepeatClass!=''|RepeatType!='')) %>% group_by(Type,RepeatClass,RepeatType) %>% count())


subsetSvData = read.csv('~/logs/SVA_SVS.csv')
View(subsetSvData %>% group_by(ResolvedType) %>% count())
View(subsetSvData %>% filter(grepl('Sate',ClusterReason)))



View(links %>% filter(ClusterCount<=50) %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
       summarise(Links=n(),
                 MaxOverlaps=max(OverlapCount),
                 MedOverlaps=median(OverlapCount)))

View(links %>% filter(ClusterCount>=5&ClusterCount<=50) %>% group_by(LinkReason) %>% count())




totalClusters = nrow(clusters)
totaSVs = sum(clusters$ClusterCount)
View(clusters %>% group_by(SuperType,ResolvedType,Synthetic) %>% 
       summarise(Clusters=n(),
                 Percent=round(n()/totalClusters,3),
                 HasSGLs=sum(SglCount>0|InfCount>0),
                 SvCount=sum(ClusterCount),
                 SvPercent=round(sum(ClusterCount)/totaSVs,3)))

View(clusters %>% filter(ResolvedType=='SimpleSV') %>% group_by(ClusterCount) %>% count())
View(clusters %>% filter(ClusterCount>1&ResolvedType=='SimpleSV') %>% group_by(ClusterDesc) %>% count())

View(clusters %>% group_by(SuperType,ResolvedType) %>% count())
View(clusters %>% group_by(SuperType,ResolvedType,Synthetic) %>% count())
View(clusters %>% filter(ResolvedType=='SIMPLE_GRP') %>% group_by(ClusterDesc) %>% count())
View(clusters %>% filter(ResolvedType=='DUP_BE') %>% group_by(ClusterDesc) %>% count())


## Short DUPs vs overlapping DBs
shortDups = svData %>% filter(ResolvedType=='DUP'&ClusterCount==1)
shortDups = shortDups %>% mutate(Length=PosEnd-PosStart,LengthBucket=5*round(Length/5)) %>% filter(Length<1e3)
nrow(shortDups)
View(shortDups)

print(plot_length_facetted(shortDups,"",'LengthBucket','LengthBucket','','Simple Short DUPS'))


knownFusions = read.csv('~/data/knownFusionPairs.csv')
promFive = read.csv('~/data/knownPromiscuousFive.csv')
promThree = read.csv('~/data/knownPromiscuousThree.csv')

View(knownFusions%>% filter(!(fiveGene %in% ensemblGeneData$GeneName)))
View(knownFusions%>% filter(!(threeGene %in% ensemblGeneData$GeneName)))
View(promFive%>% filter(!(gene %in% ensemblGeneData$GeneName)))
View(promThree%>% filter(!(gene %in% ensemblGeneData$GeneName)))

tiLinks = read.csv('~/data/sv/fusions/SVA_LINKS.csv')


# region enrichment for short TIs
shortTIs = tiLinks %>% filter(LocTopTypeStart=='TI_ONLY'&LocTopTypeEnd=='TI_ONLY')
View(shortTIs)
shortTIs = shortTIs %>% separate(ChrArm,c('Chromosome','Arm'),sep='_')
shortTIsLocations = shortTIs %>% group_by(ResolvedType,Chromosome,Arm,Location=baseWindow*round(PosStart/baseWindow)) %>% count()
View(shortTIsLocations)

nrow(shortTIs %>% filter(PosEnd-PosEnd<=5e3))

shortTIs = shortTIs %>% filter(ResolvedType!='LINE')
baseWindow=1e6
windowCount=3e9/baseWindow
totalCount=nrow(shortTIs)
expected=totalCount/windowCount

shortTIsLocations2 = shortTIs %>% filter(ResolvedType!='LINE') %>% group_by(Chromosome,Location=baseWindow*round(PosStart/baseWindow)) %>% 
  summarise(Count=n(),
            PoisProb=1-ppois(Count-1,expected)) %>% arrange(-Count)

rowIndex = data.frame(as.numeric(as.character(rownames(shortTIsLocations2))))
colnames(rowIndex) <- c("RowIndex")
View(rowIndex)
tmp = cbind(shortTIsLocations2,rowIndex)
View(tmp)
View(shortTIsLocations2)




nrow(svData)

eqvSvs = svData %>% filter(grepl('eqv', AsmbStart)|grepl('eqv', AsmbEnd))
eqvSvs = eqvSvs %>% filter(ResolvedType!='DUP_BE')
View(eqvSvs)

View(svData %>% filter(grepl('eqv', AsmbStart)|grepl('eqv', AsmbEnd)) %>% group_by(Type) %>% count())


candDupSvs = svData %>% filter(Type!='SGL'&Type!='INF'&ClusterCount>1&grepl('Prox',ClusterReason))
nrow(candDupSvs)
candDupSvs = candDupSvs %>% mutate(PS=35*round(PosStart/35),PE=35*round(PosEnd/35))

nonSglEqvs = candDupSvs %>% group_by(SampleId,ClusterId,ClusterCount,ResolvedType,ClusterDesc,Type,PS,PE,OrientStart,OrientEnd) %>% 
  summarise(Count=n(),
            Id1=first(Id),
            Id2=last(Id),
            PosDiffStart=first(PosStart)-last(PosStart),
            PosDiffEnd=first(PosEnd)-last(PosEnd),
            AsmCount=sum(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd)),
            Eqv=sum(grepl('eqv',AsmbStart)|grepl('eqv',AsmbEnd)),
            AsmbStart1=first(AsmbStart),
            AsmbStart2=last(AsmbStart),
            AsmbEnd1=first(AsmbEnd),
            AsmbEnd2=last(AsmbEnd)) %>% filter(Count>1)

nrow(nonSglEqvs)
nrow(nonSglEqvs %>% filter(grepl('eqv',AsmbStart1)|grepl('eqv',AsmbStart2)|grepl('eqv',AsmbEnd1)|grepl('eqv',AsmbEnd2)))
View(nonSglEqvs %>% filter(grepl('eqv',AsmbStart1)|grepl('eqv',AsmbStart2)|grepl('eqv',AsmbEnd1)|grepl('eqv',AsmbEnd2)))
View(nonSglEqvs)
View(nonSglEqvs %>% filter(abs(PosDiffStart)<=1&abs(PosDiffEnd)<=1))

View(svData %>% filter(Id %in% c(15655106,15655107,15655113)) %>% 
       select(SampleId,Id,ClusterId,ClusterCount,Type,ChrStart,PosStart,OrientStart,PosEnd,OrientEnd,CNChgStart,CNChgEnd,Ploidy,PloidyMin,PloidyMax,AsmbStart,AsmbEnd))


svData$ResolvedType = ifelse(svData$ResolvedType=='NONE','COMPLEX',as.character(svData$ResolvedType))
View(clusters %>% group_by(ResolvedType) %>% count())

view_sv_ids(c(15228417,15228418))

view_sv_ids(c(15228259,15228306))

view_sv_ids(c(14916689,14916690))


View(view_cluster_arm_breakends(specificSample, 'CPCT02010347TII', 68, 1, 'P'))


specificSample = read.csv('~/logs/CPCT02010347TII_SVA.csv')

View(specificSample %>% filter(ClusterId==68) %>% select(Id,Type,ChrStart,PosStart,ChrEnd,PosEnd,ClusterReason,Ploidy,PloidyMin,PloidyMax,CNChgStart,CNChgEnd))

View(specificSample %>% filter(ClusterId==68) %>% select(Id,Type,ChrStart,PosStart,ChrEnd,PosEnd,ClusterReason,Ploidy,PloidyMin,PloidyMax,CNChgStart,CNChgEnd))

# synthetic clusters and local topology
view_chromosome_sv('CPCT02330133T',3,3)

view_clusters_sv('CPCT02010347TII', 66)

view_sv_ids(c(15228255,15228256,15228279,15228265))

View(svData %>% filter(grepl('asm16-418812/480',AsmbStart)))

View(svData %>% filter(Id %in% c(15228255,15228256,15228279,15228265)) %>% 
       select(SampleId,Id,ClusterId,ClusterCount,Type,CNChgStart,CNChgEnd,Ploidy,PloidyMin,PloidyMax,AsmbStart,AsmbEnd))


View(svData %>% filter(Id %in% c(15655106,15655107,15655113)) %>% 
       select(SampleId,Id,ClusterId,ClusterCount,Type,PosStart,OrientStart,PosEnd,OrientEnd,CNChgStart,CNChgEnd,Ploidy,PloidyMin,PloidyMax,AsmbStart,AsmbEnd))

view_clusters_sv('CPCT02080019T', 944)




# LINE clusters without poly A-T support
View(svData %>% filter(IsLINE) %>% group_by(LEStart,LEEnd) %>% count())

rm(lineAll)

lineSvData = svData %>% filter(ResolvedType=='LINE')
rm(lineSvData)
nrow(lineSvData)

lineClusterData = lineSvData %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
  summarise(AnyLine=sum(grepl('Known',LEStart)|grepl('Suspect',LEStart)),
            Known=sum(grepl('Known',LEStart)),
            Suspect=sum(grepl('Suspect',LEStart)),
            PolyAT=sum(grepl('AAAAAAAAAAA',InsertSeq)|grepl('TTTTTTTTTTT',InsertSeq)),
            BndCount=sum(Type=='BND'),
            SglCount=sum(Type=='SGL')) %>% mutate(SampleClusterId=paste(SampleId, ClusterId, sep='_'))

clusters = clusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId, sep='_'))
lineClusterData = merge(lineClusterData,clusters %>% select(SampleClusterId,ArmCount),by='SampleClusterId',all.x=T)

View(lineClusterData)
write.csv(lineClusterData, 'logs/sv_line_clusters.csv', row.names = F, quote = F)
rm(lineClusterData)

print(ggplot(data = lineClusterData %>% filter(Suspect>0), aes(x=ArmCount))
      + geom_point(aes(y=Suspect, colour='Suspect'))
      + geom_point(aes(y=Known, colour='Known'))
      + geom_point(aes(y=PolyAT, colour='PolyAT'))
      + labs(title = "ArmCount vs Suspect/Known LINE Svs"))

print(ggplot(data = lineClusterData %>% filter(Suspect>0), aes(x=ArmCount,y=ClusterCount))
      + geom_point()
      + labs(title = "ArmCount vs Cluster Svs"))

print(ggplot(data = lineClusterData, aes(x=ArmCount, y=Suspect))
      + geom_point()
      + labs(title = "ArmCount vs Suspect LINE Svs"))

print(ggplot(data = lineClusterData %>% filter(SuspectCount>0), aes(x=ArmCount,y=ClusterCount))
      + geom_point()
      + labs(title = "ArmCount vs Cluster Svs"))



#######
# Subclonal SVs within larger clusters with higher ploidy

svData = read.csv('~/data/sv/fusions/SVA_SVS.csv')

csvData2 = svData %>% filter(ClusterCount>=3&ResolvedType!='LINE')

ploidyData2 = csvData2 %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
  summarise(Ploidy_LT_0_25=round(sum(PloidyMax<0.25)/n(),3),
            Ploidy_LT_0_75=round(sum(PloidyMax<0.75)/n(),3),
            Ploidy_1=round(sum(PloidyMin>0.75&PloidyMax<1.9)/n(),3),
            Ploidy_GT_1_5=round(sum(PloidyMin>1.5)/n(),3),
            Ploidy_GT_4=round(sum(PloidyMin>4)/n(),3))


csvData = svData %>% filter(ClusterCount>=3&ResolvedType!='LINE')
nrow(csvData)

ploidyData = csvData %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
  summarise(Ploidy_LT_0_25=round(sum(PloidyMax<0.25)/n(),3),
            Ploidy_LT_0_75=round(sum(PloidyMax<0.75)/n(),3),
            Ploidy_1=round(sum(PloidyMin>0.75&PloidyMax<1.9)/n(),3),
            Ploidy_GT_1_5=round(sum(PloidyMin>1.5)/n(),3),
            Ploidy_GT_4=round(sum(PloidyMin>4)/n(),3))

View(ploidyData)

nrow(ploidyData) # 40.5K
nrow(ploidyData2) # 41.7K

nrow(ploidyData %>% filter(Ploidy_1 >= 0.95))
nrow(ploidyData2 %>% filter(Ploidy_1 >= 0.95))


# number of clusters with > 20% of low and high ploidies
nrow(ploidyData %>% filter(Ploidy_LT_0_75 >= 0.2 & Ploidy_1 >= 0.2)) # 362
nrow(ploidyData2 %>% filter(Ploidy_LT_0_75 >= 0.2 & Ploidy_1 >= 0.2)) # 275

nrow(ploidyData %>% filter(Ploidy_LT_0_75 >= 0.1 & Ploidy_1 >= 0.1)) # 757
nrow(ploidyData2 %>% filter(Ploidy_LT_0_75 >= 0.1 & Ploidy_1 >= 0.1)) # 539


nrow(ploidyData %>% filter(Ploidy_LT_0_75 >= 0.2 & Ploidy_GT_1_5 >= 0.2)) # 
nrow(ploidyData2 %>% filter(Ploidy_LT_0_75 >= 0.2 & Ploidy_GT_1_5 >= 0.2)) # 

nrow(ploidyData %>% filter(Ploidy_LT_0_25 >= 0.1 & Ploidy_GT_1_5 >= 0.2)) # 
nrow(ploidyData2 %>% filter(Ploidy_LT_0_25 >= 0.1 & Ploidy_GT_1_5 >= 0.2)) # 
View(ploidyData2 %>% filter(Ploidy_LT_0_25 >= 0.1 & Ploidy_GT_1_5 >= 0.2)) # 

View(csvData %>% filter(SampleId=='CPCT02011065T'&ClusterId==315))


nrow(ploidyData %>% filter((Ploidy_LT_0_25>0|LowCNC>0))) # 8.5K clusters with subclonal or low ploidy SVs
nrow(ploidyData %>% filter((Ploidy_LT_0_25>0|LowCNC>0)&Ploidy_GT_1_5>0)) # 3K clusters with both subclonal or low ploidy and > 1.5 ploidy SV
View(ploidyData %>% filter((Ploidy_LT_0_25>0|LowCNC>0)&Ploidy_GT_1_5>0)) # 3K clusters with both subclonal or low ploidy and > 1.5 ploidy SV

# View(clusters %>% group_by(ResolvedType,Subclonal) %>% count() %>% spread(Subclonal,n) %>% mutate(SubclonalPerc=round(true/(false+true),2)))

CR_Prox=round(sum(grepl('Prox',ClusterReason))/n(),3),
CR_LongDDI=round(sum(grepl('DelDupInv',ClusterReason))/n(),3),
CR_Foldback=round(sum(grepl('Foldback',ClusterReason))/n(),3),
CR_LOH=round(sum(grepl('LOH',ClusterReason))/n(),3),
CR_ComArm=round(sum(grepl('ComArm',ClusterReason))/n(),3))

   

# problematic clusters
nrow(ploidyData %>% filter(Ploidy_LT_0_5>0&Ploidy_GT_1_5>0))


ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)


ensemblTransExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'))

# canonical trans
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'&CanonicalTranscriptId==TransId))

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'&ExonStart<7571679&ExonEnd>7571754))
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'&ExonStart<=7571678))
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'&ExonStart<=7571678&ExonEnd>=7571679))
View(ensemblTransExonData %>% filter(ExonStart<=7571678&ExonEnd>=7571754))

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000141510'&TransStart<=7571678&TransEnd>=7571754))

View(ensemblTransExonData %>% mutate(ExonLen=ExonEnd-ExonStart) %>% group_by(ExonLen=ifelse(ExonLen<500,round(ExonLen,-1),510)) %>% count())

View(svData %>% filter(SampleId=='CPCT02020378T'&(ChrStart==17|ChrEnd==17)))





linxClustersData = read.csv('~/data/sv/vis/LKCGP-P004801-275546-01-01-21-D1.linx.clusters.tsv',sep='\t')
View(linxClustersData)

linxSvData = read.csv('~/data/sv/vis/LKCGP-P004801-275546-01-01-21-D1.linx.svs.tsv',sep='\t')
View(linxSvData)

sampleSvData = read.csv('~/data/sv/vis/LKCGP-P004801-275546-01-01-21-D1.sv_data.tsv',sep='\t')
View(sampleSvData)
View(sampleSvData %>% group_by(type) %>% count())

svDriverCoc = read.csv('~/data/sv/sv_driver_gene_cooccurrence.csv')
View(svDriverCoc %>% filter(CancerType!='All'))

View(clusters %>% group_by(ResolvedType) %>% count())
View(clusters %>% group_by(ResolvedType,Subclonal) %>% count() %>% spread(Subclonal,n))

View(clusters %>% filter(ResolvedType=='NONE'|ResolvedType=='COMPLEX') %>% group_by(ResolvedType,ClusterSize=2**round(log(ClusterCount,2))) %>% count() %>% spread(ResolvedType,n))

View(clusters %>% filter((ResolvedType=='NONE'|ResolvedType=='COMPLEX')&ClusterCount<=2) %>% group_by(ClusterDesc,Subclonal) %>% count() %>% spread(Subclonal,n))


View(clusters %>% filter(ClusterCount==1&ClusterDesc=='SGL') %>% group_by(ResolvedType,Subclonal) %>% count() %>% spread(Subclonal,n))

View(svData %>% filter(ClusterCount>0) %>% group_by(SampleId,ClusterId,ChainId) %>% summarise(SvCount=n(), ChainCount=first(ChainCount))
     %>% filter(ChainCount>0&ChainCount!=SvCount))

View(clusters %>% filter(ClusterCount > 1000))

View(clusters)
View(clusters %>% filter(ClusterCount>0&(AcFbTI>0|AcFbDSB>0|AcFbPairSame>0|AcFbPairOpp>0|AcFbPairFacing>0|AcComplexFb>0|AcComplexLine>0)))

View(svData %>% filter(ResolvedType=='LOW_VAF'))
View(svData %>% filter(ResolvedType=='LOW_VAF') %>% group_by(Type,VAF=round(CNChgStart/CNStart,2)) %>% count())
View(svData %>% filter(ResolvedType=='LOW_CNC',Ploidy<0.5) %>% group_by(Type,LenBucket=ifelse(Length>0,2**round(log(Length,2)),0)) %>% count() %>% spread(Type,n))
View(svData %>% filter(ResolvedType=='LOW_CNC',Ploidy>0.8) %>% group_by(Type,LenBucket=ifelse(Length>0,2**round(log(Length,2)),0)) %>% count() %>% spread(Type,n))

View(svData %>% filter(ResolvedType=='SglPair_DEL') %>% group_by(Type) %>% count())
View(clusters %>% filter(ResolvedType=='SglPair_DEL'))

synSglPairLengths = clusters %>% filter(ResolvedType=='SGL_PAIR_DEL'|ResolvedType=='SGL_PAIR_DUP') %>% 
  mutate(SyntheticLen=2**round(log(SynDelDupLen,2)))

View(synSglPairLengths)
View(synSglPairLengths %>% group_by(ResolvedType,SyntheticLen) %>% count() %>% spread(ResolvedType,n))

print(plot_length_facetted(synSglPairLengths, '', 'ResolvedType,SyntheticLen', 'SyntheticLen', 'ResolvedType', 'Synthetic DEL/DUP Lengths'))

nrow(svData)

simpleDelsAndDups = svData %>% filter(ClusterCount==1&(Type=='DEL'|Type=='DUP')&ResolvedType=='SIMPLE')
simpleDelsAndDups = simpleDelsAndDups %>% mutate(LengthBucket=ifelse(Length>0,2**round(log(Length,2)),0))

print(plot_length_facetted(simpleDelsAndDups, '', 'Type,LengthBucket', 'LengthBucket', 'Type', 'DEL/DUP Lengths'))
print(plot_length_facetted(simpleDelsAndDups, 'Length>1e3', 'Type,LengthBucket', 'LengthBucket', 'Type', 'DEL/DUP Lengths'))

View(clusters %>% filter(ResolvedType=='SGL_PAIR_DEL'|ResolvedType=='SGL_PAIR_DUP'))

sglPairsData = read.csv('~/logs/sgl_pairs.csv')
sglPairs = sglPairsData %>% group_by(SampleId,ClusterId,ResolvedType,ClusterDesc) %>% 
       summarise(SglCount=sum(Type=='SGL'),
                 NoneCount=sum(Type=='NONE'),
                 SyntheticLength=abs(first(PosStart)-last(PosStart)))

View(sglPairs)
View(sglPairs %>% group_by(ResolvedType,NoneCount,SynLengthBucket=2**round(log(SyntheticLength,2))) %>% count() %>% spread(NoneCount,n))
View(sglPairs %>% group_by(SampleId) %>% summarise(Total=n(),NonePairs=sum(NoneCount==2)))

rm(sglPairsData) 
rm(sglPairs)

View(svData %>% filter(ClusterCount>1) %>% group_by(SampleId,ClusterId,ResolvedType) %>% 
       summarise(ClusterCount=first(ClusterCount),
                 LowCNCSupport=sum(CNChgStart<0.5&CNChgEnd<0.5)) %>% filter(LowCNCSupport>0.51*ClusterCount))



view_cluster_sv('CPCT02230120T',77)
view_chromosome_sv('CPCT02040294T','5','5')
View(view_cluster_arm_breakends('CPCT02230120T', "", '7', 'Q'))


view_cluster_sv('COLO829T',68)
view_chromosome_sv('COLO829T',3,3)


view_cluster_sv('CPCT02030286T',85)
view_chromosome_sv('CPCT02010665T','20','20')

svSpecificData = read.csv('~/logs/CPCT02020258T_SVA.csv')
svSpecificData = read.csv('~/logs/COLO829T_SVA.csv')
View(svSpecificData)
View(svSpecificData %>% filter(ClusterCount>1) %>% group_by(ClusterId,ClusterCount) %>% count())



svSpecificData = read.csv('~/logs/CPCT02040228T_SVA.csv')
View(view_cluster_arm_breakends(svSpecificData, 'CPCT02050243T', "462", '12', 'Q'))

View(svSpecificData)
View(svSpecificData %>% filter(GeneStart!=''|GeneEnd!='') %>% group_by(GeneStart,GeneEnd) %>% count())

View(svSpecificData %>% filter(ChrStart==9|ChrEnd==9))
View(svSpecificData %>% filter(ClusterId==9))
View(svSpecificData %>% filter(ClusterId==9) %>% group_by(PloidyBucket=ifelse(Ploidy<1,round(Ploidy/0.1)*0.1,round(Ploidy))) %>% count())


View(view_cluster_arm_breakends(svSpecificData, 'CPCT02050243T', "86", '9', 'Q'))

rm(svData)
rm(clusters)



# SHATTERING SIMS
simulations = read.csv('~/logs/SIM_RESULTS.csv')
nrow(simulations)
View(simulations)
View(simulations %>% group_by(SegsLinked,AdjacentSegs) %>% count() %>% spread(AdjacentSegs,n))
View(simulations %>% group_by(SegsLinked,ExactMatches) %>% count() %>% spread(ExactMatches,n))
simulations$SegsLinkedPerc = round(simulations$SegsLinked/simulations$SegCoun,5)

nrow(simulations %>% filter(SegsLinked==SegCount))
nrow(simulations %>% filter(SegsLinked==SegCount-1))
nrow(simulations %>% filter(SegsLinked==SegCount-2))
nrow(simulations %>% filter(SegsLinked==SegCount-3))

write.csv(simulations %>% group_by(SegsLinkedPerc) %>% count(), '~/logs/shattering_sim_segs_linked.csv', row.names = F, quote = F)

print(ggplot(data = simulations %>% group_by(SegsLinkedPerc) %>% count(), aes(x=SegsLinkedPerc, y=n))
      + geom_point()
      + labs(title = "Segments linked vs Total Segments"))

print(ggplot(data = simulations %>% group_by(SegsLinked,AdjacentSegs) %>% count(), aes(x=SegsLinked, y=AdjacentSegs))
      + geom_point()
      + labs(title = "Segments linked vs Adjacent Segments"))

print(ggplot(data = simulations %>% group_by(AdjacentSegs) %>% count(), aes(x=AdjacentSegs, y=n))
      + geom_point()
      + labs(title = "Adjacent Segments"))

print(ggplot(data = simulations %>% group_by(SegsLinked,ExactMatches) %>% count(), aes(x=SegsLinked, y=ExactMatches))
      + geom_point()
      + labs(title = "Segments linked vs Exact Repairs"))


# DOUBLE MINUTES FROM DUPS
simpleDups = svData %>% filter(ClusterCount==1&Type=='DUP')
simpleDups$PloidyBucket = round(simpleDups$Ploidy)
View(simpleDups %>% filter(PloidyBucket>1) %>% group_by(SampleId,PloidyBucket) %>% count())


dupDMData = read.csv('~/logs/dup_dm_data.csv')
View(dupDMData)
View(dupDMData %>% filter(PloidyMin>1.75,PloidyMin>MajorAPStart+1&PloidyMin>MajorAPEnd+1,
                          PloidyMin>ifelse(BafCountStart>1,MajorAPStart,MajorAPStart+MinorAPStart)*2.5&
                            PloidyMin>ifelse(BafCountEnd>1,MajorAPEnd,MajorAPEnd+MinorAPEnd)*2.5,
                          SVWithin<100,SVOverlapping<100) %>% mutate(length=PosEnd-PosStart)
     %>% select(SampleId,SvId,PloidyMin,SamplePurity,SamplePloidy,length,MajorAPEnd,MajorAPStart,BafCountStart,BafCountEnd,DWCStart,DWCEnd,SVOverlapping,SVWithin,everything()))

view_chromosome_sv('CPCT02070290T',7,7)


######
## DOUBLE MINUTES


potentialDMs = read.csv('~/logs/potential_dm_data.csv')
# potentialDMs = read.csv('~/data/sv/potential_dm_data.csv')
View(potentialDMs)

View(potentialDMs %>% filter(IsComplete=='false') %>% group_by(SvTypes) %>% count())
View(potentialDMs %>% group_by(Chromosome,IsComplete) %>% count())
# view_sv_ids(c(10321276,10321278,35093,35133,35088,10320965,35089,10321081))
View(potentialDMs %>% filter(IsComplete=='true') %>% group_by(HasAmpGene=AmplifiedGenes!='') %>% count())

view_sv_id(11226271)
view_chromosome_sv('CPCT02040121T',8,8)
view_chromosome_sv('CPCT02020267TII',14,14)
view_chromosome_sv('CPCT02240013T',9,9)


# PLOIDY RECALCS
ploidyData = read.csv('~/logs/CN_PLOIDY_CALC_DATA.csv')
View(ploidyData)
View(ploidyData %>% filter(MinPloidy==MaxPloidy))

View(ploidyData %>% filter(ChrStart==11&Type=='DEL')
       %>% select(SvId,Ploidy,CNChgStart,CNChgEnd,PloidyLow,PloidyHigh,PloidyUncertainty,EstPloidy,EstUncertainty,MinPloidy,MaxPloidy))

nonePloidyData = svData %>% filter(Type=='NONE'&PloidyMin==PloidyMax)
View(nonePloidyData)
View(svData %>% filter(SampleId=='CPCT02010122T'))
View(svData %>% filter(Type=='NONE'))
View(svData %>% filter(PloidyMin>=PloidyMax&Type!='INS'))

write.csv(svData %>% filter(SampleId=='CPCT02040204T'&ClusterId==1) %>% select(Id,AdjCNChgStart,AdjCNChgEnd,Ploidy), '~/logs/tmp/')

lowPloidyData = svData %>% filter(ResolvedType!='LowQual'&(Ploidy == 0|CNChgStart==0|(CNChgEnd==0&PosEnd>0)))
nrow(lowPloidyData)
View(lowPloidyData)

lowQualSVs = svData %>% filter(ResolvedType=='LowQual')
nrow(lowQualSVs) # 78K

# reasons for LowQual and types:

nrow(svData %>% filter(ResolvedType=='LowQual'))

View(svData %>% filter(Ploidy==PloidyMin&Ploidy==PloidyMax))
View(svData %>% filter(Type!='NONE'&Type!='SGL'&(AFStart==0|AFEnd==0)))
View(svData %>% filter(AFStart==0) %>% group_by(Type) %>% count())
View(svData %>% filter(Id==9754792))

View(svData %>% filter(SampleId=='CPCT02140095T'&ClusterId==87))

specSampleData = read.csv('~/logs/CPCT02010003T_SVA.csv')
View(specSampleData)
View(specSampleData %>% filter(Ploidy==PloidyMin&Ploidy==PloidyMax))

view_cluster_sv('CPCT02140095T', 49)


zeroCalcPloidyData = svData %>% filter(PloidyMin<=0&PloidyMax>5)
View(zeroCalcPloidyData)

calc_cn_change_buckets<-function(svData, sampleId, clusterId, forceIntegers = T)
{
  sampleData = svData %>% filter(SampleId==sampleId&ClusterId==clusterId) %>% select(Id,AdjCNChgStart,AdjCNChgEnd,Ploidy)
  sampleData$ConsCNChg = ifelse(sampleData$AdjCNChgStart<0.5,sampleData$AdjCNChgEnd,
                                ifelse(sampleData$AdjCNChgEnd<0.5,sampleData$AdjCNChgStart,round((sampleData$AdjCNChgStart+sampleData$AdjCNChgEnd)*0.5,2)))
  
  minRounded = round(min(sampleData$ConsCNChg),1)
  
  if(forceIntegers)
    minRounded = round(minRounded)

  sampleData$CNChgBucket = round(sampleData$ConsCNChg/minRounded)*minRounded
  sampleData$CNChgResidual = abs(sampleData$CNChgBucket-sampleData$ConsCNChg)
  
  summary = (sampleData %>% group_by(CNChgBucket) %>% summarise(Count=n(), 
                                                          Residuals=sum(CNChgResidual),
                                                          AvgResiduals=round(sum(CNChgResidual)/n(),2),
                                                          AvgResidualsPerc=round(sum(CNChgResidual)/n()/first(CNChgBucket),2)))
  
  return (summary)
}

View(calc_cn_change_buckets(svData, 'COLO829T', 68))
View(calc_cn_change_buckets(svData, 'CPCT02040204T', 1))

rm(svData)

sampleData = svData %>% filter(SampleId=='CPCT02040204T'&ClusterId==1) %>% select(Id,AdjCNChgStart,AdjCNChgEnd,Ploidy)
sampleData$ConsCNChg = ifelse(sampleData$AdjCNChgStart<0.5,sampleData$AdjCNChgEnd,
                              ifelse(sampleData$AdjCNChgEnd<0.5,sampleData$AdjCNChgStart,round((sampleData$AdjCNChgStart+sampleData$AdjCNChgEnd)*0.5,2)))

minRounded = round(min(sampleData$ConsCNChg),1)
minRounded = round(minRounded)

sampleData$CNChgBucket = round(sampleData$ConsCNChg/minRounded)*minRounded
sampleData$CNChgResidual = abs(sampleData$CNChgBucket-sampleData$ConsCNChg)
View(sampleData)


chromosomes = c(1,5,6,7,21)
specificSample = 'CPCT02010003T'
plot_cross_chr_sample_clusters(svData, chromosomes, specificSample)

lowCNChgData = svData %>% filter(Type!='SGL'&Type!='NONE') %>% filter(Ploidy>0.8&(AdjCNChgStart<0.2|AdjCNChgEnd<0.2)&(AdjCNChgStart>0.8|AdjCNChgEnd>0.8))
nrow(lowCNChgData)
View(lowCNChgData)
lowCNChgData$PloidyCNChgRatio = round(ifelse(lowCNChgData$AdjCNChgStart<0.2,lowCNChgData$Ploidy/lowCNChgData$AdjCNChgEnd,lowCNChgData$Ploidy/lowCNChgData$AdjCNChgStart),1)

View(lowCNChgData %>% group_by(PloidyCNChgRatio) %>% count())

lowCNChgData2 = svData %>% filter(Type!='SGL'&Type!='NONE') %>% filter(Ploidy>0.5&(AdjCNChgStart<0.2&AdjCNChgEnd>0.2)|(AdjCNChgStart>0.2&AdjCNChgEnd<0.2))
nrow(lowCNChgData2)
lowCNChgData2$PloidyCNChgRatio = round(ifelse(lowCNChgData2$AdjCNChgStart<0.2,lowCNChgData2$Ploidy/lowCNChgData2$AdjCNChgEnd,lowCNChgData2$Ploidy/lowCNChgData2$AdjCNChgStart),1)
View(lowCNChgData2 %>% group_by(PloidyCNChgRatio) %>% count())




# Line analysis following Suspect changes

# svData = sv_load_and_prepare('~/logs/CLUSTER.csv')
svData = sv_load_and_prepare('~/data/sv/SVA_SVS.csv')

View(svData %>% filter(SampleId=='CPCT02030459TII'))

view_chromosome_sv('CPCT02030459TII', 17, 17)
view_cluster_sv('CPCT02030459TII', 167)

sgls = svData %>% filter(Type=='SGL'|Type=='NONE')
sgls$LowVaf = sgls$AdjAFStart < 0.15
sgls$PolyC = grepl('CCCCCCCCCC',sgls$InsertSeq)
sgls$PolyG = grepl('GGGGGGGGGG',sgls$InsertSeq)
View(sgls %>% group_by(LowVaf,PolyC,PolyG) %>% count())


startShortTI = svData %>% filter(LnkLenStart>0&LnkLenStart<400)
nrow(startShortTI)
endShortTI = svData %>% filter(LnkLenEnd>0&LnkLenEnd<400)
nrow(endShortTI)
View(endShortTI)

view_cluster_sv('COLO829T',66)

View(svData %>% filter(IsLINE) %>% group_by(LEStart,LEEnd) %>% count())

lineStart = svData %>% filter(grepl('Known',LEStart)|grepl('Suspect',LEStart)|grepl('Ident',LEStart))
lineStart$LineStatus = lineStart$LEStart
lineEnd = svData %>% filter(grepl('Known',LEEnd)|grepl('Suspect',LEEnd)|grepl('Ident',LEEnd))
lineEnd$LineStatus = lineEnd$LEEnd
lineAll = rbind(lineStart,lineEnd)

View(lineAll %>% group_by(LineStatus) %>% count())

lineClusters = svData %>% filter(IsLINE) %>% group_by(SampleId,ClusterId) %>% count()
lineClusters$SampleClusterId = paste(lineClusters$SampleId, lineClusters$ClusterId, sep='_')
View(lineClusters)

lineData = svData %>% filter(paste(svData$SampleId, svData$ClusterId, sep='_') %in% lineClusters$SampleClusterId)
lineData$Known = (grepl('Known',lineData$LEStart)|grepl('Known',lineData$LEEnd))
lineData$Suspect = (grepl('Suspect',lineData$LEStart)|grepl('Suspect',lineData$LEEnd))
lineData$IndentOnly = (lineData$LEStart=='Ident'|lineData$LEEnd=='Ident')
lineData$PolyAT = (grepl('AAAAAAAA',lineData$InsertSeq)|grepl('TTTTTTTT',lineData$InsertSeq))
lineClusterData$SVtoArmRatio = lineClusterData$ClusterCount/lineClusterData$ArmCount
  
lineClusterData = (lineData %>% group_by(SampleId,ClusterId) 
                   %>% summarise(ClusterCount=first(ClusterCount),
                                 ResolvedType=first(ResolvedType),
                                 ArmCount=first(ArmCount),
                                 BndCount=sum(Type=='BND'),
                                 SglCount=sum(Type=='SGL'),
                                 BndKnownCount=sum(Type=='BND'&Known),
                                 BndSuspectCount=sum(Type=='BND'&Suspect),
                                 KnownCount=sum(Known),
                                 IdentOnlyCount=sum(IndentOnly),
                                 SuspectCount=sum(Suspect),
                                 PolyATCount=sum(PolyAT),
                                 BndPolyAT=sum(Type=='BND'&PolyAT),
                                 OtherPolyAT=sum(Type!='BND'&PolyAT)))

lineClusterData$AllKnown = (lineClusterData$KnownCount==lineClusterData$ClusterCount)

View(lineClusterData)

nrow(lineClusterData %>% filter(ClusterCount==1))
View(lineClusterData %>% filter(ClusterCount==1) %>% group_by(PolyATCount,BndCount) %>% count())
View(lineClusterData %>% filter(ClusterCount==1) %>% group_by(PolyATCount,BndSgl=(BndCount>0|SglCount>0)) %>% count())

View(lineClusterData %>% filter(PolyATCount==0&AllKnown))


print(ggplot(data = lineClusterData, aes(x=ArmCount, y=SuspectCount))
      + geom_point()
      + labs(title = "ArmCount vs Suspect LINE Svs"))

print(ggplot(data = lineClusterData %>% filter(SuspectCount>0), aes(x=ArmCount))
      + geom_point(aes(y=SuspectCount, colour='SuspectCount'))
      + geom_point(aes(y=KnownCount, colour='KnownCount'))
      + geom_point(aes(y=PolyATCount, colour='PolyATCount'))
      + labs(title = "ArmCount vs Suspect/Known LINE Svs"))

print(ggplot(data = lineClusterData %>% filter(SuspectCount>0), aes(x=ArmCount,y=ClusterCount))
      + geom_point()
      + labs(title = "ArmCount vs Cluster Svs"))



view_chromosome_sv('CPCT02020258T', 3, 3)
view_cluster_sv('CPCT02020258T', 52)
view_cluster_sv('CPCT02020258T', 423)

View(sampleSvData)


consecutiveBEs = read.csv('~/logs/consecutive_be.csv')
View(consecutiveBEs)
consecutiveBEs$HasBEBetween = (consecutiveBEs$TraversedSVs>0)
consecutiveBEs$NextClusteredSVOrientMatch = (consecutiveBEs$Orientation==consecutiveBEs$OrientNextCL)
consecutiveBEs$NextNonClusteredSVOrientMatch = (consecutiveBEs$Orientation==consecutiveBEs$OrientNext)
consecutiveBEs$NextClusteredSVCNMatch = ifelse(abs(consecutiveBEs$CNChgFront-consecutiveBEs$CNChgNextCL) < 0.5,T,F)
consecutiveBEs$NextNonClusteredSVCNMatch = ifelse(abs(consecutiveBEs$CNChgFront-consecutiveBEs$CNChgNext) < 0.5,T,F)

View(consecutiveBEs %>% group_by(HasBEBetween,NextClusteredSVOrientMatch,NextNonClusteredSVOrientMatch,NextClusteredSVCNMatch,NextNonClusteredSVCNMatch) %>% count())



view_chromosome_sv('COLO829T', 3, 3)

view_cluster_sv('CPCT02040204T',95)


bachTmp = read.csv("~/dev/bachelor/CPCT02190003T.output.csv")
View(bachTmp)
View(bachTmp %>% group_by(X22) %>% count())


sampleData = sv_load_and_prepare('~/logs/COLO829T.csv')

# check deletion bridge presence and length
View(foldbacks)

View(foldbacks %>% filter(FoldbackType=='INV'&OrientStart==1&DBLenStart>=-30&DBLenStart<=30))
View(foldbacks %>% filter(FoldbackType=='INV'&OrientStart==-1&DBLenEnd>=-30&DBLenEnd<=30))

View(foldbacks %>% filter(FoldbackType=='INV'&FoldbackLength<100))

View(tiDirectData %>% filter(TILength<8e3) %>% group_by(SampleId) %>% count() %>% arrange(-n))

-view_chromosome_sv('CPCT02010797T', 9, 9)

view_cluster_sv('CPCT02040204T',95)
view_chromosome_sv('CPCT02040204T', 11,11)

view_sv_id(4426884)

View(svData %>% filter(grepl('asm',AsmbStart)&DBLenStart>=-30&DBLenStart<0))

View(svData %>% filter(grepl('asm',AsmbStart)&DBLenStart>=-30&DBLenStart<0))

View(svData %>% filter(SampleId=='CPCT02140070T'&SubClusterId==272))

view_chromosome_sv('CPCT02140070T', 5, 5)
chartSample = 'CPCT02140070T'
chartChr = 5
print(ggplot(data = svData %>% filter(ResolvedType!='LowQual',SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) %>% mutate(
  modStart=ifelse(ChrStart==chartChr,PosStart,as.numeric(ChrStart)*-1e6),modEnd=ifelse(ChrEnd==chartChr,PosEnd,as.numeric(ChrEnd)*-1e6)),aes(modStart,modEnd)) + 
    geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())

View(svData %>% filter(SampleId=='CPCT02140070T'&ResolvedType!='LowQual'&ClusterDesc!='NONE'&ClusterDesc!='SGL'&(ChrStart==5|ChrEnd==5)) 
     %>% group_by(ClusterId) %>% summarise(ClusterCount=first(ClusterCount), 
                                           ClusterDesc=first(ClusterDesc),
                                           ResolvedType=first(ResolvedType),
                                           SvLength=first(Length),
                                           NearestType=first(NearestType)))

view_cluster_sv('CPCT02030482T',102)
View(svData %>% filter(SampleId=='CPCT02140070T'&ClusterId==197))

plot_cross_chr_sample_clusters(svData, c(5), 'CPCT02140070T', 1, T)

view_chromosome_sv('CPCT02060269T', 14, 14)

DRUP01210003T
view_sv_ids(c(4487963,4487977))

view_cluster_sv('DRUP01110005T', 132)
view_chromosome_sv('CPCT02010842T',5,5)

View(svData %>% filter(ResolvedType=='DUP_Int_TI') %>% group_by(ClusterDesc) %>% count())
View(svData %>% filter(ResolvedType=='DUP_Int_TI'&ClusterDesc=='DEL=1_DUP=1'))

View(svData %>% filter(ClusterDesc=='DUP=2'&ResolvedType=='DUP_Ext_TI'))

View(tiDirectData %>% filter(SampleId=='DRUP01110005T'&ClusterId==132))


clusters = read.csv('~/data/sv/SVA_CLUSTERS.csv')
View(clusters %>% filter(ResolvedType=='DUP_Int_TI'|ResolvedType=='DEL_Int_TI'|ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI'))


View(clusters %>% filter(ResolvedType=='SimpleChain'&ClusterCount>=3&ClusterCount<=4&FullyChained=='true'
                         &Consistency==0&OriginArms==1&HasReplicated=='false'))

syntheticDelDupChains = (clusters %>% filter(ResolvedType=='SimpleChain'&ClusterCount>=3&ClusterCount<=10&FullyChained=='true'
                                                 &Consistency==0&OriginArms==1&HasReplicated=='false'))
nrow(syntheticDelDupChains)

View(syntheticDelDupChains %>% filter(AssemblyLinks==ClusterCount-1))

syntheticDelDupChains$SampleClusterId = paste(syntheticDelDupChains$SampleId,syntheticDelDupChains$ClusterId,sep='_')

sddSvData = svData %>% filter(paste(svData$SampleId,svData$ClusterId,sep='_') %in% syntheticDelDupChains$SampleClusterId)
View(sddSvData)


clusters = read.csv('~/logs/SVA_CLUSTERS.csv')
View(clusters %>% filter(ResolvedType=='DUP_Int_TI'|ResolvedType=='DEL_Int_TI'|ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI'))

View(clusters %>% filter(ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI') %>% filter(ClusterCount>2))

View(clusters %>% filter(ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI') %>% group_by(ClusterCount) 
     %>% summarise(Count=n(), SvCount=sum(ClusterCount)))

# compare to proposed list
View(clusters %>% filter(ResolvedType=='SimpleChain'&ClusterCount>=3&ClusterCount<=10&FullyChained=='true'&ChainCount==1
                         &Consistency==0&OriginArms==1&HasReplicated=='false'))


View(clusters %>% filter(ResolvedType=='SimpleChain'&ClusterCount>=3&ClusterCount<=10&FullyChained=='true'&ChainCount==1
                         &Consistency==0&OriginArms==1&HasReplicated=='false'))


synthDelDups = clusters %>% filter(ResolvedType=='DUP_Ext_TI'|ResolvedType=='DEL_Ext_TI')
View(synthDelDups)

synthDelDups$ClusterSize = ifelse(synthDelDups$ClusterCount==2,'Pair','Longer')
synthDelDups$TiLenBucket = 2**round(log(synthDelDups$SynDelDupAvgTILen,2))
synthDelDups$LenBucket = ifelse(synthDelDups$SynDelDupLen<=0,0,2**round(log(synthDelDups$SynDelDupLen,2)))




print(ggplot(data = synthDelDups %>% group_by(TiLenBucket,ClusterSize) %>% count(), aes(x=TiLenBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterSize)
      + labs(title = "TI Lengths for Synthetic DELs and DUPs"))

print(ggplot(data = synthDelDups %>% filter(SynDelDupAvgTILen<=2e4) %>% group_by(LenBucket,ClusterSize,ResolvedType) %>% count() %>% spread(ResolvedType,n), aes(x=LenBucket, y=n))
      + geom_line(aes(y=DEL_Ext_TI, colour='DEL_Ext_TI'))
      + geom_line(aes(y=DUP_Ext_TI, colour='DUP_Ext_TI'))
      + scale_x_log10()
      + facet_wrap(~ClusterSize)
      + labs(title = "Synthetic Lengths for Synthetic DELs and DUPs"))

View(cleanClusters %>% filter(ResolvedType=='DEL'&Synthetic=='true'&ClusterCount>5))

View(dnaRnaCombinedData %>% filter(HmfId=='HMF001945A'))
View(dnaRnaCombinedData %>% filter(SampleId=='CPCT02050063T'))

View(dnaRnaCombinedData %>% filter(HmfId=='HMF001592B'))
View(dnaRnaCombinedData %>% filter(HmfId=='HMF001945A'))


View(dnaRnaCombinedData %>% filter(KnownCategory!='Unknown') %>% 
       group_by(HasRnaData,HasDnaData,KnownCategory,PhaseMatched.x,PhaseMatched.y,RnaPhaseMatched) %>% count)

# filter out unphased DNA fusions if also RNA unphased for all subsequent analysis
rnaMatchData = rnaMatchData %>% filter(PhaseMatched=='true'|RnaPhaseMatched)

View(rnaCombinedData %>% filter(SampleId=='CPCT02050063T'))
View(rnaMatchData %>% filter(SampleId=='CPCT02050063T')) 
View(svaRnaFusions %>% filter(SampleId=='DRUP01180001T')) 

View(svaRnaFusions %>% filter(SampleId=='CPCT02140030T'&GeneNameDown=='NTRK3') %>% 
       select(GeneNameUp,GeneNameDown,BreakendDistDown,DistancePrevDown,PosDown,StrandDown,TransStartDown,TransEndDown, everything())) 
colnames(svaRnaFusions)
View(svaRnaFusions)

View(dnaRnaCombinedOutputData %>% filter(is.na(BreakendExonUp)&MatchType!='RNA Only') %>% group_by(KnownCategory) %>% count)
View(dnaRnaCombinedOutputData %>% filter(SampleId=='DRUP01180001T'))
View(dnaRnaCombinedOutputData %>% filter(RnaP))

#svData = read.csv('~/logs/CLUSTER_GRIDSS.csv')
#nrow(svData)
#View(svData)

view_sv_id(4547418)
view_sv_ids(c(4546677,4546981))

# simple annotations
svData = sv_set_common_fields(svData)

# validate derived fields
View(head(svData %>% filter(InferTICount==0.5),100))
View(head(svData %>% filter(AsmbTICount==0.5),100))
View(head(svData %>% filter(InferTICount==1),100))
View(head(svData %>% filter(AsmbTICount==1),100))
View(head(svData %>% filter(IsChained > 0),100))
View(head(svData %>% filter(Length==-1),100))
View(head(svData %>% filter(RepeatedChainLinkCount > 0),100))

svResolvedSummary = (svData %>% group_by(ResolvedType,ClusterSize)
                     %>% summarise(Count=n()) %>% arrange(ResolvedType,ClusterSize))

totalSVCount = nrow(svData)
View(svData %>% group_by(ResolvedType,ClusterSize) 
     %>% summarise(Clusters=n_distinct(paste(SampleId,ClusterId,sep='_')), TotalSVs=n(), AsPerc=round(n()/totalSVCount,2))
     %>% arrange(-AsPerc))

View(svData %>% filter(ClusterCount==1&ResolvedType=='None'))
nrow(svData %>% filter(IsLINE))

svData$ChrArm = paste(svData$SampleId,svData$ChrStart,sep='_')
View(svData %>% filter(ClusterCount>20) %>% group_by(SampleId,ClusterId,ChrArm) 
     %>% summarise(SvCount=n(), LineCount=sum(IsLINE)) %>% group_by(SampleId,ClusterId) %>% summarise(SvCount=sum(SvCount),ArmCount=n(),LineCount=sum(LineCount)))


view_cluster_sv('CPCT02020670TII',3)


view_cluster_sv('COLO829T',66)

view_sv_ids(c(416884,416881)) 


view_cluster_sv('CPCT02010003T',472)

view_cluster_sv('CPCT02010003T',67)
view_chromosome_sv('CPCT02010268T','3', '')


view_sv_id(405134)


View(svData %>% filter(as.character(ChrStart)==as.character(ChrEnd),ArmStart!=ArmEnd,ResolvedType!='LowQual',ResolvedType!='Line') 
     %>% group_by(SampleId,ChrStart,ClusterId) %>% count() 
     %>% group_by(SampleId,ChrStart) %>% count())

View(svData %>% filter(IsLINE&AsmbMatchStart=='MATCH'))


View(svData %>% filter(SampleId=='CPCT02010003T'&(ChrStart==1|ChrStart==5|ChrStart==6)) %>% arrange(ChrStart,PosStart))
View(svData %>% filter(SampleId=='CPCT02010003T'&(ChrEnd==1|ChrEnd==5|ChrEnd==6)&round(AdjCNChgStart)==2) %>% arrange(ChrEnd,PosEnd))
nrow(svData %>% filter(SampleId=='CPCT02010003T'&IsLINE&(ChrStart==1|ChrStart==5|ChrStart==6|ChrEnd==1|ChrEnd==5|ChrEnd==6)) %>% arrange(ChrStart,PosStart))

nrow(svData %>% filter(SampleId=='CPCT02010003T'&IsLINE&(ChrStart==1|ChrStart==5|ChrStart==6|ChrEnd==1|ChrEnd==5|ChrEnd==6)) %>% arrange(ChrStart,PosStart))

View(svData %>% filter(SampleId=='CPCT02010003T'&(ChrStart==1|ChrStart==5|ChrStart==6|ChrEnd==1|ChrEnd==5|ChrEnd==6)) 
     %>% group_by(ChrStart,ChrEnd) %>% summarise(SvCount=n(),Clusters=n_distinct(ClusterId)))


enclosedTIs = read.csv("~/logs/enclosed_tis.csv")
View(enclosedTIs)

View(enclosedTIs %>% filter(StartGap==0|EndGap==0))

View(enclosedTIs %>% group_by(TI1Match,TI2Match) %>% count())
View(enclosedTIs %>% group_by(TI1Match,TI2Match,SameCluster=(ClusterId!='multiple')) %>% count())

startGaps = enclosedTIs
startGaps$GapLength = startGaps$StartGap
endGaps = enclosedTIs
endGaps$GapLength = endGaps$EndGap

allGapLengths = rbind(startGaps,endGaps)
allGapLengths$GapLenBucket = ifelse(allGapLengths$GapLength>0,2**round(log(allGapLengths$GapLength)),0)
allGapLengths$SameCluster = (allGapLengths$ClusterId!='multiple')
allGapLengths$TIMatch = (allGapLengths$TI1Match=='true'&allGapLengths$TI2Match=='true')

print(ggplot(data = allGapLengths %>% group_by(GapLenBucket,TIMatch) %>% count(), aes(x=GapLenBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~TIMatch))

print(ggplot(data = allGapLengths %>% group_by(GapLenBucket,SameCluster) %>% count(), aes(x=GapLenBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~SameCluster))

print(ggplot(data = enclosedTIs, aes(x=StartGap, y=EndGap))
      + geom_point()
      + scale_x_log10()
      + scale_y_log10())

print(ggplot(data = enclosedTIs %>% filter(TI1Match=='true'&TI2Match=='true'&ClusterId!='multiple'), aes(x=StartGap, y=EndGap))
      + geom_point()
      + scale_x_log10()
      + scale_y_log10())

print(ggplot(data = enclosedTIs enclosedTIsaes(modStart,modEnd)) 
      + geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())




print(ggplot(data = svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) 
             %>% mutate(modStart=ifelse(ChrStart==chartChr,PosStart,as.numeric(ChrStart)*-1e6),
                        modEnd=ifelse(ChrEnd==chartChr,PosEnd,as.numeric(ChrEnd)*-1e6)),aes(modStart,modEnd)) 
      + geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())



View(svData %>% filter(ResolvedType=='Line'&ClusterCount==1))
View(svData %>% filter(ResolvedType=='Line'&ClusterCount>1&ClusterCount<4&(LEStart=='Suspect'|LEEnd=='Suspect')))
View(svData %>% filter(ResolvedType=='Line'&ClusterCount>1&(LEStart=='Suspect'|LEEnd=='Suspect')))


sampleData = read.csv('~/logs/CPCT02070194T.csv')

View(sampleData %>% filter(SampleId=='CPCT02070194T'&ChrStart==11&ChrEnd==17))
View(sampleData %>% filter(SampleId=='CPCT02070194T'&(ChrStart==11|ChrEnd==17)))
View(svData %>% filter(SampleId=='COLO829T'&ClusterCount>5))

sampleData = read.csv('~/logs/CPCT02070176T.csv')
View(sampleData %>% filter(SampleId=='CPCT02070176T'&(Id==522800|Id==522801)))

View(svData %>% filter(SampleId=='CPCT02070194T'&SubClusterId==330))
View(svData %>% filter(SampleId=='CPCT02070194T'&SubClusterId==330) %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,ClusterId,SubClusterId))
View(svData %>% filter(SampleId=='CPCT02070194T'&SubClusterId==331) %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,ClusterId,SubClusterId))


# are line elements largely isolated?
clusterData = svData %>% group_by(SampleId,ClusterId) %>% summarise(SvCount=n(),LineCount=sum(LEStart!='None'|LEEnd!='None'))
View(clusterData %>% filter(LineCount>0))


# INV foldbacks incorrect clustering?
foldbackClusters = (svData %>% group_by(SampleId,ClusterId) 
                    %>% summarise(SvCount=n(),
                                  FBCountCount=sum(FoldbackLenStart>0|FoldbackLenStart>0),
                                  CLusterReasonFB=sum(grepl('Foldback',ClusterReason))))

View(foldbackClusters)

View(clusterData %>% filter(LineCount>0))


View(svData %>% filter(Type=='SGL'&ClusterDesc=='SGL=2'&ResolvedType=='None')
     %>% group_by(SampleId,ClusterId) %>% arrange(PosStart) 
     %>% summarise(SynLength=last(PosStart)-first(PosStart),
                   IsDUP=sum(ifelse(last(OrientStart)==1&first(OrientStart)==-1,1,0)),
                   IsDEL=sum(ifelse(last(OrientStart)==-1&first(OrientStart)==1,1,0)),
                   IsNeither=sum(ifelse(last(OrientStart)==first(OrientStart),1,0))))

View(svData)



View(svResolvedSummary)

sampleData = read.csv('~/logs/CPCT02020258T.csv')
View(sampleData %>% filter(ChrStart==13&ArmStart=='Q'))
View(sampleData %>% filter(Id==443842))
View(svData %>% filter(Id==489026))



# SAMPLE CLUSTER VIEW
chartSample = "CPCT02330035T"
chartChr = "17"

View(svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr))

print(ggplot(data = svData %>% filter(SampleId==chartSample,ChrStart==chartChr|ChrEnd==chartChr) 
             %>% mutate(modStart=ifelse(ChrStart==chartChr,PosStart,as.numeric(ChrStart)*-1e6),
                        modEnd=ifelse(ChrEnd==chartChr,PosEnd,as.numeric(ChrEnd)*-1e6)),aes(modStart,modEnd)) 
      + geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())


View(svData %>% filter(SampleId=='CPCT02020670T'&(AdjCNChgStart<0.5&AdjCNChgEnd<0.5)&ResolvedType!='LowQual'))
View(svData %>% filter(SampleId=='CPCT02020670TII'&(AdjCNChgStart<0.5&AdjCNChgEnd<0.5)&ResolvedType!='LowQual'))
View(svData %>% filter(SampleId=='CPCT02020670TII'&ClusterId==3))


View(svData %>% filter(SampleId=='CPCT02010641TII'&(ChrStart==7|ChrEnd==7)&Type=='BND'&PosStart>=2762500&PosEnd<=2763500))

View(svData %>% filter(SampleId=='CPCT02010641TII'&((ChrStart==7&PosStart>=2759000&PosStart<=2769000)|(ChrEnd==7&PosEnd>=2759000&PosEnd<=2769000)))
     %>% select (SampleId,Id,Type,ChrEnd,PosEnd,AdjCNChgEnd,ClusterId,SubClusterId))

View(svData %>% filter(SampleId=='CPCT02070064TII'&ClusterCount==21))



View(svData %>% filter(ResolvedType!='LowQual'&grepl("eqv",AsmbStart)&Type=='SGL'))


# events per arm
svArmData = svData
svArmData$Chr = svArmData$ChrStart
svArmData$Arm = svArmData$ArmStart
svArmBndData = svData %>% filter(Type=='BND')
svArmBndData$Chr = svArmBndData$ChrEnd
svArmBndData$Arm = svArmBndData$ArmEnd
combinedArmData = rbind(svArmData, svArmBndData)

svArmSummary = (combinedArmData %>% group_by(SampleId,Chr,Arm,ClusterId)
                %>% summarise(SvCount=n(),
                              ClusterCount=first(ClusterCount),
                ))
View(svArmSummary)

svArmSummary2 = svArmSummary %>% group_by(SampleId,Chr,Arm) %>% summarise(SvCount=sum(SvCount), ClusterCount=n())
View(svArmSummary2)

# arm stats
svArmStats = (combinedArmData %>% group_by(SampleId,Chr,Arm)
              %>% summarise(SvCount=n(),
                            NumOfClusters=n_distinct(ClusterId),
                            MaxClusterSize=max(ClusterCount),
                            MaxChainSize=max(ChainCount),
                            ConsistencyTotal=sum(Consistency),
                            LinkedCount=sum(IsDB|IsTI),
                            AssemblyCount=sum(AssemblyLinked),
                            ChainedCount=sum(IsChained),
                            DBCount=sum(IsDB),
                            TICount=sum(IsTI),
                            DelCount=sum(Type=='DEL'),
                            DupCount=sum(Type=='DUP'),
                            InsCount=sum(Type=='INS'),
                            IsolatedSvCount=sum(NearestType=='NHBR'&(Type=='DEL'|Type=='DUP'|Type=='INS')),
                            InvCount=sum(Type=='INV'),
                            BndCount=sum(Type=='BND'),
                            SglCount=sum(Type=='SGL'),
                            LineCount=sum(IsLINE),
                            FragileSiteCount=sum(IsFS),
                            StressedCount=sum(IsStressed))
              %>% arrange(SampleId,Chr,Arm))

View(svArmStats)

View(rbind(svData %>% unite(ChrArm,ChrStart,ArmStart) %>% filter(ClusterCount==1) %>%
             group_by(SampleId,ChrArm,Type,Id=ClusterId) %>% count(), svData %>%  unite(ChrArm,ChrEnd,ArmEnd) %>% filter(ClusterCount==1) %>%
             group_by(SampleId,ChrArm,Type,Id=ClusterId) %>% count()) %>%
       group_by(SampleId,ChrArm,Type,Id) %>% count() %>% group_by(SampleId,Type,ChrArm) %>% count() %>% spread(Type,n) %>% filter(ChrArm!='0_P'))

write.csv(svArmStats, "~/logs/sv_arm_stats.csv", row.names = F, quote = F)

nrow(svArmStats %>% filter(BndCount==0&SglCount==0))
nrow(svArmStats %>% filter(BndCount==0&SglCount==0&ConsistencyTotal==0))
nrow(svArmStats %>% filter(BndCount==0&SglCount==0&MaxClusterCount==1))

View(svArmStats %>% filter(BndCount==0&SglCount==0&ConsistencyTotal==0&MaxClusterCount>1))




View(svArmSummary %>% filter(SampleId=="CPCT02010337T"&Chr=="15"&Arm=="Q") %>% group_by(Chr,Arm,ClusterCount) %>% summarise(Clusters=n(), SvCount=sum(SvCount)))

View(svData %>% filter(SampleId=="CPCT02010359T"&ChrStart=="2"&ArmStart=="Q"))
View(svData %>% filter(ClusterCount==1205,Type=='DEL') %>% mutate(len=PosEnd-PosStart) %>% select(len,everything()))# %>% group_by(ChrStart,ArmStart,Type) %>% count()%>% spread(Type,n))
View(svData %>% filter(ClusterCount!=1205,Type=='DEL',SampleId=='CPCT02010359T') %>% mutate(len=PosEnd-PosStart) %>% select(len,everything()))

View(svData %>% filter(SampleId=='CPCT02010359TII') %>% mutate(len=PosEnd-PosStart) %>% select(len,everything()))

View(svData %>% filter(SampleId=='CPCT02010359TII'&ChrStart=='6'&PosStart>=11871212&PosEnd>0&PosEnd<=11872000))
View(svData %>% filter(ClusterCount==1,Type=='INV',Ploidy>0.15,AdjCNChgStart>0.15|AdjCNChgEnd>0.15) %>% mutate(len=PosEnd-PosStart)  %>% select(len,everything()))
View(svData %>% filter(Ploidy>0.15,AdjCNChgStart>0.15|AdjCNChgEnd>0.15) %>% group_by(cc=pmin(ClusterCount,20),Type) %>% count() %>% spread(Type,n))
View(svData %>% group_by(SampleId,ChrStart,ArmStart))

nrow(svData %>% filter(SampleId=="CPCT02020351T"&Type=='DUP'))
nrow(svData %>% filter(SampleId=="CPCT02020351T"&Type=='DUP'&ClusterCount>1))



# CHAIN ANALYSIS
chainData = svData %>% filter(ChainCount>0)
chainData$IsFoldback = grepl(';',chainData$ChainIndex)
chainData$FoldbackCount = str_count(chainData$ChainIndex,";")
View(chainData)


chainStats = (chainData %>% filter(ChainCount>1) %>% group_by(SampleId,ClusterId,ChainId)
              %>% summarise(SvCount=n(),
                            ChainCount=first(ChainCount),
                            FoldbackCount=sum(IsFoldback),
                            BndCount=sum(Type=='BND'),
                            ArmCount=first(ArmCount)))

View(chainStats)

chainStats$ChainLen = 2**round(log(chainStats$ChainCount,2))
chainStats$HasAmplification = (chainStats$FoldbackCount>0)


View(chainStats %>% group_by(ChainLen,HasAmplification) %>% count())

View(svData %>% filter(SampleId=="CPCT02020380T"&ClusterId==166&ChainId==1))
View(svData %>% filter(SampleId=="CPCT02020380T"&ClusterId==166&ChainId==1&Type=='BND'&(AsmbMatchStart=='MATCH'|AsmbMatchEnd=='MATCH')))


View(svData %>% filter(SampleId=="CPCT02010450T"&ChrStart==2))
View(svData %>% filter(SampleId=="CPCT02010450T"&ChrStart==2&ClusterCount==1))


View(svData %>% filter(SampleId=="CPCT02010450T"&ChrStart==2&ClusterCount>4&ClusterCount<100&PosStart>7800000&PosEnd<73000000))


simpleChains = chainData %>% filter()






# Single BND Analysis
singleBNDs = svData %>% filter(Type=='BND'&ClusterCount==1)
nrow(singleBNDs)
nrow(svData %>% filter(Type=='SGL'&ClusterCount==1))
nrow(svData %>% filter(Type=='NONE'&ClusterCount==1))

View(singleBNDs)

# check whether the arms have any complex activity
complexArms = svData %>% filter(ClusterSize=='Large') %>% group_by(SampleId,ChrStart,ArmStart) %>% summarise(Count=n())
View(complexArms)

singleBNDsArmInfo = merge(singleBNDs, complexArms, by.x=c('SampleId','ChrStart','ArmStart'), by.y=c('SampleId','ChrStart','ArmStart'), all.x=T)
names(singleBNDsArmInfo)[names(singleBNDsArmInfo) == 'Count'] <- 'StartArmCount'

singleBNDsArmInfo = merge(singleBNDsArmInfo, complexArms, by.x=c('SampleId','ChrEnd','ArmEnd'), by.y=c('SampleId','ChrStart','ArmStart'), all.x=T)
names(singleBNDsArmInfo)[names(singleBNDsArmInfo) == 'Count'] <- 'EndArmCount'

allStartArms = svData %>% group_by(SampleId,ChrStart,ArmStart) %>% summarise(Count=n())
allEndArms = svData %>% group_by(SampleId,ChrEnd,ArmEnd) %>% summarise(Count=n())

singleBNDsArmInfo = merge(singleBNDsArmInfo, allStartArms, by.x=c('SampleId','ChrStart','ArmStart'), by.y=c('SampleId','ChrStart','ArmStart'), all.x=T)
names(singleBNDsArmInfo)[names(singleBNDsArmInfo) == 'Count'] <- 'StartAllArmCount'
singleBNDsArmInfo = merge(singleBNDsArmInfo, allEndArms, by.x=c('SampleId','ChrEnd','ArmEnd'), by.y=c('SampleId','ChrEnd','ArmEnd'), all.x=T)
names(singleBNDsArmInfo)[names(singleBNDsArmInfo) == 'Count'] <- 'EndAllArmCount'

singleBNDsArmInfo[is.na(singleBNDsArmInfo)] <- 0
View(singleBNDsArmInfo)

singleBNDsArmInfo$PloidyBucket = ifelse(singleBNDsArmInfo$Ploidy < 1,0.2*round(singleBNDsArmInfo$Ploidy/0.2),1)
singleBNDsArmInfo$IsComplexArm = (singleBNDsArmInfo$StartArmCount>0|singleBNDsArmInfo$EndArmCount>0)
singleBNDsArmInfo$IsOnlySV = (singleBNDsArmInfo$StartAllArmCount<=1|singleBNDsArmInfo$EndAllArmCount<=0)

View(singleBNDsArmInfo %>% group_by(IsComplexArm,IsOnlySV,PloidyBucket) %>% count())

write.csv(singleBNDsArmInfo, "~/logs/sv_single_bnds.csv", row.names = F, quote = F)


View(svData %>% filter(SampleId=="CPCT02370013TII"&ChrStart==6&ArmStart=='P'))
View(svData %>% filter(SampleId=="CPCT02020186TII"&ChrEnd==3&ArmEnd=='P'))



# INV characteristics
invData = svData %>% filter(Type=='INV')
View(invData %>% group_by(ResolvedType) %>% count())

View(invData %>% group_by(ResolvedType,IsTI,ClusterSize) %>% count())

View(invData %>% filter(ResolvedType=='SimpleChain'&!IsTI))

View(svData %>% filter(SampleId=="CPCT02020245T"&ClusterId==20))
View(svData %>% filter(SampleId=="CPCT02020245T"&ChrStart==8))



# Reciprocal Translocation from 2 DUPs

recipTransData = svData %>% filter(ResolvedType=='RecipTrans')
View(recipTransData)

recipTrans = (recipTransData %>% arrange(SampleId,ClusterId,PosStart) %>% group_by(SampleId,ClusterId)
              %>% summarise(FirstPosStart=first(PosStart),
                            LastPosStart=last(PosStart),
                            FirstOrientStart=first(OrientStart),
                            LastOrientStart=last(OrientStart),
                            StartLinkType=ifelse(first(OrientStart)<last(OrientStart),'TI','DB'),
                            StartLinkLen=abs(first(PosStart)-last(PosStart)),
                            FirstPosEnd=first(PosEnd),
                            LastPosEnd=last(PosEnd),
                            FirstOrientEnd=first(OrientEnd),
                            LastOrientEnd=last(OrientEnd),
                            EndLinkType=ifelse((first(PosEnd)>last(PosEnd))==(first(OrientEnd)>last(OrientEnd)),'TI','DB'),
                            EndLinkLen=abs(first(PosEnd)-last(PosEnd))))

View(recipTrans)

write.csv(recipTrans, "~/logs/sv_recip_trans.csv", row.names = F, quote = F)


View(recipTrans %>% filter(StartLinkType=='TI'&EndLinkType=='TI'))

recipTransStarts = recipTrans %>% select


recipTrans$StartLen = ifelse()


rtDups = recipTrans %>% filter(LnkTypeStart=='TI'&LnkTypeEnd=='TI')
View(rtDups)

