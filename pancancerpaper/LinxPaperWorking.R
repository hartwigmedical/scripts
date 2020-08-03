
#### LINE


View(linePsdgSampleData)

write.csv(linePsdgSampleSummary,'~/logs/linePsdgSampleSummary.csv',row.names = F,quote = F)

# 'Count of LINE Break Junctions'
# add n = sample count

linePsdgSampleSummary = linePsdgSampleSummary %>% group_by(CancerType) %>% mutate(MedianLineBE=median(LineBEs),
                                                                                  PseudogeneTotal=sum(Pseudogenes)) %>% ungroup()
linePsdgSampleSummary = merge(linePsdgSampleSummary,cancerTypeTotals,by='CancerType',all.x=T)
View(linePsdgSampleSummary)

print(ggplot(linePsdgSampleSummary %>% mutate(LineBEs=LineBEs+0.5,CancerTypeLabel=sprintf('%s (n=%d)',CancerType,SampleCount)), 
             aes(x=reorder(CancerTypeLabel,-MedianLineBE),y=LineBEs))
      + geom_violin(scale="area",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "LINE Break Junctions per Sample",y='LINE Break Junctions',x=''))

print(ggplot(linePsdgSampleSummary %>% mutate(Pseudogenes=Pseudogenes+0.5,CancerTypeLabel=sprintf('%s (n=%d)',CancerType,SampleCount)), 
             aes(x=reorder(CancerTypeLabel,-PseudogeneTotal),y=Pseudogenes))
      + geom_violin(scale="area",fill="#6baed6")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "Pseudogenes per Sample",y='Pseudogenes',x=''))

## LINE by type and location
# CIRCOS showing count of breakends in LINE source elements and insertion sites by 100k bucket.   Known LINE source regions in different colour and novel in a 3rd colour.
# Table of most activated regions? Known vs novel?? 




#####
## BFB without foldbacks
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
svData = svData %>% filter(ClusterCount>1)
clusters = read.csv('~/data/sv/cohort/LNX_CLUSTERS.csv')
clusters = clusters %>% filter(ClusterCount>1)

bfbnf = clusters %>% filter(Foldbacks==0&Replication=='true'&ClusterCount<=5&ResolvedType=='COMPLEX')
nrow(bfbnf)
View(bfbnf)

bfbnfChained = bfbnf %>% filter(FullyChained=='true'&ChainCount==1)
nrow(bfbnfChained)

View(bfbnfChained %>% filter(BndCount>=1&SglCount==0&InfCount==0) %>%
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,MinJcn,MaxJcn,AssemblyTIs,OriginArms,FragmentArms,everything()))

write.csv(bfbnfChained %>% filter(BndCount>=1&SglCount==0&InfCount==0) %>% filter(MaxJcn==2) %>% select(SampleId,ClusterId),
          '~/logs/bfbnf_ids.csv',row.names = F,quote = F)

View(bfbnfChained %>% filter(BndCount>=2&SglCount==0&InfCount==0&ClusterCount==3) %>%
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,MinJcn,MaxJcn,AssemblyTIs,OriginArms,FragmentArms,everything()))

View(bfbnf %>% filter(ResolvedType=='COMPLEX'&FullyChained=='true'&ChainCount==1&BndCount>=2))

      



## cluster size vs arms involved by resolved type
View(clusters)


armSvCounts = clusters %>% group_by(SampleId,ClusterType)
sampleClusterTypes = clusters %>% group_by(SampleId,ClusterType) %>% summarise(ClusterTypeCount=sum(ClusterCount))

print(ggplot(clusters %>% filter(ClusterType=='LINE'|ClusterType=='COMPLEX'), aes(x=ClusterCount,y=ArmCount))
      + geom_point()
      + scale_x_log10()
      + facet_wrap(~ClusterType)
      + labs(title = "SV Count vs Arm Count by Cluster Type"))




#####
## Complex Clusters

complexClusters = clusters %>% filter(ResolvedType=='COMPLEX'&ClusterCount>=3)
complexClusters = complexClusters %>% mutate(ClusterSize=ifelse(ClusterCount<=4,'3-4',ifelse(ClusterCount<=10,'5-10',ifelse(ClusterCount<=20,'11-20','>50'))))

ccSummary = complexClusters %>% group_by(SampleId,ClusterSize) %>% summarise(CountOfClusters=n(),ClusterSizeN=min(ClusterCount)) %>% 
  group_by(ClusterSize,CountOfClusters=pmin(CountOfClusters,20)) %>% summarise(Count=n(),ClusterSizeN=first(ClusterSizeN))

View(ccSummary)

print(ggplot(ccSummary) +
        geom_bar(aes(x=CountOfClusters,y=Count,fill=reorder(ClusterSize,-ClusterSizeN)),stat='identity',alpha=1) +
        labs(title='Frequency of Complex Clusters by Size',x='Clusters per Sample',y='# Clusters',fill='SVs per Cluster'))


knowFusionData = read.csv('~/data/sv/known_fusion_data.csv')
View(knowFusionData)

# 1. Frequency of complex events per sample
ccPerSample = complexClusters %>% group_by(SampleId) %>% summarise(SampleCC=n())
View(ccPerSample)

ccZeroSamples = samplesDD %>% filter(!(SampleId %in% ccPerSample$SampleId)) %>% mutate(SampleCC=0) %>% select(SampleId,SampleCC)
nrow(ccZeroSamples)
View(ccZeroSamples)
ccPerSample = rbind(ccPerSample,ccZeroSamples)

ccPerSample = ccPerSample %>% mutate(CCBucket=ifelse(SampleCC<=50,as.character(SampleCC),'>50'))

View(ccPerSample %>% group_by(CCBucket) %>% summarise(SampleCC=first(SampleCC),Count=n()))

print(ggplot(ccPerSample %>% group_by(CCBucket) %>% summarise(SampleCC=first(SampleCC),Count=n()), aes(x=reorder(CCBucket,SampleCC), y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Number of Complex Clusters per Sample', y='Sample Count',x=''))



ggplot(df) +
  geom_tile(aes(x = MajorAllele, y = MinorAllele,  color = Penalty)) +
  scale_colour_gradientn(colours=csColours, limits = c(0, 15)) +
  theme( panel.border = element_blank()) +
  xlab("Major Allele Ploidy") + ylab("Minor Allele Ploidy")







## LINE and Pseudogenes
print(ggplot(linePsdgSampleData %>% filter(Pseudogenes>=0&LineBucket>=0),aes(x=reorder(LineBucket2,LineBEs),y=Pseudogenes))
      + geom_point()
      # + scale_x_log10()
      + labs(title = 'LINE elements vs Pseudogenes per Sample', x = 'LINE Break Junctions', y = 'Pseudogenes'))

linePsdgSampleData = merge(linePsdgSampleData,samplesDD %>% select(SampleId,CancerType),by='SampleId',all.x=T)

linePsdgSampleSummary = linePsdgSampleData %>% group_by(CancerType) %>% mutate(MedianLine=median(LineBEs),
                                                                               MedianPseudogenes=median(Pseudogenes)) %>% ungroup()

View(linePsdgSampleSummary)

linePsdgColours = c("#a6611a","#80cdc1") # "#dfc27d", "#018571"
linePsdgColours = setNames(linePsdgColours, c('LINE Break Junctions','Pseudogenes'))
linePsdgLinetypes = c('solid','solid')
linePsdgLinetypes = setNames(linePsdgLinetypes, c('LINE Break Junctions','Pseudogenes'))

print(ggplot(linePsdgSampleSummary %>% 
               # filter(CancerType %in% c('Esophagus','Breast','Prostate')) %>%
               mutate(LineBEs=LineBEs+0.5,Pseudogenes=Pseudogenes+0.5) %>%
               arrange(-MedianLine) %>% mutate_at(vars(CancerType), funs(factor(., levels=unique(.)))))
      + stat_ecdf(size = 0.3, aes(LineBEs,color='LINE Break Junctions',linetype='LINE Break Junctions'), geom = "step", pad = FALSE) 
      + geom_segment(size = 0.3, aes(x = MedianLine, xend = MedianLine, y = 0.25, yend = 0.75, color='LINE Break Junctions'), show.legend = F)
      + stat_ecdf(size = 0.3, aes(Pseudogenes,color='Pseudogenes',linetype='Pseudogenes') ,geom = "step", pad = FALSE) 
      # + geom_segment(size = 0.3, aes(x = MedianPseudogenes, xend = MedianPseudogenes, y = 0.25, yend = 0.75, color='HMF MNV'), show.legend = F) +
      + scale_x_log10(labels = c('0','10', '100', '1000', '10000'), breaks = c(0, 10, 100, 1000, 10000))  
      + facet_grid(~CancerType) + theme(panel.spacing = unit(1, "pt"))
      + scale_colour_manual(name = "Combined", values=linePsdgColours)
      + scale_linetype_manual(name = "Combined", values = linePsdgLinetypes)
      + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
              axis.title.y = element_text(vjust = 0),
              # strip.text=element_text(vjust=10),
              strip.text.x = element_text(angle=90, hjust=1,size=8),
              # panel.border = element_rect(size = 0.3),
              panel.grid.minor.x = element_blank(),
              # panel.grid.major = element_blank(),
              panel.grid.major = element_line(size = 0.3),
              # panel.background = element_blank(), 
              strip.background = element_blank(), 
              # strip.text = element_blank(), 
              legend.position=c(0.5, 0.02), 
              legend.title = element_blank(), 
              legend.margin= margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"),
              plot.margin = margin(t = 3, b = 0, l = 3, r = 3, unit = "pt"),
              legend.background=element_blank(), legend.key=element_blank())
      + xlab("LINE Break Junctions vs Pseudogenes")
      + guides(colour = guide_legend(nrow = 1))
      + coord_flip())


plotIndex = 1
cancerPlots = list()

for(cancerType in c('Prostate','Uterus'))
  #for(cancerType in unique(sampleResolvedTypes$CancerType))
{
  plot = (ggplot(sampleResolvedTypes %>% filter(CancerType==cancerType),
                 aes(x=reorder(SampleLabel,-SampleSvTotal),y=ClusterTypePercent,fill=ClusterType))
          + geom_bar(stat="identity",width=1)
          + labs(title=cancerType, x='', y='')
          + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                  axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank()))
  
  cancerPlots[[plotIndex]] = plot
  plotIndex = plotIndex + 1
}

plot_grid(cancerPlots[[1]],cancerPlots[[2]],nrow=4,ncol=1)

plot_grid(cancerPlots[[1]],cancerPlots[[2]],cancerPlots[[3]],cancerPlots[[4]],nrow=4,ncol=1)


write.csv(ensemblGeneData %>% filter(GeneName %in% c('COL1A1','COL1A2','COL3A1','COL4A1','COL6A3','FN1','IGFBP7')) %>% select(GeneId,GeneName),
          '~/logs/col_gene_ids.csv',row.names = F,quote = F)

colExp = read.csv('~/data/sv/cohort/isofox_exp_col_genes.csv')
View(colExp)
View(colExp %>% filter(SampleId=='WIDE01010564T'|SampleId=='WIDE01010575T'))



View(homDisClusters)
View(homDisClusters %>% filter(SvMatchType=='DB'&ClusterCount==1))
write.csv(homDisClusters %>% filter(SvMatchType=='DB'&ClusterCount==1) %>% group_by(SampleId) %>% count %>% select(SampleId),
          '~/logs/hom_dis_del_sample_ids.csv',row.names = F,quote = F)

tmpDrivers = read.csv('~/logs/LNX_DRIVERS.csv')
View(tmpDrivers)

View(svData %>% filter(SampleId=='CPCT02010110TIII') %>% filter(ChrStart==2))


View(svData %>% filter(SampleId=='CPCT02340020T'))


clusterJcn = clusterJcn %>% mutate(ClusterSize=apply(clusterJcn[,'ClusterCount'],1,function(x) get_cluster_size(x[1])))

clusterJcn$ClusterSize = apply(clusterJcn,1,get_cluster_size,clusterCount=ClusterCount)
tmp = apply(clusterJcn,1,get_cluster_size,clusterCount=ClusterCount)

View(clusterJcn)

View(clusterJcn %>% filter(RelativeJcn<1) %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,RelativeJcn,MaxJcn,MinJcn,Ploidy,everything()))
View(clusterJcn %>% group_by(MaxJcn=ifelse(RelativeJcn<=0,0,2**round(log(RelativeJcn,2)))) %>% count)

get_cluster_size<-function(clusterCount)
{
  clusterSize = ''
  
  if(clusterCount <= 2)
    clusterSize = as.character(clusterCount)
  else if(clusterCount <= 5)
    clusterSize = '3-5'
  else if(clusterCount <= 10)
    clusterSize = '6-10'
  else if(clusterCount <= 25)
    clusterSize = '11-25'
  else
    clusterSize = '>25'
  
  return (clusterSize)
}

print(get_cluster_size(1))
print(get_cluster_size(4))
print(get_cluster_size(18))





topHomDisGenes = head(topDelHomDisGenes %>% arrange(-HomDisCount), 20) %>% select(GeneName=Gene)
topDelGenes = head(topDelHomDisGenes %>% arrange(-DelCount), 20) %>% select(GeneName=Gene)
View(topHomDisGenes)
View(topDelGenes)
topGenes = rbind(topHomDisGenes,topDelGenes %>% filter(!(GeneName %in% topHomDisGenes$GeneName)))
topGenes = merge(ensemblGeneData %>% select(GeneId,GeneName),topGenes,all.y=T)
View(topGenes)
write.csv(topGenes %>% select(GeneId,GeneName),'~/data/sv/cohort/top_del_hom_dis_gene_ids.csv',row.names = F,quote = F)
