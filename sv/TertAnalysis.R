


svData = read.csv('~/data/sv/fusions/LNX_SVS.csv') # want to include genes 10K upstream
clusters = read.csv('~/data/sv/drivers/LNX_CLUSTERS.csv')
drivers = read.csv('~/data/sv/drivers/LNX_DRIVERS.csv')
fusions = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv')

View(clusters)
View(svData %>% filter(Type=='INF'|Type=='SGL'))

View(fusions %>% filter(GeneNameDown=='TERT'))


# BND coming into TERT promotor region for CNS

tertBnds = svData %>% filter(Type=='BND'&ChrStart==5&ChrEnd==10&PosStart>1e6&PosStart<1.25e6)
View(tertBnds)

# TERT: chr5 strand-1, pos: 1,253,262 to 1,295,184



oncoGenes = read.csv('~/data/oncogenes.csv')
sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
oncogeneBreakends = read.csv('~/data/sv/oncogene_breakends.csv')

oncogeneBreakends = oncogeneBreakends %>% mutate(CancerType=ifelse(CancerType=='-1','Unknown',as.character(CancerType)),
                                                 OrientEnd=ifelse(Type=='INF'|Type=='SGL',0,OrientEnd))
View(oncogeneBreakends %>% group_by(CancerType) %>% count)
View(oncogeneBreakends)
View(oncogeneBreakends %>% group_by(GeneName) %>% count)


oncogeneBreakends = merge(oncogeneBreakends,
                          svData %>% select(SampleId,ClusterId,ClusterCount,ResolvedType,Type,PosStart,PosEnd,OrientStart,OrientEnd),
                          by=c('SampleId','Type','PosStart','PosEnd','OrientStart','OrientEnd'),all.x=T)
oncogeneBreakends = oncogeneBreakends %>% filter(!is.na(ClusterId))
nrow(oncogeneBreakends %>% filter(is.na(ClusterId)))
View(oncogeneBreakends)

oncoCopyNumber = read.csv('~/data/sv/prod_oncogene_cn.csv',sep='\t')
View(oncoCopyNumber)
oncogeneSummary = merge(oncogeneBreakends,oncoCopyNumber,by=c('SampleId','GeneName'),all.x=T)
View(oncogeneSummary)
View(oncoCopyNumber %>% group_by(GeneName) %>% count)

ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)
oncogeneSummary = merge(oncogeneSummary,ensemblGeneData %>% select(GeneName,GeneStart,GeneEnd,Strand),by='GeneName',all.x=T)
View(oncogeneSummary)
nrow(oncogeneSummary)

oncogeneSummary = oncogeneSummary %>% filter(MinCopyNumber/SamplePloidy<1.5)

oncogeneSummary = oncogeneSummary %>% mutate(BePosition=ifelse(IsStart==1,PosStart,PosEnd),
                                             OtherBePosition=ifelse(IsStart==1,PosEnd,PosStart),
                                             PreGeneDistance=ifelse(Strand==1,GeneStart-BePosition,BePosition-GeneEnd))

View(oncogeneSummary)
View(oncogeneSummary %>% filter(GeneName=='TERT'))
oncogeneSummary = oncogeneSummary %>% filter(PreGeneDistance>=0)
colnames(oncogeneSummary)

# take the closest breakend per sample to the gene
oncogeneSampleData = oncogeneSummary %>% arrange(SampleId,GeneName,PreGeneDistance) %>% group_by(SampleId,GeneName) %>%
  summarise(CancerType=first(CancerType),
            Type=first(Type),
            BePosition=first(BePosition),
            OtherBePosition=first(OtherBePosition),
            ClusterCount=first(ClusterCount),
            ResolvedType=first(ResolvedType),
            PreGeneDistance=first(PreGeneDistance))

oncogeneSampleData = oncogeneSampleData %>% mutate(PreGeneBucket=round(PreGeneDistance,-4))

View(oncogeneSampleData)
View(oncogeneSampleData %>% group_by(GeneName,PreGeneBucket) %>% summarise(SampleCount=n()) %>% spread(PreGeneBucket,SampleCount))
View(oncogeneSampleData %>% group_by(GeneName,CancerType) %>% summarise(SampleCount=n()) %>% spread(CancerType,SampleCount))
View(oncogeneSampleData %>% group_by(GeneName,Type) %>% summarise(SampleCount=n()) %>% spread(Type,SampleCount))
View(oncogeneSampleData %>% group_by(GeneName,ResolvedType) %>% summarise(SampleCount=n()) %>% spread(ResolvedType,SampleCount))

View(oncogeneSampleData %>% filter(GeneName=='TERT') %>% group_by(CancerType,PreGeneBucket) %>% summarise(SampleCount=n()) %>% spread(PreGeneBucket,SampleCount))

foxa1Drivers = read.csv('~/data/sv/foxa1_drivers.csv',sep='\t')
View(foxa1Drivers)

View(oncogeneSampleData %>% filter(GeneName=='TERT'&PreGeneDistance<25e3))
View(oncogeneSampleData %>% filter(GeneName=='FOXA1'))

oncogeneSampleData = merge(oncogeneSampleData,foxa1Drivers,by=c('SampleId','GeneName'),all.x=T)

View(drivers)


View(svData %>% filter(grepl('HIST1H3B',GeneStart)|grepl('HIST1H3B',GeneEnd)))


View()


tertChr=5
tertStart=1253262
tertEnd=1295184
tertProx=5e5

tertSVs = svData %>% filter((ChrStart==tertChr&PosStart>tertStart&PosStart<tertEnd+tertProx)|
                              (ChrEnd==tertChr&PosEnd>tertStart&PosEnd<tertEnd+tertProx))


tertSVs = tertSVs %>% mutate(StartInGene=grepl('TERT',GeneStart),
                             EndInGene=grepl('TERT',GeneEnd),
                             StartPreGene=(ChrStart==tertChr&PosStart>tertEnd&PosStart<tertEnd+tertProx),
                             EndPreGene=(ChrEnd==tertChr&PosEnd>tertEnd&PosEnd<tertEnd+tertProx),
                             PreGeneLength=ifelse(StartPreGene&EndPreGene,pmin(PosStart,PosEnd)-tertEnd,
                                                  ifelse(StartPreGene,PosStart-tertEnd,ifelse(EndPreGene,PosEnd-tertEnd,0))),
                             PreGeneLengthBucket=2**round(log(PreGeneLength,2)))

View(tertSVs %>% group_by(Type,StartInGene,EndInGene,StartPreGene,EndPreGene) %>% count)
View(tertSVs %>% group_by(Type,PreGeneLengthBucket) %>% count)
View(tertSVs %>% group_by(ResolvedType) %>% count)


print(ggplot(data=tertSVs %>% group_by(PreGeneLengthBucket,Type) %>% count(), aes(x=PreGeneLengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Type))

print(ggplot(data=tertSVs %>% filter(ResolvedType %in% c('COMPLEX','DEL','DUP','RECIP_TRANS')) %>% 
               group_by(PreGeneLengthBucket,ResolvedType) %>% count(), aes(x=PreGeneLengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ResolvedType))


View(fusions %>% filter(GeneNameDown=='TERT') %>% select(SampleId,GeneNameUp,DistancePrevDown,BreakendExonDown,PhaseDown,RegionTypeDown,
                                                         PhaseMatched,ChrUp,PosUp,TypeUp,ResolvedType,everything()))

colnames(fusions)





