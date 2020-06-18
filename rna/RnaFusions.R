

svFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191205.csv')
View(svFusions)

# initial 630 RNA samples vs LINX comparison
rnaLinxComp = read.csv('~/data/rna/fusions/LINX_dna_rna_combined_data_20191211.csv')
View(rnaLinxComp)

rnaLinxFusions = read.csv('~/data/rna/fusions/LNX_FUSIONS.csv')
View(rnaLinxFusions)
View(rnaLinxFusions %>% group_by(SampleId) %>% summarise(KnownPom=sum(KnownType!=''),
                                                         Reported=sum(Reportable=='true')) %>% filter(KnownPom>0|Reported>0))

View(rnaLinxFusions %>% filter(ClusterCount==1&TypeUp=='DUP'&RegionTypeUp=='Intronic'&CodingTypeUp=='Coding') %>% 
       mutate(Distance=abs(PosUp-PosDown)) %>%
       select(GeneNameUp,GeneNameDown,Distance,ChrUp,PosUp,PosDown,everything()))



dnaLinxSamples = rnaLinxFusions %>% group_by(SampleId) %>% summarise(KnownPom=sum(KnownType!=''),Reported=sum(Reportable=='true'))
dnaLinxSamples = merge(dnaLinxSamples,rnaSampleData %>% select(SampleId,ReadLength),by='SampleId',all.x=T)
nrow(dnaLinxSamples)
View(dnaLinxSamples)
write.csv(dnaLinxSamples %>% filter(KnownPom>0|Reported>0) %>% select(SampleId,ReadLength), '~/data/rna/samples/dna_fusion_samples.csv',row.names = F, quote = F)



## RNA prod fusion genes
rnaSampleFusions = read.csv('~/data/rna/fusions/rna_sample_fusion_genes.csv')
View(rnaSampleFusions)
View(rnaSampleFusions %>% group_by(SampleId) %>% count)
View(rnaSampleFusions %>% group_by(GeneNameUp) %>% count)
rnaSampleGenes = rbind(rnaSampleFusions %>% select(GeneName=GeneNameUp),rnaSampleFusions %>% select(GeneName=GeneNameDown))
View(rnaSampleGenes %>% group_by(GeneName) %>% count)
rnaSampleGenes = rnaSampleGenes %>% group_by(GeneName) %>% count
rnaSampleGenes = merge(ensemblGeneData %>% select(GeneId,GeneName),rnaSampleGenes,by='GeneName',all.y=T)
View(rnaSampleGenes)
write.csv(rnaSampleGenes %>% select(GeneId,GeneName), '~/data/rna/dna_fusion_gene_ids.csv',row.names = F, quote = F)

## Fusions and Chimeric Reads

sampleId = 'CPCT02020378T'

View(chimReads %>% filter(PosStart==42870046|PosEnd==39795483))
                            
starFusion = read.csv('~/data/rna/fusions/LNX_RNA_DATA.csv')
View(starFusion)

rnaFusions = read.csv(formFilename('~/logs/',sampleId,'fusions.csv'))
View(rnaFusions)

View(rnaFusions %>% filter(Valid=='true'&NonSupp=='true') %>% 
       select(NonSupp,GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,SplitFrags,RealignedFrags,DiscordantFrags,everything()))

# all
sampleId = 'CPCT02020378T'
allRnaFusions = read.csv(formFilename('~/data/rna/runs/',sampleId,'fusions.csv'))
nrow(allRnaFusions)
View(allRnaFusions)

View(allRnaFusions %>% filter(Valid=='true'&SplitFrags>=3&JuncTypeUp %in% c('KNOWN','CANONICAL')&JuncTypeDown %in% c('KNOWN','CANONICAL')))
View(allRnaFusions %>% filter(Valid=='true'&SplitFrags>=10))
View(allRnaFusions %>% filter(Valid=='true'&NonSupp=='true'))
View(allRnaFusions %>% filter(Valid=='true'&SVType=='DEL'&as.character(GeneNameUp)==as.character(GeneNameDown)))
View(allRnaFusions %>% filter(SVType=='INV'&SplitFrags>1&PosUp-PosDown<1e3))
View(allRnaFusions %>% filter(SVType=='INV'&PosUp-PosDown<1e3))
View(allRnaFusions %>% filter(SVType=='DEL'&PosUp-PosDown<5e5))




# duplicates
View(allRnaFusions %>% group_by(ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count() %>% filter(n>1))

# counts by chr pair
View(allRnaFusions %>% group_by(ChrUp,ChrDown) %>% count())
print(sum(allRnaFusions$TotalFragments))

View(allRnaFusions %>% filter(Valid=='true'&NonSupp=='true') %>% 
       select(GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,SplitFrags,RealignedFrags,DiscordantFrags,
              CoverageUp,CoverageDown,everything()))

View(allRnaFusions %>% filter(Valid=='true'&GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG') %>% 
       select(NonSupp,GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,SplitFrags,RealignedFrags,DiscordantFrags,
              CoverageUp,CoverageDown,everything()))


View(rnaFusions %>% filter(SplitFrags>1))
View(rnaFusions %>% group_by(GeneNameUp) %>% count %>% filter(n>5))
View(rnaFusions %>% filter(SVType=='INV'&PosUp==PosDown))
View(rnaFusions %>% filter(SVType=='DEL'&PosUp-PosDown<5e5))
View(rnaFusions %>% group_by(ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count() %>% filter(n>1))

# restricted genes
rnaFusions = read.csv(formFilename('~/logs/',sampleId,'fusions.csv'))

View(rnaFusions %>% select(GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,
                           SplitFrags,RealignedFrags,DiscordantFrags,CoverageUp,CoverageDown,everything()))
View(rnaFusions)

chimFrags = read.csv(formFilename('~/logs/',sampleId,'chimeric_frags.csv'))
View(chimFrags) 

chimReads = read.csv(formFilename('~/logs/',sampleId,'chimeric_reads.csv'))
View(chimReads) # 2272

View(chimFrags %>% group_by(FusionGroup,Type) %>% count %>% spread(Type,n,fill=0))
View(chimFrags %>% group_by(Type) %>% count)
View(chimFrags %>% group_by(Type,ReadCount) %>% count)

# specific chromosome
sampleId = 'CPCT02020378T'
rnaFusionsChr = read.csv(formFilename('~/logs/',sampleId,'spec_chr.fusions.csv'))
View(rnaFusionsChr)
View(rnaFusionsChr %>% group_by(ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count() %>% filter(n>1))
View(rnaFusionsChr %>% filter(SVType=='INV'))
View(rnaFusionsChr %>% filter(SVType=='INV'&PosUp==PosDown))
View(rnaFusionsChr %>% filter(SVType=='DEL'&PosUp-PosDown<5e5))

chimReadsChr = read.csv(formFilename('~/logs/',sampleId,'spec_chr.chimeric_reads.csv'))
View(chimReadsChr)
View(chimReadsChr %>% group_by(FusionGroup) %>% count)

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000233864'))
View(ensemblTransExonData %>% filter(TransName=='ENST00000518476'))
View(ensemblGeneData %>% filter(Chromosome==8&GeneStart>88145481-1e5&GeneEnd>88145481+1e5))


chimFragsChr = read.csv(formFilename('~/logs/',sampleId,'spec_chr.chimeric_frags.csv'))
View(chimFragsChr)
View(chimFragsChr %>% filter(ChrStart!=ChrEnd))

chimReadsAll = read.csv(formFilename('~/logs/',sampleId,'all.chimeric_reads.csv'))
View(tmpChm)

View(ensemblGeneData %>% filter(GeneName=='MIR3648'|GeneName=='MIR3687'))

chimReads = merge(chimReads, chimFrags %>% select(ReadId,Type,SameGene),by='ReadId',all.x=T)
View(chimReads %>% filter(FusionGroup!='SINGLE_GENE') %>% select(ReadId,Type,everything()))

View(chimReads %>% group_by(FusionGroup) %>% count)
View(chimReads %>% filter(FirstInPair==NegStrand))

chimReads2 = read.csv(formFilename('~/logs/',sampleId,'chimeric_reads2.csv'))
View(chimReads2)

discordantFrags = chimFrags %>% filter(Type=='DISCORDANT')
View(chimReads %>% filter(ReadId %in% discordantFrags$ReadId))

max(ensemblGeneData$GeneEnd - ensemblGeneData$GeneStart)
mean(ensemblGeneData$GeneEnd - ensemblGeneData$GeneStart)
median(ensemblGeneData$GeneEnd - ensemblGeneData$GeneStart)

dsRnaFusions = read.csv(formFilename('~/data/rna/runs/',sampleId,'ds.fusions.csv'))
View(dsRnaFusions)

View(dsRnaFusions %>% filter(MultiMapFrags>0) %>% select(GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,
                                                         SplitFrags,RealignedFrags,DiscordantFrags,MultiMapFrags,CoverageUp,CoverageDown,everything()))


# short and local DELs, DUPs and INVs
sampleId = 'CPCT02440010T'
sampleId = 'CPCT02060066T'
rnaFusions = read.csv(formFilename('~/logs/',sampleId,'fusions.csv'))
View(rnaFusions)

chimFrags = read.csv(formFilename('~/logs/',sampleId,'chimeric_frags.csv'))
View(chimFrags) 

chimReads = read.csv(formFilename('~/logs/',sampleId,'chimeric_reads.csv'))
View(chimReads) 
View(chimReads %>% filter(FusionGroup=='Id_0'&Supplementary=='false') %>% group_by(ReadId) %>% 
       summarise(Count=n(),
                 Orient1=ifelse((first(FirstInPair)=='true')==(first(ReadReversed)=='false'), 1,-1),
                 Orient2=ifelse((last(FirstInPair)=='true')==(last(ReadReversed)=='false'), 1,-1))) 

View(chimReads %>% group_by(ReadId) %>% summarise(Count=n(),MinPosStart=min(PosStart),MaxPosEnd=max(PosEnd)) %>% mutate(ReadSpan=MaxPosEnd-MinPosStart))

readData = read.csv(formFilename('~/logs/',sampleId,'read_data.csv'))
readData = read.csv(formFilename('~/logs/','CPCT02060066T','read_data.csv'))

View(readData) 
rm(readData) 

View(geneOverlaps)

# short DEL in CPCT02440010T - RBBP4 (ENSG00000162521) 1:33M,123-290K and S100PBP (ENSG00000116497)
View(ensemblGeneData %>% filter(GeneId %in% c('ENSG00000162521','ENSG00000116497')))
View(ensemblTransExonData %>% filter(GeneId %in% c('ENSG00000162521','ENSG00000116497')))


# short INV in CPCT02440010T - ERBB2, PGAP3 - 17,842-868K
View(ensemblGeneData %>% filter(GeneName %in% c('ERBB2', 'PGAP3')))

# CPCT02010802T, CPCT02080141T - GOPC (ENSG00000047932) ROS1 (ENSG00000047936) - overlapping 
View(ensemblGeneData %>% filter(GeneId %in% c('ENSG00000162521','ENSG00000116497')))

# CPCT02080164T - SLC45A3 (ENSG00000158715) ELK4 (ENSG00000158711), interesting case as the SV DEL is only 71 bases long
View(ensemblGeneData %>% filter(GeneId %in% c('ENSG00000047932','ENSG00000047936'))) 


# for DUPs lots of examples of FGFR3 (ENSG00000068078) TACC3 (ENSG00000013810) in the database
# for INV please take a look at NAB2 (ENSG00000166886) STAT6 (ENSG00000166888) examples



exonDelDupSamples = read.csv('~/data/sv/fusions/exon_del_dup_samples.csv')
nrow(exonDelDupSamples)
View(exonDelDupSamples %>% group_by(SampleId) %>% count)

eddFusions = read.csv('~/logs/LNX_FUSIONS.csv')
View(eddFusions)
View(eddFusions %>% filter(GeneNameUp=='EGFR'&GeneNameDown=='EGFR'))
View(eddFusions %>% group_by(KnownType,Reportable) %>% count)

igFusions = read.csv('~/logs/LNX_FUSIONS.csv')
View(igFusions)

write.csv(igFusions %>% filter(KnownType=='IG_KNOWN_PAIR'| KnownType=='IG_PROMISCUOUS'),'~/logs/LNX_IG_FUSIONS.csv', row.names = F, quote = F)

View(igFusions %>% group_by(CodingTypeDown,RegionTypeDown) %>% count())
View(igFusions %>% filter(CodingTypeDown=='Coding'&RegionTypeDown=='Upstream'&OrientDown=='-1') %>% 
       select(CodingTypeDown,RegionTypeDown,PosDown,CodingStartDown,CodingEndDown,TotalCodingDown,CodingBasesDown,
              TransStartDown,TransEndDown,OrientDown,BreakendExonDown,ExonicPhaseDown,everything()))
colnames(igFusions)

cciaSVs = read.csv('~/logs/ccia_p007204_svs.csv')
View(cciaSVs)

View(ensemblTransData %>% filter(GeneId=='ENSG00000173327')) # MAP3K11
View(ensemblTransData %>% filter(GeneId=='ENSG00000173039')) # RELA

View(ensemblGeneData %>% filter(Chromosome==2&Strand==1&GeneStart>=89890568&GeneEnd<=90274235))
View(ensemblGeneData %>% filter(Chromosome==14&Strand==-1&GeneStart>=106032614&GeneEnd<=107288051))
View(ensemblGeneData %>% filter(Chromosome==22&Strand==1&GeneStart>=22380474&GeneEnd<=23265085))


# Prod Linx 5' genes for RNA expression
upstreamFusionGenes = newFusions %>% select(SampleId,GeneIdUp,GeneNameUp)
View(upstreamFusionGenes %>% group_by(SampleId,GeneIdUp,GeneNameUp) %>% count)
upstreamFusionGenes = upstreamFusionGenes %>% filter(SampleId %in% rnaSampleIds_2131$SampleId)
View(upstreamFusionGenes %>% group_by(SampleId) %>% count)

write.csv(upstreamFusionGenes %>% group_by(SampleId,GeneIdUp,GeneNameUp) %>% mutate(CancerType='Unknown') %>% ungroup() %>% select(SampleId,CancerType,GeneId=GeneIdUp),
          '~/data/rna/fusions/linx_upstream_sample_gene_ids.csv',row.names = F,quote = F)

nrow(rnaSampleData2131)


# Sample investigations
specSampleId = 'CPCT02010386T'
specSampleId = 'CPCT02020310T'
specSampleId = 'CPCT02020380T'
specSampleId = 'CPCT02140038T'
specSampleId = 'CPCT02020419T'
sampleFusions = read.csv(formFilename('~/logs/',specSampleId,'fusions.csv'))
View(sampleFusions)

sampleFrags = read.csv(formFilename('~/logs/',specSampleId,'chimeric_frags.csv'))
View(sampleFrags)

sampleReads = read.csv(formFilename('~/logs/',specSampleId,'chimeric_reads.csv'))
View(sampleReads)
View(sampleReads %>% filter(ReadId %in% 
                              c('NB500901:87:HYLCHBGX7:1:11210:17290:2382',
                                'NB500901:87:HYLCHBGX7:4:11611:16992:3375',
                                'NB500901:87:HYLCHBGX7:3:22503:17184:13635',
                                'NB500901:87:HYLCHBGX7:1:13206:19326:9633')))


View(rnaFusionsChr %>% group_by(ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count() %>% filter(n>1))
View(rnaFusionsChr %>% filter(SVType=='INV'))
View(rnaFusionsChr %>% filter(SVType=='INV'&PosUp==PosDown))
View(rnaFusionsChr %>% filter(SVType=='DEL'&PosUp-PosDown<5e5))

arribaTmp=read.csv('~/data/rna/runs/CPCT02020684T.fusions.tsv',sep='\t')
View(arribaTmp)

chimReadsChr = read.csv(formFilename('~/logs/',sampleId,'spec_chr.chimeric_reads.csv'))
View(chimReadsChr)
View(chimReadsChr %>% group_by(FusionGroup) %>% count)

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000233864'))
View(ensemblTransExonData %>% filter(TransName=='ENST00000518476'))
View(ensemblGeneData %>% filter(Chromosome==8&GeneStart>88145481-1e5&GeneEnd>88145481+1e5))


chimFragsChr = read.csv(formFilename('~/logs/',sampleId,'spec_chr.chimeric_frags.csv'))
View(chimFragsChr)
View(chimFragsChr %>% filter(ChrStart!=ChrEnd))





# CPCT02020378T matching
View(starFusion %>% filter(SampleId=='CPCT02020378T'&FusionName=='TMPRSS2_ERG'))



View(ensemblTransExonData %>% filter(TransName=='ENST00000332149')) # TMP
View(ensemblTransExonData %>% filter(TransName=='ENST00000497881')) # TMP
View(ensemblTransExonData %>% filter(TransName=='ENST00000398897')) # ERG
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000184012')) # TMPRSS2
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000157554')) # ERG
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000175567')) 
View(ensemblTransExonData %>% filter(TransName=='ENST00000341446'))



# NCOR1 - FO



knownPairs = read.csv('~/data/sv/knownFusionPairs.csv')
View(knownPairs)

kpBedInfo = merge(knownPairs,ensemblGeneData %>% select(GeneName,UpChr=Chromosome,UpStrand=Strand,UpGeneStart=GeneStart,UpGeneEnd=GeneEnd),by.x='fiveGene',by.y='GeneName',all.x=T)
kpBedInfo = merge(kpBedInfo,ensemblGeneData %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),by.x='threeGene',by.y='GeneName',all.x=T)
View(kpBedInfo)

preGeneBuffer=10e3
kpBedInfo = kpBedInfo %>% mutate(Name=paste(fiveGene,threeGene,sep='-'),
                                 Start1=ifelse(UpStrand==1,UpGeneStart-preGeneBuffer,UpGeneStart),
                                 End1=ifelse(UpStrand==1,UpGeneEnd,UpGeneEnd+preGeneBuffer),
                                 Start2=ifelse(DownStrand==1,DownGeneStart-preGeneBuffer,DownGeneStart),
                                 End2=ifelse(DownStrand==1,DownGeneEnd,DownGeneEnd+preGeneBuffer),
                                 Strand1=ifelse(UpStrand==1,'+','-'),
                                 Strand2=ifelse(DownStrand==1,'+','-'),
                                 Info='.')
View(kpBedInfo)

write.csv(kpBedInfo %>% select(UpChr,Start1,End1,DownChr,Start2,End2,Name,Strand1,Strand2,Info),'~/data/sv/known_pairs_bedpe',row.names = F,quote = F)






