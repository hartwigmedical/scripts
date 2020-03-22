

## Ensembl ref data
ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)

ensemblTransExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000181163'))

ensemblTransData = read.csv('~/data/sv/ensembl_trans_data.csv')
View(ensemblTransData)
View(ensemblTransData %>% group_by(StableId) %>% count %>% filter(n>1))

genePanel = read.csv('~/data/sv/rna_exp/gene_panel_ids.csv')
View(genePanel)

gpCoords = merge(genePanel,ensemblGeneData %>% select(GeneId,GeneName,Chromosome,GeneStart,GeneEnd),by=c('GeneId','GeneName'),all.x=T)
View(gpCoords)
write.csv(gpCoords %>% mutate(ChrStr=paste('Chr',Chromosome,sep=''),GeneStart=GeneStart-1000,GeneEnd=GeneEnd+1000) %>% select(ChrStr,GeneStart,GeneEnd),'~/logs/gene_panel.bed',row.names = F, quote = F)

sampleId = 'CPCT02020378T'

formFilename<-function(dir,sampleId,type)
{
   return (paste(dir,sampleId,'.',type,sep=''))
}

#####
## Gene RNA Expression
transData = read.csv(formFilename('~/logs/',sampleId,'transcript_data.csv'))
transData = read.csv(formFilename('~/logs/',sampleId,'gp_overs.transcript_data.csv'))
View(transData)
View(transData %>% filter(is.na(EmFitAllocation)))
View(transData %>% filter(GeneName=='HIST1H1C'))

geneData = read.csv(formFilename('~/logs/',sampleId,'gp_overs.gene_data.csv'))
geneData = read.csv(formFilename('~/logs/',sampleId,'all_gene_data.csv'))
View(geneData)

exonData = read.csv(formFilename('~/logs/',sampleId,'exon_data.csv'))
View(exonData)

View(rnaSampleData)
View(rnaSampleData %>% filter(SampleId %in% c('CPCT02020723T','CPCT02040216T','CPCT02080180T','CPCT02140007T','CPCT02110046T',
                       'CPCT02120030T','CPCT02210029T','DRUP01010047T','CPCT02070066T','CPCT02120162T','CPCT02040100T','CPCT02220059T','DRUP01090008T')))


## Alternative Splice Sites
altSJs = read.csv(formFilename('~/logs/',sampleId,'alt_splice_junc.csv'))
altSJs = read.csv(formFilename('~/logs/',sampleId,'gp.alt_splice_junc.csv'))
altSJs = read.csv(formFilename('~/logs/',sampleId,'gp_overs.alt_splice_junc.csv'))
View(altSJs)
nrow(altSJs)

# validation
tmpAltSJs = read.csv(formFilename('~/logs/',sampleId,'alt_splice_junc.csv'))
View(tmpAltSJs)

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000127616'))
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000127616'&(ExonEnd==11071850|ExonStart==11094800)))

View(altSJs %>% filter(FragCount>StartDepth|FragCount>EndDepth))

View(altSJs %>% group_by(Type) %>% count)
View(altSJs %>% group_by(GeneName) %>% count)
View(altSJs %>% group_by(Distance=5*round(NearestStartExon/5),StartContext) %>% count %>% spread(StartContext,n))


print(ggplot(altSJs %>% filter(StartContext=='INTRONIC') %>% group_by(Distance=50*round(NearestStartExon/50)) %>% count, aes(x=Distance, y=n))
      + geom_bar(stat='identity',colour='black')
      + xlim(0,10000)
      + labs(title='Distance to next exon for Intronic SJs'))

print(ggplot(altSJs %>% filter(StartContext=='EXONIC') %>% group_by(Distance=5*round(NearestStartExon/5)) %>% count, aes(x=Distance, y=n))
      + geom_bar(stat='identity',colour='black')
      + xlim(-500,0)
      + labs(title='Distance to next exon for Exonic SJs'))

View(ensemblTransExonData %>% filter(TransName=='ENST00000311189'))
View(ensemblGeneData %>% filter(GeneName=='SMAD4'|GeneName=='ELAC1'))

# comparison between 1 and 2 pass BAMs
altSJs1Pass = read.csv(formFilename('~/logs/',sampleId,'1pass.alt_splice_junc.csv'))
altSJs2Pass = read.csv(formFilename('~/logs/',sampleId,'2pass.alt_splice_junc.csv'))
nrow(altSJs1Pass)
nrow(altSJs2Pass)

altMerged = merge(altSJs1Pass,altSJs2Pass,by=c('GeneId','GeneName','Chromosome','Strand','SjStart','SjEnd'),all=T)
View(altMerged %>% mutate(FragDiff=FragCount.x-FragCount.y,FragDiffAbs=abs(FragDiff)) %>% filter(FragDiff!=0) %>%
        select(GeneId,GeneName,Chromosome,Strand,SjStart,SjEnd,FragDiff,FragDiffAbs,FragCount.x,FragCount.y,Type.x,Type.y,
               StartContext.x,StartContext.y,EndContext.x,EndContext.x,everything()))

## Read Data

readData = read.csv(formFilename('~/logs/',sampleId,'read_data.csv'))
View(readData)

View(readData %>% filter(GeneName=='SDHD'&GeneClass=='ALT'))
View(readData %>% filter(GeneName=='TOP2A'))


View(readData %>% filter(TransId=='ENST00000518797'&TransClass=='SPLICE_JUNCTION'&ExonEnd=='149782188'&ValidTrans==1))


matchedReads = readData %>% filter(ExonStart==18288659&GeneClass=='TRANS_SUPPORTING') %>% group_by(ReadId) %>% 
   summarise(Reads=n(),
             MinReadPos=min(PosStart),MaxReadPos=max(PosEnd),
             ExonStart=first(ExonStart),ExonEnd=first(ExonEnd)) %>% 
   mutate(ExonLength=ExonEnd-ExonStart+1,
          BaseOverlap=pmin(ExonEnd,MaxReadPos)-pmax(ExonStart,MinReadPos)+1)

print(sum(matchedReads$BaseOverlap)/nrow(matchedReads))
mean(matchedReads$BaseOverlap)
sum(matchedReads$BaseOverlap)

readsGrouped = readData %>% group_by(GeneId,GeneName,TransId,ReadId,ExonStart,ExonEnd) %>% 
   summarise(Reads=n(),InsertSize=abs(first(InsertSize)),
             Cigar1=first(Cigar),Cigar2=last(Cigar),
             HasNInsert=grepl('N',first(Cigar))|grepl('N',last(Cigar)),
             PosDiff1=abs(first(PosEnd)-last(PosStart))+1,PosDiff2=abs(last(PosEnd)-first(PosStart))+1) %>%
   mutate(InsertPosMatched=(InsertSize==PosDiff1|InsertSize==PosDiff2))

View(readsGrouped %>% filter(Reads==2))
View(readsGrouped %>% filter(Reads>2))


## Gene stats
print(ggplot(geneExpData %>% filter(TotalFragments>0) %>% group_by(ReadCounts=2**round(log(TotalFragments))) %>% count, aes(x=ReadCounts, y=n))
      + geom_line()
      # + scale_y_log10()
      + scale_x_log10())

geneExpData = geneExpData %>% mutate(ReadsPer1KB=TotalFragments/(GeneLength/1e3),
                                     TramsReadsPer1KB=SupportingTrans/(GeneLength/1e3))
View(geneExpData)

print(ggplot(geneExpData %>% filter(SupportingTrans>0&GeneLength>1e4&ReadsPer1KB<1e5) %>% filter(GeneName %in% genePanel$GeneName) %>%
                group_by(TramsReadsPer1KB=2**round(log(TramsReadsPer1KB))) %>% count, 
             aes(x=TramsReadsPer1KB, y=n))
      + geom_line()
      # + scale_y_log10()
)

print(ggplot(geneExpData %>% filter(SupportingTrans>0&GeneLength>1e4) %>% filter(GeneName %in% genePanel$GeneName) %>%
                group_by(SupportingTrans=2**round(log(SupportingTrans))) %>% count, 
             aes(x=SupportingTrans, y=n))
      + geom_line()
      # + scale_y_log10()
)

print(ggplot(geneExpData %>% filter(TotalFragments>0&GeneLength>1e4&ReadsPer1KB<1500) %>% group_by(ReadsPer1KB=2**round(log(ReadsPer1KB))) %>% count, 
             aes(x=ReadsPer1KB, y=n))
      + geom_line()
      # + scale_y_log10()
      )


## transcript fit validation
transData = read.csv(formFilename('~/logs/',sampleId,'gp_overs.transcript_data.csv'))
transData = read.csv(formFilename('~/logs/',sampleId,'transcript_data.csv'))
View(transData)
View(transData %>% filter(is.na(EmFitAllocation)))
View(transData %>% filter(GeneName=='HIST1H1C'))

transCatCounts = read.csv(formFilename('~/logs/',sampleId,'category_counts.csv'))
View(transCatCounts)

geneData = read.csv(formFilename('~/logs/',sampleId,'gp_overs.gene_data.csv'))
geneData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
View(geneData)



# Gene exon overlaps
geneExonOverlaps = read.csv('~/logs/gene_exon_overlaps.csv')
nrow(geneExonOverlaps)
View(geneExonOverlaps %>% group_by(GeneName1) %>% count)

View(geneExonOverlaps %>% filter(GeneName1 %in% genePanel$GeneName | GeneName2 %in% genePanel$GeneName) %>% group_by(GeneName1) %>% count %>% arrange(-n))

# Gene overlaps
geneOverlaps = read.csv('~/logs/gene_overlaps.csv')
nrow(geneOverlaps)
View(geneOverlaps)
View(geneOverlaps %>% group_by(GeneCount) %>% count)
View(geneExonOverlaps %>% group_by(GeneName1) %>% count)
View(genePanel)

View(geneOverlaps %>% filter(grepl('RBFOX1',GeneNames)))

View(geneOverlaps %>% filter(grepl('LRRFIP2',GeneNames)))
View(geneOverlaps %>% filter(grepl('MLH1',GeneNames)))

View(ensemblGeneData %>% filter(Chromosome==10))

# nrow(genePanel)
View(genePanel)
View(geneOverlaps)

View(ensemblGeneData %>% filter(grepl(GeneName,'ZNRF3;ZNRF3-IT1;ZNRF3-AS1')))

grepl('ZNRF3','ZNRF3;ZNRF3-IT1;ZNRF3-AS1')
tmp = str_split_fixed('ZNRF3;ZNRF3-IT1;ZNRF3-AS1',';',100)
tmp = strsplit('H3F3A=EIF2AK4;snoU13;H3F3AP1',';')
tmp = strsplit('H3F3A=EIF2AK4;snoU13;H3F3AP1',';')
View(ensemblGeneData %>% filter(GeneName %in% tmp[[1]]))
View(ensemblGeneData %>% filter(GeneName=='Y_RNA'))
View(tmp[[1]])

genesList = data.frame(matrix(ncol=2,nrow=0))
colnames(genesList) = c('GeneId','GeneName')

for(i in 1:nrow(genePanel))
{
   geneData = genePanel[i,]
   geneId = as.character(geneData$GeneId)
   geneName = as.character(geneData$GeneName)
   
   overlaps = geneOverlaps %>% filter(grepl(geneId,GeneNames))
   
   index = nrow(genesList)+1
   genesList[index,1] = geneId
   genesList[index,2] = geneName

   if(nrow(overlaps) > 0)
   {
      overlapData = overlaps[1,]
      overlappingGenes = as.character(overlapData$GeneNames)
      print(paste(geneName,'=',overlappingGenes,sep=''))
      
      tmp = strsplit(overlappingGenes,';')

      ensData = ensemblGeneData %>% filter(GeneId!=geneId&GeneId %in% tmp[[1]])
      print(paste(geneName,' matched ensembl records = ',nrow(ensData),sep=''))
      
      if(nrow(ensData) > 0)
      {
         for(j in 1:nrow(ensData))
         {
            ensGeneData = ensData[j,]
            newGeneId = as.character(ensGeneData$GeneId)
            newGeneName = as.character(ensGeneData$GeneName)
   
            index = nrow(genesList)+1
            genesList[index,1] = newGeneId
            genesList[index,2] = newGeneName
         }
      }
   }
}

View(genesList)
write.csv(genesList,'~/data/sv/rna_exp/gene_panel_plus_overlaps.csv',row.names = F, quote = F)




# 1 vs 2 pass BAMs from STAR

gene1PassData = read.csv(formFilename('~/logs/',sampleId,'gene_data_1pass.csv'))
View(gene1PassData)
gene2PassData = read.csv(formFilename('~/logs/',sampleId,'gene_data_2pass.csv'))
View(gene2PassData)

# total diffs
print(sum(gene1PassData$TotalFragments) - sum(gene2PassData$TotalFragments))
print(sum(gene1PassData$Alt) - sum(gene2PassData$Alt))
print(sum(gene1PassData$SupportingTrans) - sum(gene2PassData$SupportingTrans))
print(sum(gene1PassData$Chimeric) - sum(gene2PassData$Chimeric))
print(sum(gene1PassData$Unspliced) - sum(gene2PassData$Unspliced))
print(sum(gene1PassData$UnsplicedAlloc) - sum(gene2PassData$UnsplicedAlloc))

passMerge = merge(gene1PassData,gene2PassData,by=c('GeneId','GeneName','Chromosome'),all.x=T)
View(passMerge %>% filter(Alt.x!=Alt.y) %>% select(GeneId,GeneName,Alt.x,Alt.y,SupportingTrans.x,SupportingTrans.y,TotalFragments.x,TotalFragments.y,everything()))

View(passMerge %>% mutate(AltNew=Alt.y-Alt.x,SuppDiff=SupportingTrans.y-SupportingTrans.x) %>% 
        select(GeneId,GeneName,AltNew,Alt.x,Alt.y,SuppDiff,SupportingTrans.x,SupportingTrans.y,TotalFragments.x,TotalFragments.y,everything()))

readDataSlice1 = read.csv(formFilename('~/logs/',sampleId,'read_data.csv'))
View(readDataSlice1)

readDataSlice2 = read.csv(formFilename('~/logs/',sampleId,'read_data.csv'))
View(readDataSlice2)

readCompare = merge(readDataSlice2,readDataSlice1,by=c('GeneId','GeneName','ReadIndex','ReadId','TransId','Chromosome','PosStart','PosEnd'),all=T)
readCompare = merge(readDataSlice2,readDataSlice1,by=c('GeneId','GeneName','ReadIndex','ReadId','TransId','Chromosome'),all=T)
View(readCompare)
View(readCompare %>% filter(is.na(GeneClass.x)))
View(readCompare %>% filter(!is.na(GeneClass.x)&!is.na(GeneClass.y)&as.character(GeneClass.x)!=as.character(GeneClass.y)) %>%
        select(GeneId,PosStart.x,PosEnd.x,PosStart.y,PosEnd.y,GeneClass.x,GeneClass.y,Cigar.x,Cigar.y,everything()))

gene1PassData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
View(gene1PassData)
gene2PassData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
View(gene2PassData)




## Exon read depth study
exonData = read.csv(formFilename('~/logs/',sampleId,'exon_data.csv'))
View(exonData)
colnames(exonData)

singleTransGene = ensemblTransData %>% group_by(GeneId) %>% count %>% filter(n==1)
View(singleTransGene)

View(ensemblTransData)

exonSummary = exonData %>% mutate(Length=ExonEnd-ExonStart) %>% group_by(GeneId,GeneName,TransId) %>%
   summarise(Exons=n(),TotalExonBases=sum(Length),
             TransTotalCoverage=sum(TotalCoverage),
             TotalDepth=sum(AvgDepth*Length),
             FirstExonLength=(first(Length)),FirstExonAvgDepth=first(AvgDepth),FirstExonCovPerc=first(round(TotalCoverage/Length,1)),
             LastExonLength=(last(Length)),LastExonAvgDepth=last(AvgDepth),LastExonCovPerc=last(round(TotalCoverage/Length,1))) %>%
   mutate(TransAvgDepth=round(TotalDepth/TotalExonBases,1))

View(exonSummary %>% filter(GeneName=='SMAD4'))

View(exonSummary %>% filter(GeneId %in% singleTransGene$GeneId))

# depth by exon length
transExonCounts = exonData %>% group_by(TransId) %>% summarise(ExonCount=n())

exonData = merge(exonData,transExonCounts %>% select(TransId,ExonCount),by='TransId',all.x=T)

exonData = exonData %>% mutate(Length=ExonEnd-ExonStart)
rm(exonData)

# freq dist of exon lengths by first/last or not
print(ggplot(exonData %>% filter(ExonCount>2&Length<2000) %>% filter(GeneId %in% singleTransGene$GeneId) %>% 
                group_by(LengthBucket=round(Length,-2),IsFirst=(ExonRank==1)) %>% count, 
             aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~IsFirst))

print(ggplot(exonData %>% group_by(LengthBucket=round(ExonEnd-ExonStart,-2)) %>% summarise(Depth=median(AvgDepth)) %>% filter(Depth<500), aes(x=LengthBucket, y=Depth))
      + geom_line()
      + scale_x_log10()
      )

nrow(singleTransGene)

View(exonData %>% filter(ExonCount>2&Length<2000&AvgDepth<300) %>% filter(GeneId %in% singleTransGene$GeneId) %>%
        group_by(LengthBucket=round(Length,-1),IsFirst=(ExonRank==1)) %>% count)

print(ggplot(exonData %>% filter(ExonCount>2&Length<2000&AvgDepth<300) %>% filter(ExonCount==3) %>% filter(GeneId %in% singleTransGene$GeneId) %>% 
                group_by(LengthBucket=round(Length,-1),IsFirst=(ExonRank==1)) %>% 
                summarise(Depth=mean(AvgDepth)), aes(x=LengthBucket, y=Depth))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~IsFirst)
)

print(ggplot(exonData %>% filter(ExonCount>2&Length<2000&AvgDepth<300)  %>% filter(GeneId %in% singleTransGene$GeneId) %>% 
                group_by(DepthBucket=round(AvgDepth,-1),IsFirst=(ExonRank==1)) %>% count, aes(x=DepthBucket, y=n))
      + geom_line()
      # + scale_x_log10()
      + facet_wrap(~IsFirst)
)

View(exonData %>% filter(ExonCount>2&Length<2000) %>% filter(GeneId.x %in% singleTransGene$GeneId) %>% group_by(GeneId.x) %>% count)

exonLengthData = exonData %>% mutate(Unique=(SharedTrans==1),
                         Length=ExonEnd-ExonStart,
                         FirstOrLast=(ExonRank==1|ExonRank==ExonCount),
                         LengthBucket=round(Length,-2),
                         CovPerc=round(TotalCoverage/Length,2)) %>% 
        group_by(Unique,FirstOrLast,LengthBucket) %>% 
        summarise(Depth=median(AvgDepth),Coverage=median(CovPerc))


print(ggplot(exonLengthData %>% filter(Depth<500), aes(x=LengthBucket, y=Depth))
      + geom_line()
      # + scale_x_log10()
      + facet_wrap(~FirstOrLast)
      )


## Transcript prioritisation

transData = read.csv('~/logs/RNA_EXP_TRANS_DATA.csv')
View(transData)

tpData = transData %>% mutate(GeneName,TransId,ExonCount,
                              ExonsMatched,ExonPerc=round(ExonsMatched/ExonCount,1),
                              SJSupport=SpliceJuncSupported,SJSupportPerc=round(SpliceJuncSupported/(ExonCount-1),1),
                              UnqSJs=UniqueSpliceJunc,UnqSJPerc=ifelse(UniqueSpliceJunc>0,round(UniqueSpliceJuncSupported/UniqueSpliceJunc,1),0),
                              SJSupport=SpliceJuncSupported,SJSupportPerc=round(SpliceJuncSupported/(ExonCount-1),1),
                              UnqBases=UniqueBases,UnqBaseSupportPerc=ifelse(UniqueBases>0,round(UniqueBaseCoverage/UniqueBases,1),0))
                              
                              
View(tpData %>% select(GeneName,TransId,ExonCount,ExonsMatched,ExonPerc,SJSupport,SJSupportPerc,UnqSJs,UnqSJPerc,SJSupport,SJSupportPerc,
                       UnqBases,UnqBaseSupportPerc,UnqBaseDepth=UniqueBaseAvgDepth,ShortUnqFrags=ShortUniqueFragments,LongUnqeFrags=LongUniqueFragments,
                       SJFrags=SpliceJuncFragments,UnqSJFrags=UniqueSpliceJuncFragments,everything()))

colnames(transData)

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000174775'))

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000174775'&(ExonStart==533453|ExonEnd==533358)))



rsemResults = read.csv('~/data/sv/rna_exp/rsem_output.isoforms.results',sep='\t')
View(rsemResults)



#####
## Fragment Lengths

rnaSampleData = read.csv('~/data/sv/rna_exp/rna_exp_samples.csv')

load_frag_length_data<-function(sampleId)
{
   filename = formFilename('~/logs/',sampleId,'frag_length.csv')
   #sprintf('loading data for sample(%s) from file(%s)', sampleId,filename)
   print(paste('loading data for sample(',sampleId,') from file:',filename,sep = ''))
   
   fragLengths = read.csv(filename)
   fragLengths$SampleId = sampleId
   return (fragLengths %>% select(SampleId,FragmentLength,Count))
}

sampleFraglengths = load_frag_length_data('CPCT02010944T')
View(sampleFraglengths)
View(sampleFraglengths %>% filter(Count==0))

sampleFraglengths = data.frame()
for(sampleId in unique(rnaSampleData$SampleId))
{
   sampleFraglengths = rbind(sampleFraglengths,load_frag_length_data(sampleId))
}

View(sampleFraglengths)

write.csv(sampleFraglengths,'~/data/sv/rna_exp/sample_frag_lengths.csv',row.names = F,quote = F)

print(ggplot(sampleFraglengths, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,1000)
      + facet_wrap(~SampleId))

print(ggplot(sampleFraglengths %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count)) %>% filter(FragmentLength>=50&FragmentLength<100), aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000))

View(sampleFraglengths %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count)) %>% filter(FragmentLength>=50&FragmentLength<100))

# 151bp samples
samples151=c('CPCT02110046T','CPCT02120030T','CPCT02210029T','CPCT02070066T','CPCT02120162T','CPCT02220059T','DRUP01090008T')
sample151Fraglengths = data.frame()
for(sampleId in samples151)
{
   sample151Fraglengths = rbind(sample151Fraglengths,load_frag_length_data(sampleId))
}

View(sample151Fraglengths)
write.csv(sample151Fraglengths,'~/data/rna/exp/sample_151_frag_lengths.csv',row.names = F,quote = F)

print(ggplot(sample151Fraglengths, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,500)
      + facet_wrap(~SampleId))


fragLengths151Adj = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.033.csv'))
fragLengths151Adj$SampleId = 'CPCT02210029T_param_0.33'
fragLengths151Adj = fragLengths151Adj %>% select(SampleId,FragmentLength,Count)

fragLengths151_033 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.033.csv'))
fragLengths151_033$SampleId = 'CPCT02210029T_param_0.33'
fragLengths151_033 = fragLengths151_033 %>% select(SampleId,FragmentLength,Count)

fragLengths151_033_n = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.033_n.csv'))
fragLengths151_033_n$SampleId = 'CPCT02210029T_param_0.33_n'
fragLengths151_033_n = fragLengths151_033_n %>% select(SampleId,FragmentLength,Count)

fragLengths151_000 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.000.csv'))
fragLengths151_000$SampleId = 'CPCT02210029T_param_0.00'
fragLengths151_000 = fragLengths151_000 %>% select(SampleId,FragmentLength,Count)

View(fragLengthsCompare %>% group_by(SampleId) %>% summarise(Reads=sum(Count)))

fragLengthsCompare = rbind(fragLengths151_000,fragLengths151_033,fragLengths151_033_n,sample151Fraglengths %>% filter(SampleId=='CPCT02210029T'))
print(ggplot(fragLengthsCompare, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,300)
      + ylim(0,10000)
      + facet_wrap(~SampleId))

print(ggplot(fragLengthsCompare, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,100)
      + ylim(0,10000)
      + facet_wrap(~SampleId))


fragLengths151_p002 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length_p002.csv'))
fragLengths151_p002$SampleId = 'CPCT02210029T_param_002'
fragLengths151_p002 = fragLengths151_p002 %>% select(SampleId,FragmentLength,Count)

fragLengths151_p003 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.p003.csv'))
fragLengths151_p003$SampleId = 'CPCT02210029T_param_003'
fragLengths151_p003 = fragLengths151_p003 %>% select(SampleId,FragmentLength,Count)
View(fragLengths151_p003)

fragLengths151_004 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.p004.csv'))
fragLengths151_004$SampleId = 'CPCT02210029T_param_004'
fragLengths151_004 = fragLengths151_004 %>% select(SampleId,FragmentLength,Count)
View(fragLengths151_004)


fragLengths151_005 = read.csv(formFilename('~/logs/','CPCT02210029T','frag_length.p005.csv'))
fragLengths151_005$SampleId = 'CPCT02210029T_param_005'
fragLengths151_005 = fragLengths151_005 %>% select(SampleId,FragmentLength,Count)
View(fragLengths151_005)


View(fragLengthsCompare %>% group_by(SampleId) %>% summarise(Reads=sum(Count)))

print(ggplot(fragLengths151_005 %>% filter(FragmentLength<=300), aes(x=FragmentLength, y=Count))
      + geom_line())

print(ggplot(fragLengths151_p003, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,300)
      + ylim(0,10000))


## Retained Introns
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','retained_intron.csv'))
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','gp.retained_intron.csv'))
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','gp_overs.retained_intron.csv'))
View(retIntrons)
View(retIntrons %>% filter(FragCount>TotalDepth))
View(retIntrons %>% filter(GeneId=='ENSG00000244203'))

retIntronsTmp = read.csv(formFilename('~/logs/','CPCT02020378T','retained_intron.csv'))
View(retIntronsTmp)

View(ensemblTransExonData %>% filter(TransName=='ENST00000374690'))

View(geneOverlaps %>% group_by(GeneCount) %>% count)
sum(geneOverlaps$GeneCount)

matchingExonStarts = ensemblTransExonData %>% group_by(ExonStart,Strand) %>% summarise(GeneCount=n()) %>% group_by(ExonStart) %>%
   summarise(GeneCount=sum(GeneCount),StrandCount=n()) 
View(matchingExonStarts %>% filter(GeneCount>1))
View(matchingExonStarts %>% filter(GeneCount>1&StrandCount==2))

View(ensemblTransExonData %>% filter(ExonStart==230452))
nrow(ensemblTransExonData %>% filter(ExonStart==230452))

altSJs = read.csv(formFilename('~/logs/','CPCT02020378T','gp_overs.alt_splice_junc.csv'))
altSJs = altSJs %>% mutate(StartSeq=stri_sub(StartBases,3,4),EndSeq=stri_sub(EndBases,9,10))
View(altSJs %>% filter(Strand==1))
View(altSJs %>% filter(Strand==1) %>% group_by(Seq=paste(StartSeq,EndSeq,sep='_'),Type) %>% count %>% spread(Type,n,fill=0))
View(altSJs %>% filter(FragCount>StartDepth|FragCount>EndDepth))

View(altSJs %>% filter(StartSeq=='CT'&EndSeq=='AC'))


expRates = read.csv(formFilename('~/logs/','CPCT02020378T','gp.exp_rates.csv'))
View(expRates)
View(expRates %>% group_by(GeneSetId) %>% count)


## Summary of 151b reads samples

# 1x extreme 151b samples
sampleId = 'CPCT02210029T'
geneData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
geneData$SampleId = sampleId

topGenes = head(geneData %>% arrange(-TotalFragments),10)

geneDataSum = geneData %>% filter(!(GeneName %in% topGenes$GeneName)) %>% group_by(SampleId) %>% summarise(GeneName='OTHER',TotalFragments=sum(TotalFragments))
geneDataSum = rbind(geneDataSum,geneData %>% filter(GeneName %in% topGenes$GeneName) %>% select(SampleId,GeneName,TotalFragments))
View(geneDataSum)
allGeneData = geneDataSum

# 1x normal 151b samples
sampleId = 'CPCT02110046T'
geneData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
geneData$SampleId = sampleId
geneDataSum = geneData %>% filter(!(GeneName %in% topGenes$GeneName)) %>% group_by(SampleId) %>% summarise(GeneName='OTHER',TotalFragments=sum(TotalFragments))
geneDataSum = rbind(geneDataSum,geneData %>% filter(GeneName %in% topGenes$GeneName) %>% select(SampleId,GeneName,TotalFragments))
View(geneDataSum)
allGeneData = rbind(allGeneData,geneDataSum)

# 2x 76b samples
sampleId = 'CPCT02020378T'
geneData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
geneData$SampleId = sampleId
geneDataSum = geneData %>% filter(!(GeneName %in% topGenes$GeneName)) %>% group_by(SampleId) %>% summarise(GeneName='OTHER',TotalFragments=sum(TotalFragments))
geneDataSum = rbind(geneDataSum,geneData %>% filter(GeneName %in% topGenes$GeneName) %>% select(SampleId,GeneName,TotalFragments))
View(geneDataSum)
allGeneData = rbind(allGeneData,geneDataSum)

sampleId = 'CPCT02010944T'
geneData = read.csv(formFilename('~/logs/',sampleId,'gene_data.csv'))
geneData$SampleId = sampleId
geneDataSum = geneData %>% filter(!(GeneName %in% topGenes$GeneName)) %>% group_by(SampleId) %>% summarise(GeneName='OTHER',TotalFragments=sum(TotalFragments))
geneDataSum = rbind(geneDataSum,geneData %>% filter(GeneName %in% topGenes$GeneName) %>% select(SampleId,GeneName,TotalFragments))
View(geneDataSum)
allGeneData = rbind(allGeneData,geneDataSum)

sampleTotals = allGeneData %>% group_by(SampleId) %>% summarise(SampleFragCount=sum(TotalFragments))
allGeneData = merge(allGeneData,sampleTotals,by='SampleId',all.x=T)
allGeneData = allGeneData %>% mutate(GenePerc=round(TotalFragments/SampleFragCount,4))
View(allGeneData)

colors = c("wheat2", "hotpink", "darkorange", "seagreen3", "gray", "thistle2", "steelblue2", "darkgreen", "indianred", "honeydew2",
                          "turquoise3", "lightpink2", "goldenrod2", "cornsilk3", "yellowgreen", "wheat2", "violetred2", "ivory3", "coral1", "springgreen2")

print(ggplot(allGeneData, aes(x=SampleId, y=GenePerc, fill=GeneName))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='', y='% of Sample SFragments', title='% of Total Fragment Counts for Top Genes vs All Other Genes')
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + scale_fill_manual(values = colors))

# fragment length for 2x 151b samples
sampleId = 'CPCT02210029T'
fragLengths = read.csv(formFilename('~/logs/',sampleId,'frag_length.p004.csv'))
fragLengths$SampleId = sampleId
fragLengths151_004 = fragLengths151_004 %>% select(SampleId,FragmentLength,Count)
View(fragLengths)

print(ggplot(fragLengths151_004, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,300)
      + ylim(0,10000))


## GC-Bias
# 2x 151b samples
sampleId = 'CPCT02110046T'
gcData = read.csv(formFilename('~/logs/',sampleId,'read_gc_ratios.csv'))
gcData$SampleId = sampleId
allGcData = gcData %>% filter(GeneName=='ALL') %>% select(SampleId,GeneName,GcRatio,Count)
View(allGcData)

sampleId = 'CPCT02210029T'
gcData = read.csv(formFilename('~/logs/',sampleId,'read_gc_ratios.csv'))
gcData$SampleId = sampleId
allGcData = rbind(allGcData,gcData %>% filter(GeneName=='ALL') %>% select(SampleId,GeneName,GcRatio,Count))

# 2x 76b samples
sampleId = 'CPCT02020378T'
gcData = read.csv(formFilename('~/logs/',sampleId,'read_gc_ratios.csv'))
gcData$SampleId = sampleId
allGcData = rbind(allGcData,gcData %>% filter(GeneName=='ALL') %>% select(SampleId,GeneName,GcRatio,Count))

sampleId = 'CPCT02010944T'
gcData = read.csv(formFilename('~/logs/',sampleId,'read_gc_ratios.csv'))
gcData$SampleId = sampleId
allGcData = rbind(allGcData,gcData %>% filter(GeneName=='ALL') %>% select(SampleId,GeneName,GcRatio,Count))

print(ggplot(allGcData, aes(x=GcRatio, y=Count))
      + geom_line()
      + facet_wrap(~SampleId))



View(fragLengthsCompare %>% group_by(SampleId) %>% summarise(Reads=sum(Count)))

print(ggplot(fragLengths151_004, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,300)
      + ylim(0,10000))




## Fragment Length Read Data
fragLengthReads76 = read.csv('~/logs/CPCT02020378T.frag_length_reads.csv')
fragLengthReads = read.csv('~/logs/CPCT02210029T.frag_length_reads.csv')
fragLengthReadsPrev = fragLengthReads
fragLengthReads = read.csv('~/logs/CPCT02210029T.frag_length_reads.p003.csv')
View(fragLengthReads)
View(fragLengthReads %>% mutate(InsertSize=abs(InsertSize)) %>% filter(InsertSize>0&InsertSize<76))
View(fragLengthReads %>% filter(abs(InsertSize)<=76,ReadLength==76) %>% group_by(Cigar) %>% count)
View(fragLengthReads %>% mutate(InsertSize=abs(InsertSize)) %>% filter(InsertSize==100))

tmpFLReads = fragLengthReads %>% mutate(InsertSize=abs(InsertSize)) %>% filter(InsertSize<180)
View(tmpFLReads %>% filter(InsertSize>0) %>% group_by(InsertSize) %>% count)

print(ggplot(tmpFLReads %>% filter(InsertSize>0) %>% group_by(InsertSize) %>% count, aes(x=InsertSize, y=n))
      + geom_line())


rnaSampleReads = read.csv('~/logs/rna_sample_reads.csv')
View(rnaSampleReads)
rnaSampleReads = rnaSampleReads %>% mutate(ReadLength=pmax(stri_length(ReadBases),76))
View(rnaSampleReads %>% group_by(ReadLength) %>% count)

sampleCancerTypes = read.csv('~/data/hpc_sample_cancer_types.csv')
rnaCTs = read.csv('~/data/rna/rna_sample_cancer_types.csv')
rnaCTs = rbind(rnaCTs,sampleCancerTypes %>% select(SampleId,CancerType))

rnaSampleReads = merge(rnaSampleReads,rnaCTs,by='SampleId',all.x=T)
View(rnaSampleReads)
View(rnaSampleReads %>% filter(is.na(CancerType)))
View(rnaSampleReads %>% filter(is.na(CancerType)) %>% select(SampleId))
View(rnaSampleReads %>% group_by(CancerType,ReadLength) %>% count %>% spread(ReadLength,n,fill=0))
View(rnaSampleReads %>% group_by(CancerType) %>% count)

rnaSampleData = rnaSampleReads %>% mutate(FastqDir=paste(Dir,SampleId,sep='_')) %>% select(-ReadBases,-Dir)
View(rnaSampleData)
write.csv(rnaSampleData,'~/logs/rna_sample_data.csv',row.names = F,quote = F)
write.csv(rnaSampleData %>% filter(CancerType=='Urinary tract') %>% select(FastqDir),'~/logs/rna_urinary_tract_fastq_dir.txt',row.names = F,quote = F)


tmpFLReads76 = fragLengthReads76 %>% mutate(InsertSize=abs(InsertSize)) %>% filter(InsertSize<160)

print(ggplot(tmpFLReads76 %>% filter(InsertSize>0&InsertSize<100) %>% group_by(InsertSize) %>% count, aes(x=InsertSize, y=n))
      + geom_line())

View(tmpFLReads76 %>% filter(InsertSize>=70&InsertSize<=77) %>% group_by(InsertSize,Cigar) %>% count %>% spread(InsertSize,n,fill=0))



## GC Bias

gcRatios = read.csv('~/data/rna/exp/gcbias_ratios.csv')
View(gcRatios)
colnames(gcRatios) = c('Chromosome','Region','GcRatio')
gcRatioFreq = gcRatios %>% filter(GcRatio>0) %>% group_by(GcRatio) %>% count
View(gcRatioFreq)

print(ggplot(gcRatioFreq, aes(x=GcRatio, y=n))
      + geom_line())

geneData = read.csv(formFilename('~/logs/','CPCT02210029T','gene_data.csv'))
View(geneData)

fastqDirs = read.csv('~/logs/fastq_dir.txt')
View(fastqDirs)
write.csv(fastqDirs %>% select(SampleId,FastqDir),'~/logs/rna_fastq_dir.csv',row.names = F,quote = F)
fastqDirs = fastqDirs %>% mutate(FastqDir=paste(Dir,SampleId,sep='_')) %>% select(SampleId,FastqDir)

fastqDirs = merge(fastqDirs,sampleCancerTypes %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(fastqDirs %>% filter(CancerType=='Urinary Tract'))
View(fastqDirs %>% group_by(CancerType) %>% count)

View(fastqDirs %>% select(SampleId))

View(rnaSampleData)


readGcRatios = read.csv(formFilename('~/logs/','CPCT02210029T','read_gc_ratios.csv'))
readGcRatios = read.csv(formFilename('~/logs/','CPCT02210029T','dups.read_gc_ratios.csv'))
View(readGcRatios)
View(readGcRatios %>% filter(GcRatio>=0.9))
View(readGcRatios %>% filter(GeneName=='RPS29'))
View(readGcRatios %>% group_by(GeneName) %>% summarise(Total=sum(Count)))


filteredReads = readGcRatios %>% filter(GeneName!='ALL'&GeneName!='RN7SL2'&GeneName!='RN7SL1'&GeneName!='RPS29'&GcRatio>0.57&GcRatio<0.63) %>% 
   group_by(GeneName) %>% summarise(Total=sum(Count))

highReadGenes = head(geneData %>% arrange(-Duplicates),100)
View(highReadGenes)
View(readGcRatios %>% filter(!(GeneName %in% highReadGenes$GeneName)))
filteredReads = readGcRatios %>% filter(!(GeneName %in% highReadGenes$GeneName)&GeneName!='ALL') %>% group_by(GcRatio) %>% summarise(Total=sum(Count))
View(filteredReads)
print(ggplot(filteredReads, aes(x=GcRatio, y=Total))
      + geom_line())


View(highReadGenes)

print(ggplot(readGcRatios %>% filter(GeneName=='ALL'), aes(x=GcRatio, y=Count))
      + geom_line())

print(ggplot(readGcRatios %>% filter(GeneName=='RPS29'), aes(x=GcRatio, y=Count))
      + geom_line())

print(ggplot(readGcRatios %>% filter(GeneId %in% highReadGenes$GeneId), aes(x=GcRatio, y=Count))
      + geom_line()
      + facet_wrap(~GeneName))

tmpGcRatios = read.csv(formFilename('~/logs/',sampleId,'gp.read_gc_ratios.csv'))
View(tmpGcRatios)

print(ggplot(tmpGcRatios %>% filter(GeneName=='ALL'), aes(x=GcRatio, y=Count))
      + geom_line())



View(fragLengths)
nrow(fragLengths)
View(fragLengths %>% filter(is.na(Count)))
sum(fragLengths %>% filter(!is.na(Count)))
nrow(fragLengths %>% filter(!is.na(Count)) %>% filter(FragmentLength<400))

fragLengthSummary = fragLengths %>% filter(!is.na(Count)) %>% group_by(FragLengthBucket=round(FragmentLength,-2)) %>% summarise(Total=sum(Count))
View(fragLengthSummary)

print(ggplot(fragLengths, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,1000)
)

View(fragLengths %>% group_by(FragLengthBucket=50*round(FragmentLength/50)) %>% summarise(Total=sum(Count)))

print(ggplot(fragLengths %>% group_by(FragLengthBucket=50*round(FragmentLength/50)) %>% summarise(Total=sum(Count)), aes(x=FragLengthBucket, y=Total))
      + geom_bar(stat='identity',colour='black')
      + xlim(0,1000))



print(ggplot(fragLenByGene %>% group_by(FragmentLength,IsRPS29=(GeneName %in% c('RPS29'))) %>% 
                summarise(FragCount=sum(Count)), aes(x=FragmentLength, y=FragCount))
      + geom_line()
      # + xlim(140,160)
      + xlim(0,500)
      + facet_wrap(~IsRPS29))

print(ggplot(fragLenByGene %>% filter(GeneName=='RPS29') %>% group_by(FragmentLength) %>% 
                summarise(FragCount=sum(Count)), aes(x=FragmentLength, y=FragCount))
      + geom_line()
      # + xlim(140,160)
      + xlim(0,500))

print(ggplot(fragLengths, aes(x=FragmentLength, y=Count))
      + geom_line()
      # + scale_y_log10()
      # + scale_x_log10()
      + xlim(0,1000)
)


# clean-up
rm(fragLenByGene)
rm(fragLenByGene75b)
rm(fragLenByGene150)
rm(fragLengths)

fragLenByGene = read.csv('~/logs/RNA_EXP_FRAG_LENGTHS.csv')

# 150b reads
fragLenByGene = read.csv('~/logs/CPCT02120030T_FRAG_LENGTHS_BY_GENE.csv')
fragLenByGene = read.csv('~/logs/CPCT02110046T_FRAG_LENGTHS_BY_GENE.csv')
fragLenByGene = read.csv('~/logs/CPCT02210029T_FRAG_LENGTHS_BY_GENE.csv')
nrow(fragLenByGene)

fragLengthsAll = fragLenByGene %>% filter(GeneName!='RPS29') %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count))
View(fragLengthsAll)
View(fragLenByGene)

print(ggplot(fragLengthsAll, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000))

print(ggplot(fragLengthsAll, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(100,400))

fragLenByGene75b = read.csv('~/logs/CPCT02020378T_FRAG_LENGTHS_BY_GENE.csv')
fragLenByGene75b = read.csv('~/logs/CPCT02020834T_FRAG_LENGTHS_BY_GENE.csv')
fragLenByGene75b = read.csv('~/logs/CPCT02010963T_FRAG_LENGTHS_BY_GENE.csv')

fragLengthsAll75b = fragLenByGene75b %>% filter(GeneName!='RPS29') %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count))
View(fragLengthsAll75b)
View(fragLenByGene75b)

print(ggplot(fragLengthsAll75b, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000))

print(ggplot(fragLengthsAll75b, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(50,500))


fragLenByGene150 = read.csv('~/logs/CPCT02110046T_FRAG_LENGTHS_BY_GENE.csv')
fragLengthsAll150 = fragLenByGene150 %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count))
View(fragLengthsAll150)

View(fragLenByGene150 %>% group_by(GeneName) %>% summarise(LengthBuckets=n(),
                                                           Fragments=sum(Count)))

print(ggplot(fragLengthsAll150, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000)
)


print(ggplot(fragLenByGene150 %>% group_by(FragmentLength,IsRPS29=(GeneName %in% c('RPS29','RP5-857K21.4'))) %>% 
                summarise(FragCount=sum(Count)), aes(x=FragmentLength, y=FragCount))
      + geom_line()
      # + xlim(140,160)
      + xlim(0,500)
      + facet_wrap(~IsRPS29))


print(ggplot(fragLengthsAll150, aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000))

print(ggplot(fragLenByGene150 %>% filter(GeneName=='TRAPPC9'), aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,1000))


View(ensemblTransExonData %>% filter(GeneId=='ENSG00000131747'))

# gene investigations
View(fragLenByGene150 %>% filter(Count>100))


flSample = read.csv('~/logs/CPCT02020894T_FRAG_LENGTHS.csv')
flSample150bp = read.csv('~/logs/CPCT02110046T_FRAG_LENGTHS.csv')
flSample150bp2 = read.csv('~/logs/CPCT02210029T_FRAG_LENGTHS.csv')
flSample75bp = read.csv('~/logs/CPCT02020378T_FRAG_LENGTHS.csv')
View(flSample75bp)
# CPCT02110046T_FRAG_LENGTHS.csv
View(flSample)
View(flSample150bp)

print(ggplot(flSample, aes(x=FragmentLength, y=Count))
      + geom_line()
      # + scale_y_log10()
      # + scale_x_log10()
      + xlim(0,1000)
)

print(ggplot(flSample75bp, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,200)
)

print(ggplot(flSample150bp, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,500))

print(ggplot(flSample150bp2, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,500))

# by chromomsome
flSampleChr150bp = read.csv('~/logs/CPCT02110046T_FRAG_LENGTHS_CHR20.csv')
flSampleChr150bp$Chromosome = 20
flSampleChr = flSampleChr150bp
flSampleChr = rbind(flSampleChr,flSampleChr150bp)

View(flSampleChr)

print(ggplot(flSampleChr, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,250)
      + facet_wrap(~Chromosome))

View(ensemblGeneData %>% filter(Chromosome==1&GeneEnd-GeneStart>1e4&GeneEnd-GeneStart<1e6))

write.csv(ensemblGeneData %>% filter(Chromosome==17&GeneEnd-GeneStart>1e4&GeneEnd-GeneStart<1e6) %>% select(GeneId,GeneName),
          '~/logs/genes_chr17.csv', row.names = F, quote = F)

write.csv(ensemblGeneData %>% filter(Chromosome==1&GeneEnd-GeneStart>1e4&GeneEnd-GeneStart<1e6) %>% select(GeneId,GeneName),
          '~/logs/genes_chr1.csv', row.names = F, quote = F)


genePanel = read.csv('~/data/full_gene_panel.csv')
View(genePanel)
nrow(genePanel)

oncoGenes = read.csv('~/data/oncogenes.csv')
View(oncoGenes)
nrow(oncoGenes)

genePanel = rbind(genePanel,oncoGenes)
genePanel = merge(ensemblGeneData %>% select(GeneId,GeneName),genePanel,by='GeneName',all.y = T)
write.csv(genePanel %>% select(GeneId,GeneName),'~/data/sv/rna/gene_panel_ids.csv',row.names = F, quote = F)




geneTransData = ensemblTransExonData %>% group_by(GeneId,Trans) %>% summarise(ExonCount=n()) %>% group_by(GeneId) %>% summarise(TransCount=n())
View(geneTransData)

transCodingLengths = ensemblTransExonData %>% group_by(GeneId,TransId,Trans,IsCanonical=CanonicalTranscriptId==TransId) %>% summarise(Exons=n(),
                                                                                    ExonTotalLength=sum(ExonEnd-ExonStart))
nrow(transCodingLengths)
View(transCodingLengths %>% filter(IsCanonical))
View(transCodingLengths %>% filter(Trans=='ENST00000398905'))

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000164362') %>% arrange(ExonStart))

View(ensemblTransExonData %>% filter(Trans=='ENST00000576024'))

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000184012') %>% filter((ExonStart<=42839661&ExonEnd>=42839661)|(ExonStart<=42845345&ExonEnd>=42845345)))
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000157554'))





