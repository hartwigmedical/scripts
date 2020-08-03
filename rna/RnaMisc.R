library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
library(stringr)
library(hexbin)


## Ensembl ref data
ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)

ensemblTransExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000106546'))
View(ensemblTransExonData %>% filter(TransName=='ENST00000603077'))

ensemblTransData = read.csv('~/data/sv/ensembl_trans_data.csv')
View(ensemblTransData)
View(ensemblTransData %>% group_by(StableId) %>% count %>% filter(n>1))

genePanel = read.csv('~/data/sv/rna_exp/gene_panel_ids.csv')
View(genePanel)

sampleId = 'CPCT02020378T'

formFilename<-function(dir,sampleId,type)
{
   return (paste(dir,sampleId,'.isf.',type,sep=''))
}


####################
## Isofox data types
sampleId = 'CPCT02020378T'
sampleId = 'CPCT02210029T'
sourceDir='~/logs/'
sourceDir='~/data/rna/runs/'
sourceDir='~/data/rna/gcp_data/'

transData = read.csv(formFilename(sourceDir,sampleId,'transcript_data.csv'))
View(transData)

geneSetData = read.csv(formFilename(sourceDir,sampleId,'gene_collection_data.csv'))
View(geneSetData)

geneData = read.csv(formFilename(sourceDir,sampleId,'gene_data.csv'))
View(geneData)

summaryData = read.csv(formFilename(sourceDir,sampleId,'summary.csv'))

altSJ = read.csv(formFilename(sourceDir,sampleId,'alt_splice_junc.csv'))
View(altSJ)

gcRatioData = read.csv(formFilename(sourceDir,sampleId,'gc_ratio_data.csv'))
View(gcRatioData)

readData = read.csv(formFilename('~/logs/',sampleId,'read_data.csv'))
View(readData)

exonData = read.csv(formFilename('~/logs/',sampleId,'exon_data.csv'))
View(exonData)

enrichedGeneIds = c('ENSG00000265150','ENSG00000258486','ENSG00000202198','ENSG00000266037','ENSG00000263740','ENSG00000265735')


# Cohort Stats
starStats = read.csv('~/data/rna/star_stats.csv')
View(starStats)
colnames(starStats)
View(starStats %>% group_by(SampleId) %>% count)
View(rnaSampleData %>% group_by(SampleId) %>% count)

View(rnaSampleData %>% filter(SampleId %in% starStats$SampleId))
write.csv(rnaSampleData %>% filter(SampleId %in% starStats$SampleId) %>% select(SampleId,ReadLength),'~/data/rna/samples/batch_isofox_completed_samples.csv',row.names=F,quote=F)


# PGBD5 cohort
pgbd5Tpm = read.csv('~/data/rna/samples/isfox_pgbd5_all.transcript_consolidated.csv')
pgbd5Tpm = pgbd5Tpm %>% mutate(LogTPM=ifelse(TPM>0,log(TPM),0))
View(pgbd5Tpm)
View(pgbd5Tpm %>% group_by(TransName) %>% summarise(Samples=n(),
                                                    EffectiveLength=first(EffectiveLength),
                                                    MinTpm=min(TPM),
                                                    AvgTpm=mean(TPM),
                                                    MedianTpm=median(TPM),
                                                    MaxTpm=max(TPM)))

print(ggplot(pgbd5Tpm %>% filter(TPM>0.001), aes(x=TransName,y=log(TPM,10)))
      + geom_violin(scale="count",fill="#6baed6"))

pgbd5Tpm = merge(pgbd5Tpm,rnaSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(pgbd5Tpm %>% filter(CancerType=='Urinary tract'))


rnaFastqFiles = read.csv('~/data/rna/rna_fastq_files.csv')
View(rnaFastqFiles)
View(rnaFastqFiles %>% group_by(FastqDir) %>% count)



# Depth Comparison
oldAltSJ = read.csv(formFilename(sourceDir,sampleId,'alt_splice_junc.csv'))
View(oldAltSJ)

newAltSJ = read.csv(formFilename(sourceDir,sampleId,'alt_splice_junc.csv'))
View(newAltSJ)

altMerged = merge(oldAltSJ %>% select(GeneId,SjStart,SjEnd,DepthStart,DepthEnd),newAltSJ %>% select(GeneId,SjStart,SjEnd,DepthStart,DepthEnd),
                  by=c('GeneId','SjStart','SjEnd'),all=T)

View(altMerged)

## Summary data
utSummaryData = read.csv('~/data/rna/logs/isofox_summary.csv')
View(utSummaryData)
colnames(utSummaryData)

utSummaryData$CancerType = 'Urinary tract'

print(ggplot(utSummaryData, aes(x=SampleId))
      + geom_point(aes(y=FragLength5th, colour='FragLength5th')))

print(ggplot(utSummaryData, aes(x=SampleId, y=MedianGCRatio))
      + geom_point())

print(ggplot(utSummaryData, aes(x=CancerType,y=MedianGCRatio))
      + geom_violin(scale="count",fill="#6baed6"))

print(ggplot(utSummaryData, aes(x=CancerType,y=EnrichedGenePercent))
      + geom_violin(scale="count",fill="#6baed6"))

utSummaryFragLengths = utSummaryData %>% filter(FragLength95th<1000) %>% select(SampleId,FragLength5th,FragLength50th,FragLength95th) %>%
   gather('Percentile','Count', 2:4)
View(utSummaryFragLengths)

print(ggplot(utSummaryFragLengths, aes(x=Percentile,y=Count))
      + geom_violin(scale="count",fill="#6baed6"))

print(ggplot(utSummaryData, aes(x=CancerType,y=FragLength50th))
      + geom_violin(scale="count",fill="#6baed6"))

print(ggplot(utSummaryData, aes(x=CancerType,y=FragLength5th))
      + geom_violin(scale="count",fill="#6baed6"))

print(ggplot(utSummaryData %>% filter(FragLength95th<1000), aes(x=CancerType,y=FragLength95th))
      + geom_violin(scale="count",fill="#6baed6"))


qcStats = cohortCombinedStats
qcStats = read.csv('~/hmf/analyses/RNA/star_isofox_summary.csv')

# some samples have tiny amounts of total fragments
ggplot(qcStats, aes(TotalFragments)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
View(qcStats %>% filter(TotalFragments<2e6))

#

ggplot(qcStats, aes(FragLength5th)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(FragLength50th)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(FragLength95th)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(EnrichedGenePercent)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(MedianGCRatio)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(TooShortReads)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)

#151 base have fewer too many loci but similar multiple loci
ggplot(qcStats, aes(TooManyLociReads)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)
ggplot(qcStats, aes(MultipleLoci)) + stat_ecdf(geom = "step") +  facet_wrap(~ReadLength)

# too short reads and multiple loci reads not strongly associated
ggplot(qcStats, aes(TooShortReads,TooManyLociReads)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(MultipleLoci,TooManyLociReads)) + geom_point()+  facet_wrap(~ReadLength)

#High enriched gene percent associated with high GC, long 95th%, short 50th%, short 5th% and too shortReads
ggplot(qcStats, aes(EnrichedGenePercent,MedianGCRatio)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(EnrichedGenePercent,FragLength95th)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(EnrichedGenePercent,FragLength50th)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(EnrichedGenePercent,FragLength5th)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(EnrichedGenePercent,TooShortReads)) + geom_point()+  facet_wrap(~ReadLength)
ggplot(qcStats, aes(FragLength50th,TooShortReads)) + geom_point()+  facet_wrap(~ReadLength)


# Samples with few fragments highly associated with multiple loci
ggplot(qcStats, aes(TotalFragments,MultipleLoci)) + geom_point()+  facet_wrap(~ReadLength)


## Validation of test (151p read) samples

testSamples = c('CPCT02010409T','CPCT02120066T','CPCT02010419T','CPCT02210029T','CPCT02010414T')
testSamples2 = c('CPCT02010409T','CPCT02120066T','CPCT02010419T','CPCT02010414T')

transData = read.csv('~/data/rna/runs/CPCT02020378T.isf.transcript_data.csv')
transData = read.csv('~/data/rna/gcp_data/CPCT02210029T.isf.transcript_data.csv')
transData = transData %>% mutate(GcFitDiff=abs(FitAllocation-PreGcFit),GcFitDiffPerc=ifelse(PreGcFit>0,round(GcFitDiff/PreGcFit,3),0))
View(transData)


transDataCombined = data.frame()
for(sampleId in testSamples)
{
   transData = read.csv(formFilename('~/data/rna/gcp_data/',sampleId,'transcript_data.csv'))
   transData$SampleId=sampleId

   transData = transData %>% filter(FitAllocation>10|PreGcFit>10)
   transData = transData %>% mutate(GcFitDiff=abs(FitAllocation-PreGcFit),GcFitDiffPerc=ifelse(PreGcFit>0,round(GcFitDiff/PreGcFit,3),0))
   transData = transData %>% select(GeneId,GeneName,TransName,Canonical,ExonCount,PreGcFit,FitAllocation,GcFitDiff,GcFitDiffPerc,EffectiveLength,TPM,everything())
   transDataCombined = rbind(transDataCombined,transData)
}

View(transDataCombined)
write.csv(transDataCombined,'~/logs/trans_data_5_151b_samples.csv',row.names = F,quote = F)

geneDataCombined = data.frame()
for(sampleId in testSamples)
{
   geneData = read.csv(formFilename('~/data/rna/gcp_data/',sampleId,'gene_data.csv'))
   geneData$SampleId=sampleId
   
   geneData = geneData %>% filter(SplicedFragments>1e5|UnsplicedFragments>1e5)
   geneDataCombined = rbind(geneDataCombined,geneData)
}

View(geneDataCombined)

geneSetDataCombined = data.frame()
for(sampleId in testSamples)
{
   geneSetData = read.csv(formFilename('~/data/rna/gcp_data/',sampleId,'gene_collection_data.csv'))
   geneSetData$SampleId=sampleId
   
   # geneSetData = geneSetData %>% filter(TotalFragments>1e5)
   geneSetDataCombined = rbind(geneSetDataCombined,geneSetData)
}

View(geneSetDataCombined)
View(geneSetDataCombined %>% filter(Duplicates>TotalFragments))


#####
## Enriched Genes

## confirmation that no other genes' exons overlap these at all
enrChroms = c(14,14,14,14,3,6,9)
enrExonStarts = c(50053297, 50053298, 50329271, 50320336, 15780022, 52860418, 9442060) 
enrExonEnds = c(50053596, 50053594, 50329567, 50320632, 15780315, 52860748, 9442347)

for(i in 1:7)
{
   print(sprintf("chr(%s) exons(%d -> %d)", enrChroms[i], enrExonStarts[i], enrExonEnds[i]))
   
   chrGenes = ensemblGeneData %>% filter(Chromosome==enrChroms[i])
   
   overlaps = ensemblTransExonData %>% 
      filter(GeneId %in% chrGenes$GeneId) %>%
      filter(!(GeneId %in% enrichedGeneIds)) %>% 
      filter(!(ExonStart>enrExonEnds[i]|ExonEnd<enrExonStarts[i]))
   
   if(nrow(overlaps))
   {
      for(j in 1:nrow(overlaps))
      {
         transData = overlaps[j,]
         print(sprintf("trans(%s:%s) overlaps", transData$GeneId, transData$TransName))
      }
   }
}

View(geneData %>% filter(GeneSetId %in% c('14_462','14_471','14_472')))


######
## GC Bias and Ratios

## GC Expected Counts
transGcExpRatios = read.csv('~/data/rna/read_100_exp_gc_ratios.csv')
View(head(transGcExpRatios,100))
View(transGcExpRatios)

gcRatioData = read.csv(formFilename(sourceDir,sampleId,'gc_ratio_data.csv'))
gcRatioData = read.csv('~/logs/CPCT02020378T.isf.gc_ratio_data.csv')


gcRatioDataCombined = data.frame()
for(sampleId in testSamples)
{
   gcRatioData = read.csv(formFilename('~/data/rna/gcp_data/',sampleId,'gc_ratio_data.csv'))
   gcRatioData$SampleId=sampleId
   gcRatioData = gcRatioData %>% select(SampleId,everything())
   gcRatioDataCombined = rbind(gcRatioDataCombined,gcRatioData)
}

View(gcRatioDataCombined)
write.csv(gcRatioDataCombined,'~/logs/gc_ratio_data_5_151b_samples.csv',row.names = F,quote = F)

# Plot GC ratios
gcRatioData2 = gcRatioDataCombined %>% gather('Ratio','Count',3:ncol(gcRatioData))
gcRatioData2 = gcRatioData2 %>% mutate(RatioDec=stri_replace_all_fixed(Ratio,'Gcr_',''))
gcRatioData2 = gcRatioData2 %>% mutate(Count=as.numeric(as.character(Count)))
gcRatioData2$Category = gcRatioData2$GeneName
gcRatioData2 = within(gcRatioData2,rm(GeneName))
View(gcRatioData2)

View(gcRatioData2 %>% group_by(SampleId,Category) %>% summarise(Total=sum(Count)) %>% filter(Category=='ALL'|Category=='TRANS_FIT_EXPECTED'))


gcRatioData3 = gcRatioData2 %>% spread(GeneName,Count,fill=0)
View(gcRatioData3)

print(ggplot(gcRatioData2 %>% 
                filter(Category=='ALL_PERC'|Category=='TRANS_FIT_EXPECTED_PERC'), aes(x=RatioDec, y=Count, fill=Category))
      + geom_bar(stat='identity', position='dodge')
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + facet_wrap(~SampleId))

# p <- ggplot(d, aes(lab, perc, fill = name)) + geom_bar(position = "dodge")

print(ggplot(gcRatioData2 %>% 
                filter(Category=='ALL_PERC'|Category=='TRANS_FIT_EXPECTED_PERC'), aes(x=RatioDec, y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + facet_wrap(~Category))

print(ggplot(gcRatioData2 %>% filter(SampleId=='CPCT02210029T') %>% 
                filter(GeneName=='ALL'|GeneName=='NON_ENRICHED'), aes(x=RatioDec, y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + facet_wrap(~GeneName))

print(ggplot(gcRatioData2 %>% filter(GeneName=='ADJUSTMENTS'), aes(x=RatioDec, y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9)))

print(ggplot(gcRatioData3 %>% group_by(RatioDec) %>% summarise(Count=first(ALL)), aes(x=RatioDec, y=Count))
      + geom_line())

print(ggplot(gcRatioData3, aes(x=RatioDec, y=ALL_PERC))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9)))





## Read Data

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
max(geneOverlaps$RangeEnd - geneOverlaps$RangeStart)

View(geneOverlaps %>% mutate(Length=RangeEnd-RangeStart) %>% filter(Length>1e6))

View(geneOverlaps %>% filter(grepl('RBFOX1',GeneNames)))


View(ensemblGeneData %>% filter(Chromosome==10))

# nrow(genePanel)
View(genePanel)
View(geneOverlaps)

View(ensemblGeneData %>% filter(grepl(GeneName,'ZNRF3;ZNRF3-IT1;ZNRF3-AS1')))

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




#####
## Fragment Lengths


load_frag_length_data<-function(sampleId)
{
   filename = formFilename('~/data/rna/logs/',sampleId,'frag_length.csv')
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

fragLengthsScData = read.csv('~/logs/CPCT02020378T.isf.frag_length.csv')
View(fragLengthsScData)

fragLengthsScData2 = fragLengthsScData %>% gather('ScLength','Count',3:ncol(fragLengthsScData))
View(fragLengthsScData2)

print(ggplot(fragLengthsScData2 %>% filter(FragmentLength>1000), aes(x=FragmentLength,y=Count))
      + geom_line()
      + facet_wrap(~ScLength))
      

print(ggplot(fragLengthsScData %>% filter(FragmentLength<1000), aes(x=FragmentLength))
      # + geom_line(aes(y=Count,color='Count'))
      + geom_line(aes(y=Sc1,color='Sc1'))
      + geom_line(aes(y=Sc2,color='Sc2'))
      + geom_line(aes(y=Sc3,color='Sc3'))
      + geom_line(aes(y=Sc5,color='Sc5'))
      + geom_line(aes(y=Sc10,color='Sc10'))
      + geom_line(aes(y=Sc25,color='Sc25')))

print(ggplot(fragLengthsScData %>% filter(FragmentLength<500), aes(x=FragmentLength))
      # + geom_line(aes(y=Count,color='Count'))
      + geom_bar(stat='identity',aes(y=Sc1),color='red',position='dodge')
      + geom_bar(stat='identity',aes(y=Sc2),color='orange',position='dodge')
      + geom_bar(stat='identity',aes(y=Sc25),color='blue',position='dodge'))

View(tmp)
fragLengthsNoLimit = read.csv('~/logs/CPCT02020378T.isf.frag_length.csv')
fragLengthsScLimit = read.csv('~/logs/CPCT02020378T.isf.frag_length.csv')
View(fragLengthsNoLimit)
View(fragLengthsScLimit)



fragLengthsMerged = merge(fragLengthsNoLimit %>% mutate(NoLimit=Count) %>% select(-Count),fragLengthsScLimit %>% mutate(ScLimited=Count) %>% select(-Count),
                          by='FragmentLength',all.x=T)

View(fragLengthsMerged)

print(ggplot(fragLengthsMerged %>% filter(FragmentLength<600), aes(x=FragmentLength))
      + geom_line(aes(y=NoLimit,color='NoLimit'))
      + geom_line(aes(y=ScLimited,color='ScLimited')))


print(ggplot(fragLengthsNoLimit %>% filter(FragmentLength<300), aes(x=FragmentLength, y=Count))
      + geom_line())

print(ggplot(fragLengthsScLimit %>% filter(FragmentLength<300), aes(x=FragmentLength, y=Count))
      + geom_line())

fragLengthsByGene = read.csv('~/logs/CPCT02020378T.isf.frag_length_by_gene.csv')
View(fragLengthsByGene)

longFragments = read.csv('~/logs/long_fragments.csv')
View(longFragments)

flSamples = c('CPCT02210041T','CPCT02050181T','CPCT02020560T','CPCT02060015T','CPCT02230006T')

fragLengthsCombinedNew = fragLengthsCombined

fragLengthsCombined = data.frame()
for(sampleId in flSamples)
{
   flData = read.csv(formFilename('~/data/rna/tmp/',sampleId,'frag_length.csv'))
   flData$SampleId=sampleId
   fragLengthsCombined = rbind(fragLengthsCombined,flData)
}

View(fragLengthsCombined)
write.csv(fragLengthsCombined,'~/logs/test_samples_frag_lengths.csv',row.names = F,quote = F)

fragLengthsCombinedBoth = merge(fragLengthsCombined,fragLengthsCombinedNew %>% mutate(ScLimitedCount=Count) %>% select(-Count),by=c('SampleId','FragmentLength'),all.x=T)

View(fragLengthsCombinedBoth)
View(fragLengthsCombinedBoth %>% filter(SampleId=='CPCT02020560T'|SampleId=='CPCT02230006T') %>% filter(FragmentLength>70&FragmentLength<80))

print(ggplot(fragLengthsCombined, aes(x=FragmentLength, y=Count))
      + geom_line()
      + xlim(0,1000)
      + facet_wrap(~SampleId))

print(ggplot(fragLengthsCombined %>% filter(FragmentLength<300), aes(x=FragmentLength, y=Count))
      + geom_line()
      + facet_wrap(~SampleId))

print(ggplot(fragLengthsCombinedBoth %>% filter(FragmentLength<300), aes(x=FragmentLength))
      + geom_line(aes(y=Count,color='Count'))
      + geom_line(aes(y=ScLimitedCount,color='ScLimitedCount'))
      + facet_wrap(~SampleId))



geneIds = ensemblGeneData %>% filter(Chromosome==11) # % &GeneStart>=6662885&GeneEnd<=6676967)
View(geneIds)
# 7760576	7760651	7742086	7742159
View(ensemblTransExonData %>% filter(GeneId %in% geneIds$GeneId) %>% filter(!(ExonStart>7760576|ExonEnd<7760651)&!(ExonStart>7742086|ExonEnd<7742159)))
View(ensemblTransExonData %>% filter(GeneId %in% geneIds$GeneId) %>% filter(ExonStart>=6662885&ExonEnd<6676967))
View(ensemblTransExonData %>% filter(ExonStart>=6662885&ExonEnd<6676967))

altSJs = read.csv('~/data/rna/runs/CPCT02020378T.isf.alt_splice_junc.csv')
View(altSJs %>% filter(Chromosome==11&SjStart>6e6&SjEnd<7e6))

# 6676892	6676967	6662885	6662960

View(fragLengthsByGene %>% group_by(Chromosome) %>% count)
View(fragLengthsByGene %>% group_by(Genes) %>% count)
fragLengthsGrouped = fragLengthsByGene %>% group_by(FragmentLength) %>% summarise(TotalCount=sum(Count))

print(ggplot(fragLengths, aes(x=FragmentLength, y=TotalCount))
      + geom_line()
      + scale_x_log10())

print(ggplot(fragLengths, aes(x=FragmentLength, y=TotalCount))
      + geom_line()
      + xlim(0,1000))

print(ggplot(sampleFraglengths %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count)) %>% filter(FragmentLength>=50&FragmentLength<100), aes(x=FragmentLength, y=FragCount))
      + geom_line()
      + xlim(0,1000))

View(sampleFraglengths %>% group_by(FragmentLength) %>% summarise(FragCount=sum(Count)) %>% filter(FragmentLength>=50&FragmentLength<100))



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










