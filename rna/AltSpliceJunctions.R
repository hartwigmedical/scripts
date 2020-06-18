## RNA Alternative Splice Junction investigations



## Isofox basic data validation
altSJs = read.csv(formFilename('~/logs/',sampleId,'alt_splice_junc.csv'))
altSJs = read.csv(formFilename('~/logs/',sampleId,'gp_overs.alt_splice_junc.csv'))
View(altSJs)
nrow(altSJs)

tmpAltSJs = read.csv(formFilename('~/logs/',sampleId,'alt_splice_junc.csv'))
View(tmpAltSJs)

View(ensemblTransExonData %>% filter(TransName=='ENST00000242057'))
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000106546'))
View(ensemblTransExonData %>% filter(TransName=='ENST00000242057'))

altSJs = altSJs %>% mutate(StartSeq=stri_sub(StartBases,3,4),EndSeq=stri_sub(EndBases,9,10))
View(altSJs %>% filter(Strand==1))
View(altSJs %>% filter(Strand==1) %>% group_by(Seq=paste(StartSeq,EndSeq,sep='_'),Type) %>% count %>% spread(Type,n,fill=0))
View(altSJs %>% filter(FragCount>StartDepth|FragCount>EndDepth))

View(altSJs %>% filter(StartSeq=='CT'&EndSeq=='AC'))






## AHR fusions of exons 7 to 10
utSampleData = read.csv('~/data/urinary_tract_sample_data.csv')

ahrFusionSamples = utSampleData %>% filter(AHRMutationType=='FUSION')
View(ahrFusionSamples)
View(ahrFusionSamples %>% group_by(SampleId) %>% count)
View(ahrFusionSamples %>% filter(SampleId %in% ahrFusionSampleIds))

utSamples = rnaSampleData %>% filter(CancerType=='Urinary tract'&ReadLength==76)

ahrFusionSampleIds = c('CPCT02020787T','CPCT02330011T','CPCT02020840T','CPCT02020867T','CPCT02070272T','CPCT02020642TII','CPCT02020894T',
                       'CPCT02020927T','DRUP01030012T','CPCT02050124T','CPCT02010944T',
                       'CPCT02010944TII','CPCT02010944TIII','CPCT02050103T','CPCT02020619T','CPCT02020618T','CPCT02020629T','CPCT02070044T')

utSamples = utSamples %>% mutate(HasAHRFusion=(SampleId %in% ahrFusionSampleIds))
View(utSamples %>% group_by(HasAHRFusion) %>% count)


## Cohort Analysis
write.csv(utSamples %>% filter(!grepl('TII',SampleId)&!grepl('TIII',SampleId)) %>% 
            mutate(CohortName=ifelse(HasAHRFusion,'CohortA','CohortB')) %>% select(SampleId,CancerType,CohortName),
          '~/data/rna/dl_ut_138_samples.csv', row.names = F,quote = F)

View(utSamples %>% filter(SampleId %in% c('CPCT02010944T','CPCT02010944TII','CPCT02010944TIII','CPCT02020618T','CPCT02020619T','CPCT02020629T','CPCT02020642T','CPCT02020701T',
                                          'CPCT02020787T','CPCT02020840T','CPCT02020867T','CPCT02020894T','CPCT02050103T',
                                          'CPCT02050124T','CPCT02070044T','CPCT02070272T','CPCT02070294T','CPCT02120168T','CPCT02330011T','CPCT02330103T','DRUP01030012T')))

View(completedSamples)

write.csv(rnaSampleData %>% filter(SampleId %in% completedSamples$SampleId) %>% mutate(CohortName=ifelse(CancerType=='Urinary tract','CohortA','CohortB')) %>% 
       select(SampleId,CancerType,CohortName),'~/data/rna/samples/dl_ut_vs_rest.csv',row.names = F,quote = F)

altCohortUtAhr = read.csv('~/data/rna/logs/isfox_ut_ahr.alt_sj_cohort.csv')

altCohortUtAhr = merge(ensemblGeneData %>% select(GeneId,GeneName,Strand),altCohortUtAhr,by='GeneId',all.y=T)
View(altCohortUtAhr)
write.csv(altCohortUtAhr,'~/data/rna/logs/isfox_ut_ahr.alt_sj_cohort.csv',row.names = F,quote = F)


## CRC with APC or not
apcSamples = c('CPCT02010730T','CPCT02020695T','CPCT02120162T','CPCT02220059T','CPCT02040165T','CPCT02050181T','DRUP01070017T','CPCT02020521T',
               'CPCT02300017T','CPCT02230027T','CPCT02080131T','CPCT02070066T','CPCT02040100T','DRUP01090008T')


View(crcSamples %>% filter(SampleId %in% apcSamples))

crcApcSamples = crcSamples %>% mutate(CohortName=ifelse(SampleId %in% apcSamples,'CohortA','CohortB'))
View(crcApcSamples)
write.csv(crcApcSamples %>% select(SampleId,CancerType,CohortName),'~/data/rna/dl_crc_apc.csv',quote = F, row.names = F)

altCohortCrcApc = read.csv('~/data/rna/logs/isfox_crc_apc.alt_sj_cohort.csv')

altCohortCrcApc = merge(ensemblGeneData %>% select(GeneId,GeneName,Strand),altCohortCrcApc,by='GeneId',all.y=T)
View(altCohortCrcApc)
write.csv(altCohortCrcApc,'~/data/rna/logs/isfox_crc_vs_ut.alt_sj_cohort.csv',row.names = F,quote = F)
  

# CRC vs UT
utSamples = rnaSampleData %>% filter(CancerType=='Urinary tract'&ReadLength==76)
View(utSamples)
View(rnaSampleData)


crcSamples = read.csv('~/data/rna/gcp_batch_isofox_crc_24.txt',header = F)
crcSamples = rbind(crcSamples,read.csv('~/data/rna/gcp_batch_isofox_crc_100_b2.txt',header = F))
colnames(crcSamples) = c('SampleId','ReadLength')
crcSamples$CancerType = 'Colon/Rectum'
View(crcSamples)

crcVsUtSamples = rbind(crcSamples %>% select(SampleId,CancerType),
                       utSamples %>% filter(!grepl('TII',SampleId)&!grepl('TIII',SampleId)) %>% select(SampleId,CancerType))

crcVsUtSamples = crcVsUtSamples %>% mutate(CohortName=ifelse(CancerType=='Urinary tract','CohortA','CohortB'))
View(crcVsUtSamples)
write.csv(crcVsUtSamples, '~/data/rna/dl_crc_vs_ut.csv',quote = F, row.names = F)

altCohortCrcVsUt = read.csv('~/data/rna/logs/isfox_crc_vs_ut.alt_sj_reoccurring.csv')

altCohortUtVsRest = read.csv('~/data/rna/samples/isfox_ut_vs_rest.alt_sj_cohort.csv')
View(altCohortUtVsRest)

altCohortCrcVsUt = merge(ensemblGeneData %>% select(GeneId,GeneName,Strand),altCohortCrcVsUt,by='GeneId',all.y=T)
View(altCohortCrcVsUt)
write.csv(altCohortCrcVsUt,'~/data/rna/logs/isfox_crc_vs_ut.alt_sj_reoccurring.csv',row.names = F,quote = F)

View(altCohortCrcVsUt %>% filter(CohortAWithAltSJ>ExpVal))
View(altCohortCrcVsUt %>% filter(FetProb<0.0001))
colnames(altCohortCrcVsUt)
write.csv(altCohortCrcVsUt %>% filter(FetProb<0.01),'~/data/rna/logs/isfox_crc_vs_ut.alt_sj_prob_0.01.csv',row.names = F, quote = F)


altSJReo = read.csv('~/data/rna/logs/isofox_alt_sj_reoccurring.csv')
altSJReo = read.csv('~/logs/isofox_alt_sj_reoccurring.csv')
View(altSJReo)
View(altSJReo %>% group_by(SampleCount) %>% count)
View(altSJReo %>% filter(WithMutationCount>2*NoMutationCount))
View(altSJReo %>% filter(WithMutationCount>2*NoMutationCount))

write.csv(altSJReo %>% filter(FetProb<0.01), '~/data/rna/logs/isofox_alt_sj_reoccurring_prob_0.01.csv',quote = F, row.names = F)


altSJProbs = altSJReo %>% filter(WithMutationCount>NoMutationCount)

altSJProbs = altSJReo %>% filter(WithMutationCount>=2&WithMutationCount>NoMutationCount) %>% 
  mutate(ScAll=138,
         ScWithMut=16,
         ScWithAltSJ=SampleCount,
         ScWithAltSJNoMut=NoMutationCount,
         ScNoAltSJWithMut=ScWithMut-WithMutationCount,
         SvNoAltSJNoMut=ScAll-WithMutationCount-ScNoAltSJWithMut-ScWithAltSJNoMut,
         ExpCount=round(ScWithMut * ScWithAltSJ / ScAll, 2))

altSJProbs$FetProb = 1
colIndex = ncol(altSJProbs)
for(i in 1:nrow(altSJProbs))
{
  asjData = altSJProbs[i,]
  fMatrix = rbind(c(asjData$WithMutationCount,asjData$ScNoAltSJWithMut),c(asjData$ScWithAltSJNoMut,asjData$SvNoAltSJNoMut))
  altSJProbs[i,colIndex] = fisher.test(fMatrix, alternative="greater")$p.value
}

altSJProbs = merge(altSJProbs,ensemblGeneData %>% select(GeneId,GeneName),by='GeneId',all.x=T)
View(altSJProbs %>% select(GeneName,FetProb,ExpCount,SampleCount,WithMutationCount,Type,SjStart,SjEnd,MaxFragCount,AvgFragCount,everything()))


g1Name = 'AltSJ'
g2Name = 'AHRFusion'
scAll = 138
scWithG1 = 7 # with alt SJ
scWithG2 = 16 # with AHR 7-10 fusion
scWithG1WithG2 = 6

scWithG1NoG2 = scWithG1 - scWithG1WithG2
scNoG1WithG2 = scWithG2 - scWithG1WithG2
scNoG1NoG2 = scAll - scWithG1NoG2 - scNoG1WithG2 - scWithG1WithG2

fishMatrix = rbind(c(scWithG1WithG2,scNoG1WithG2), c(scWithG1NoG2,scNoG1NoG2))

expected = round(scWithG1 * scWithG2 / scAll, 2)
fetProb = fisher.test(fishMatrix, alternative="greater")$p.value

print(sprintf("all=%d with%s=%d with%s=%d both=%d neither=%d with%sNo%s=%d no%sWith%s=%d", 
              scAll,g1Name,scWithG1,g2Name,scWithG2,scWithG1WithG2,scNoG1NoG2,
              g1Name,g2Name,scWithG1NoG2,g1Name,g2Name,scNoG1WithG2))

print(sprintf("expected=%.2f prob=%f",expected, fetProb))



## Manual analysis from raw files

asjSampleData = data.frame()

for(sampleId in sampleList)
{
  asjData = read.csv(formFilename('~/data/rna/logs/',sampleId,'alt_splice_junc.csv'))
  sprintf("loaded %d records for sample(%s)", nrow(asjData), sampleId)
  asjData$SampleId = sampleId
  asjSampleData = rbind(asjSampleData,asjData)
}

write.csv(rnaSampleData %>% filter(CancerType=='Urinary tract'&ReadLength==76) %>% select(SampleId),'~/data/rna/ut_76_samples.txt',row.names = F,quote = F)

View(asjSampleData)
View(head(asjSampleData,100))
nrow(asjSampleData)
View(asjSampleData %>% group_by(SampleId) %>% count)

# repeated locations
asjLocations = asjSampleData %>% group_by(GeneId,GeneName,Chromosome,SjStart,SjEnd,Type,StartBases,EndBases) %>% 
  summarise(SampleCount=n(),
            MaxFrags=max(FragCount),
            AvgFrags=mean(FragCount),
            TotalFrags=sum(FragCount),
            AvgVAF=mean(FragCount/(pmax(StartDepth,EndDepth))))


View(asjLocations %>% filter(SampleCount>=5))
nrow(asjLocations)

# unique SJs with high support
View(asjLocations %>% filter(SampleCount==1) %>% group_by(FragCountBucket=2**round(log(TotalFrags,2))) %>% count)
View(asjLocations %>% group_by(SampleCount) %>% count)
View(asjLocations %>% filter(SampleCount==1&TotalFrags>=10) %>% group_by(Type) %>% count)
View(asjLocations %>% filter(SampleCount==1&TotalFrags>=50))

# intronixc SJs by context
intronicSJs = asjSampleData %>% filter(Type=='INTRONIC'&FragCount>=10) %>% mutate(StartSeq=stri_sub(StartBases,3,4),EndSeq=stri_sub(EndBases,9,10))
nrow(intronicSJs)
View(intronicSJs %>% group_by(Strand,Seq=paste(StartSeq,EndSeq,sep='_'),GeneId,GeneName) %>% 
       summarise(SampleCount=n(),
                 MaxFrags=max(FragCount),
                 AvgFrags=mean(FragCount),
                 TotalFrags=sum(FragCount),
                 AvgVAF=mean(FragCount/(pmax(StartDepth,EndDepth)))))


#####% 
## Retained Introns
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','retained_intron.csv'))
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','gp.retained_intron.csv'))
retIntrons = read.csv(formFilename('~/logs/','CPCT02020378T','gp_overs.retained_intron.csv'))
View(retIntrons)
View(retIntrons %>% filter(FragCount>TotalDepth))
View(retIntrons %>% filter(GeneId=='ENSG00000244203'))



#####
## Splice Variant Matching
svMatches = read.csv('~/data/rna/logs/isofox_splice_variant_matching.csv')
View(svMatches)
View(svMatches %>% filter(SampleId=='CPCT02010663T'))

tmpSvMatches = read.csv('~/data/rna/logs/isofox_splice_variant_matching.csv')
View(tmpSvMatches)

tmpAltSJs = read.csv(formFilename('~/logs/','CPCT02010663T','alt_splice_junc.csv'))
View(tmpAltSJs)

tmpAltSJs2 = read.csv(formFilename('~/data/rna/gcp_data/','CPCT02010663T','alt_splice_junc.csv'))
View(tmpAltSJs2)

tmpAltSJs3 = read.csv(formFilename('~/data/rna/runs/','CPCT02010663T','alt_splice_junc.csv'))
View(tmpAltSJs3)

View(tmpAltSJs %>% filter(Chromosome==1&(abs(SjStart-1248167)<1e4|abs(SjEnd-1248167)<1e4)))
View(tmpAltSJs %>% filter(GeneId=='ENSG00000127054'))

View(ensemblGeneData)
View(ensemblTransExonData)


## Validation
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
