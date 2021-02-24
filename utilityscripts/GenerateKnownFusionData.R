library(dplyr)
library(tidyr)

ensemblPath37 = '~/data/ensembl_hg37/'
ensemblPath38 = '~/data/ensembl_hg38/'
outputDir = '~/data/fusion_ref/'
fileVersion = 'v3'
rawFusionFile = '~/Downloads/HMF Fusion Knowledgebase - Main.tsv'

ensemblGeneData37 = read.csv(paste0(ensemblPath37,'ensembl_gene_data.csv'))
ensemblGeneData38 = read.csv(paste0(ensemblPath38,'ensembl_gene_data.csv'))

View(ensemblGeneData37)
View(ensemblGeneData38)

rawFusionData = read.csv(rawFusionFile,sep='\t')
rawFusionData = rawFusionData %>% arrange(Type)
View(rawFusionData)
View(rawFusionData %>% group_by(Type) %>% count)
View(rawFusionData %>% filter(as.character(FiveGene)!=as.character(FiveGeneRef38)|as.character(ThreeGene)!=as.character(ThreeGeneRef38)) %>%
       select(Type,FiveGene,FiveGeneRef38,ThreeGene,ThreeGeneRef38,everything()))

# write each file version of known_fusion_data for Linx and Isofox
write.csv(rawFusionData %>% select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange,HighImpactPromiscuous,Overrides),
          paste0(outputDir,'known_fusion_data.37_',fileVersion,'.csv'),row.names = F,quote = F)

write.csv(rawFusionData %>% select(Type,FiveGene=FiveGeneRef38,ThreeGene=ThreeGeneRef38,CancerTypes,PubMedId,KnownExonTranscript,
                                   KnownExonUpRange,KnownExonDownRange,HighImpactPromiscuous,Overrides=OverridesRef38),
          paste0(outputDir,'known_fusion_data.38_',fileVersion,'.csv'),row.names = F,quote = F)

# Type,FiveGene,ThreeGene,CancerTypes,PubMedId,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange,HighImpactPromiscuous,Overrides

## GRIPSS BED-PE files (37 and 38)

create_bedpe_file<-function(fusionData,ensemblGeneData,isRef38,outputFile)
{
  # extract Ensembl gene data for known pair genes
  kpBedInfo = merge(fusionData %>% filter(Type == 'KNOWN_PAIR') %>% select(-Overrides),
                    ensemblGeneData %>% select(GeneName,UpChr=Chromosome,UpStrand=Strand,UpGeneStart=GeneStart,UpGeneEnd=GeneEnd),
                    by.x='FiveGene',by.y='GeneName',all.x=T)
  
  # use specified IG gene ranges
  igKnownGenes = fusionData %>% filter(Type == 'IG_KNOWN_PAIR')
  
  # example: IG_RANGE=-1;2;89890568;90274235, or IG_RANGE=-1;2;89890568;90274235 DOWN_GENE_DOWNSTREAM_DISTANCE=500000
  igKnownGenes = fusionData %>% filter(Type == 'IG_KNOWN_PAIR') %>% separate(Overrides,c('IgRangeData','OtherInfo'),sep=' ') %>% select(-OtherInfo)
  igKnownGenes = igKnownGenes %>% mutate(IgRangeData=stri_replace_all_fixed(IgRangeData,'IG_RANGE=',''))
  igKnownGenes = igKnownGenes %>% separate(IgRangeData,c('UpStrand','UpChr','UpGeneStart','UpGeneEnd'),sep=';')
  igKnownGenes = igKnownGenes %>% mutate(UpStrand=as.numeric(UpStrand),UpGeneStart=as.numeric(UpGeneStart),UpGeneEnd=as.numeric(UpGeneEnd))
  
  kpBedInfo = rbind(kpBedInfo,igKnownGenes)
  
  kpBedInfo = merge(kpBedInfo,
                    ensemblGeneData %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),
                    by.x='ThreeGene',by.y='GeneName',all.x=T)
  
  View(kpBedInfo)
  # str(kpBedInfo)
  
  preGeneBuffer=10e3
  
  kpBedInfo = kpBedInfo %>%
    mutate(
      UpChr = factor(UpChr, levels = c(1:22,'X','Y'), ordered = T),
      DownChr = factor(DownChr, levels = c(1:22,'X','Y'), ordered = T),
      Name=paste(FiveGene,ThreeGene,sep='-'),
      ## NOTE THAT WE SUBTRACT 1 FROM START FOR BEDPE FORMAT
      Start1=ifelse(UpStrand==1,UpGeneStart-preGeneBuffer,UpGeneStart) - 1,
      End1=ifelse(UpStrand==1,UpGeneEnd,UpGeneEnd+preGeneBuffer),
      Start2=ifelse(DownStrand==1,DownGeneStart-preGeneBuffer,DownGeneStart) - 1,
      End2=ifelse(DownStrand==1,DownGeneEnd,DownGeneEnd+preGeneBuffer),
      Strand1=ifelse(UpStrand==1,'+','-'),
      ## NOTE WE REVERSE STRAND2
      Strand2=ifelse(DownStrand==1,'-','+'),
      Score=0)
  
  # pre-pend 'chr' for HG38
  if(isRef38)
  {
    kpBedInfo = kpBedInfo %>% mutate(UpChr=paste('chr',UpChr,sep=''),DownChr=paste('chr',DownChr,sep=''))
  }

  bedpe = kpBedInfo %>% mutate(
    StartIsUp = UpChr < DownChr | (UpChr == DownChr & Start1 < Start2),
    StartChr = if_else(StartIsUp, UpChr, DownChr),
    StartPositionStart = ifelse(StartIsUp, Start1, Start2),
    StartPositionEnd = ifelse(StartIsUp, End1, End2),
    StartStrand = ifelse(StartIsUp, Strand1, Strand2),
    EndChr = if_else(!StartIsUp, UpChr, DownChr),
    EndPositionStart = ifelse(!StartIsUp, Start1, Start2),
    EndPositionEnd = ifelse(!StartIsUp, End1, End2),
    EndStrand = ifelse(!StartIsUp, Strand1, Strand2)) %>%
    select(StartChr, StartPositionStart, StartPositionEnd, EndChr, EndPositionStart, EndPositionEnd, Name, Score, StartStrand, EndStrand, UpChr, DownChr) %>%
    arrange(StartChr, StartPositionStart)
  
  # View(bedpe)
  write.table(bedpe,outputFile,row.names=F,col.names=F,sep="\t",quote=F)  
}

create_bedpe_file(rawFusionData %>% select(Type,FiveGene,ThreeGene,Overrides),ensemblGeneData37,F,
                  paste0(outputDir,'known_fusions.37_',fileVersion,'.bedpe'))

create_bedpe_file(rawFusionData %>% select(Type,FiveGene=FiveGeneRef38,ThreeGene=ThreeGeneRef38,Overrides=OverridesRef38),ensemblGeneData38,T,
                  paste0(outputDir,'known_fusions.38_',fileVersion,'.bedpe'))

## GRIPSS BED files (37 and 38) for high-confidence promiscuous gene exon ranges
ensemblTransExonData37 = read.csv(paste0(ensemblPath37,'ensembl_trans_exon_data.csv'))
ensemblTransExonData38 = read.csv(paste0(ensemblPath38,'ensembl_trans_exon_data.csv'))
View(ensemblTransExonData37 %>% filter(TransName=='ENST00000288602'))


create_bed_file<-function(fusionData,ensemblGeneData,ensemblTransExonData,isRef38,outputFile)
{
  highConfRegions = fusionData %>% filter(Type %in% c('PROMISCUOUS_3','PROMISCUOUS_5')&KnownExonTranscript!='')
  
  highConfRegions = highConfRegions %>% mutate(ExonInfo=ifelse(Type=='PROMISCUOUS_5',as.character(KnownExonUpRange),as.character(KnownExonDownRange)),
                                               GeneName=ifelse(Type=='PROMISCUOUS_5',as.character(FiveGene),as.character(ThreeGene)))
  
  highConfRegions = merge(highConfRegions,ensemblGeneData %>% select(GeneName,Chromosome),by='GeneName',all.x=T)
  
  highConfRegions = highConfRegions %>% separate(ExonInfo,c('ExonRankUp','ExonRankDown'),sep=';')
  
  highConfRegions = highConfRegions %>% mutate(ExonRankUpPrev=ifelse(Type=='PROMISCUOUS_3',as.numeric(ExonRankUp)-1,as.numeric(ExonRankUp)),
                                               ExonRankDownNext=ifelse(Type=='PROMISCUOUS_5',as.numeric(ExonRankDown)+1,as.numeric(ExonRankDown)))
  
  # View(highConfRegions %>% select(Type,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange,ExonRankUp,ExonRankUpPrev,ExonRankDown,ExonRankDownNext,everything()))
  
  
  highConfRegions2 = merge(highConfRegions,ensemblTransExonData %>% select(KnownExonTranscript=TransName,Strand,ExonRankUp=ExonRank,ExonUpStart=ExonStart,ExonUpEnd=ExonEnd),
                           by=c('KnownExonTranscript','ExonRankUp'),all.x=T)
  
  highConfRegions2 = merge(highConfRegions2,ensemblTransExonData %>% select(KnownExonTranscript=TransName,ExonRankUpPrev=ExonRank,ExonUpPrevStart=ExonStart,ExonUpPrevEnd=ExonEnd),
                           by=c('KnownExonTranscript','ExonRankUpPrev'),all.x=T)
  
  highConfRegions2 = merge(highConfRegions2,ensemblTransExonData %>% select(KnownExonTranscript=TransName,ExonRankDown=ExonRank,ExonDownStart=ExonStart,ExonDownEnd=ExonEnd),
                           by=c('KnownExonTranscript','ExonRankDown'),all.x=T)
  
  highConfRegions2 = merge(highConfRegions2,ensemblTransExonData %>% select(KnownExonTranscript=TransName,ExonRankDownNext=ExonRank,ExonDownNextStart=ExonStart,ExonDownNextEnd=ExonEnd),
                           by=c('KnownExonTranscript','ExonRankDownNext'),all.x=T)
  
  # View(highConfRegions2 %>% select(Type,Strand,ExonRankUp,ExonRankDown,ExonUpStart,ExonUpEnd,ExonUpPrevStart,ExonUpPrevEnd,ExonDownStart,ExonDownEnd,ExonDownNextStart,ExonDownNextEnd,everything()))
  
  # correct for regions ending at the last exon
  highConfRegions2 = highConfRegions2 %>% mutate(ExonDownNextStart=ifelse(is.na(ExonDownNextStart),ExonDownEnd,ExonDownNextStart),
                                                 ExonDownNextEnd=ifelse(is.na(ExonDownNextEnd),ExonDownEnd,ExonDownNextEnd))
  
  highConfRegions2 = highConfRegions2 %>% mutate(RangeStart=ifelse(Type=='PROMISCUOUS_5'&Strand==1,ExonUpEnd,ifelse(Type=='PROMISCUOUS_5'&Strand==-1,ExonDownNextEnd,
                                                                                                                    ifelse(Type=='PROMISCUOUS_3'&Strand==1,ExonUpPrevEnd,ExonDownEnd))),
                                                 RangeEnd=ifelse(Type=='PROMISCUOUS_5'&Strand==1,ExonDownNextStart,ifelse(Type=='PROMISCUOUS_5'&Strand==-1,ExonUpStart,
                                                                                                                          ifelse(Type=='PROMISCUOUS_3'&Strand==1,ExonDownStart,ExonUpPrevStart))))
  
  # View(highConfRegions2 %>% select(Type,Strand,ExonRankUp,ExonRankDown,RangeStart,RangeEnd,ExonUpStart,ExonUpEnd,ExonUpPrevStart,ExonUpPrevEnd,ExonDownStart,ExonDownEnd,ExonDownNextStart,ExonDownNextEnd,everything()))
  bedData = highConfRegions2 %>% select(Chromosome,RangeStart,RangeEnd,GeneName,Strand) %>% 
    mutate(ChrSorted=factor(Chromosome,levels=c(1:22,'X','Y'), ordered = T)) %>% arrange(ChrSorted,RangeStart)
  
  if(isRef38)
  {
    bedData = bedData %>% mutate(Chromosome=paste('chr',Chromosome,sep=''))
  }
  
  bedData = bedData %>% mutate(Score=0,Strand=ifelse(Strand=='1','+','-'))
  
  write.table(bedData %>% select(Chromosome,RangeStart,RangeEnd,GeneName,Score,Strand),outputFile,row.names=F,col.names=F,sep="\t",quote=F)  
}

create_bed_file(rawFusionData %>% select(Type,FiveGene,ThreeGene,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange),
                ensemblGeneData37,ensemblTransExonData37,F,paste0(outputDir,'known_fusions.37_',fileVersion,'.bed'))

create_bed_file(rawFusionData %>% select(Type,FiveGene=FiveGeneRef38,ThreeGene=ThreeGeneRef38,KnownExonTranscript=KnownExonTranscriptRef38,KnownExonUpRange,KnownExonDownRange),
                ensemblGeneData38,ensemblTransExonData38,T,paste0(outputDir,'known_fusions.38_',fileVersion,'.bed'))

View(highConfRegions2 %>% select(Type,GeneName,Chromosome,Strand,ExonRankUp,ExonRankDown,RangeStart,RangeEnd,ExonUpStart,ExonUpEnd,ExonUpPrevStart,ExonUpPrevEnd,ExonDownStart,ExonDownEnd,ExonDownNextStart,ExonDownNextEnd,everything()))
write.csv(highConfRegions2 %>% select(Type,GeneName,Chromosome,Strand,ExonRankUp,ExonRankDown,RangeStart,RangeEnd,ExonUpStart,ExonUpEnd,ExonUpPrevStart,ExonUpPrevEnd,ExonDownStart,ExonDownEnd,ExonDownNextStart,ExonDownNextEnd,everything()),
          '~/logs/fusion_bed_working.csv',row.names = F,quote = F)
