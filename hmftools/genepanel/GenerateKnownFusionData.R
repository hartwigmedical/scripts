library(dplyr)
library(tidyr)

ensemblPath37 = '~/data/ensembl_hg37/'
ensemblPath38 = '~/data/ensembl_hg38/'

# following switch to use HGNC for gene names
ensemblPath37 = '~/data/ensembl_db/ensembl_37_104/'
ensemblPath38 = '~/data/ensembl_db/ensembl_38_104/'
outputDir = '~/data/fusion_ref/'
fileVersion = 'v7'

# export from Google sheet and move & rename
# mv ~/Downloads/HMF\ Fusion\ Knowledgebase\ -\ Main.tsv ~/data/fusion_ref/HmfKnownFusionSheet_20220801.tsv
#rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet.tsv')
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20210419.tsv')
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20210927.tsv') # v3
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20211129.tsv') # v4
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20220606.tsv') # v5
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20220801.tsv') # v6
rawFusionFile = paste0(outputDir,'HmfKnownFusionSheet_20220901.tsv') # v7 - fixes up IG 'chr' prefix for HG38

ensemblGeneData37 = read.csv(paste0(ensemblPath37,'ensembl_gene_data.csv'))
ensemblGeneData38 = read.csv(paste0(ensemblPath38,'ensembl_gene_data.csv'))

View(ensemblGeneData37)
View(ensemblGeneData38)

rawFusionData = read.csv(rawFusionFile,sep='\t')
rawFusionData = rawFusionData %>% arrange(Type,FiveGene,ThreeGene)
View(rawFusionData)
View(rawFusionData %>% group_by(Type) %>% count)
View(rawFusionData %>% filter(as.character(FiveGene)!=as.character(FiveGeneRef38)|as.character(ThreeGene)!=as.character(ThreeGeneRef38)) %>%
       select(Type,FiveGene,FiveGeneRef38,ThreeGene,ThreeGeneRef38,everything()))

# confirm all genes are in Ensembl

nrow(rawFusionData %>% filter(FiveGene!=''&!(FiveGene %in% c('IGH','IGL','IGK')) & !(FiveGene %in% ensemblGeneData37$GeneName))) # none
nrow(rawFusionData %>% filter(ThreeGene!='' & !(ThreeGene %in% ensemblGeneData37$GeneName))) # none

nrow(rawFusionData %>% filter(FiveGene!=''&!(FiveGene %in% c('IGH','IGL','IGK')) & !(FiveGene %in% ensemblGeneData38$GeneName))) # none
nrow(rawFusionData %>% filter(ThreeGene!='' & !(ThreeGene %in% ensemblGeneData38$GeneName))) # none


# write each file version of known_fusion_data for Linx and Isofox
write.csv(rawFusionData %>% select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange,HighImpactPromiscuous,Overrides),
          paste0(outputDir,'known_fusion_data.37_',fileVersion,'.csv'),row.names = F,quote = F)

write.csv(rawFusionData %>% select(Type,FiveGene=FiveGeneRef38,ThreeGene=ThreeGeneRef38,CancerTypes,PubMedId,KnownExonTranscript,
                                   KnownExonUpRange,KnownExonDownRange,HighImpactPromiscuous,Overrides=OverridesRef38),
          paste0(outputDir,'known_fusion_data.38_',fileVersion,'.csv'),row.names = F,quote = F)

## GRIPSS BED-PE files (37 and 38)

create_bedpe_file<-function(fusionData,ensemblGeneData,isRef38,outputFile)
{
  # extract Ensembl gene data for known pair genes
  kpBedInfo = merge(fusionData %>% filter(Type == 'KNOWN_PAIR') %>% select(-Overrides),
                    ensemblGeneData %>% select(GeneName,UpChr=Chromosome,UpStrand=Strand,UpGeneStart=GeneStart,UpGeneEnd=GeneEnd),
                    by.x='FiveGene',by.y='GeneName',all.x=T)
  
  kpBedInfo = merge(kpBedInfo,
                    ensemblGeneData %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),
                    by.x='ThreeGene',by.y='GeneName',all.x=T)
  
  # use specified IG gene ranges and write both orientations
  igKnownGenes = fusionData %>% filter(Type == 'IG_KNOWN_PAIR')
  
  # example: IG_RANGE=-1;2;89890568;90274235, or IG_RANGE=-1;2;89890568;90274235 DOWN_GENE_DOWNSTREAM_DISTANCE=500000
  igKnownGenes = fusionData %>% filter(Type == 'IG_KNOWN_PAIR') %>% separate(Overrides,c('IgRangeData','DownstreamDistance'),sep=' ')
  igKnownGenes = igKnownGenes %>% mutate(IgRangeData=stri_replace_all_fixed(IgRangeData,'IG_RANGE=',''))
  igKnownGenes = igKnownGenes %>% separate(IgRangeData,c('UpStrand','UpChr','UpGeneStart','UpGeneEnd'),sep=';')
  igKnownGenes = igKnownGenes %>% mutate(UpStrand=as.numeric(UpStrand),UpGeneStart=as.numeric(UpGeneStart),UpGeneEnd=as.numeric(UpGeneEnd))
  
  # strip chr prefix from IG genes for 38
  igKnownGenes = igKnownGenes %>% mutate(UpChr=stri_replace_all_fixed(UpChr,'chr',''))
  
  igKnownGenes = merge(igKnownGenes,
                       ensemblGeneData %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),
                       by.x='ThreeGene',by.y='GeneName',all.x=T)
  
  igKnownGenes = igKnownGenes %>% mutate(DownstreamDistance=ifelse(is.na(DownstreamDistance),'0',stri_replace_all_fixed(DownstreamDistance,'DOWN_GENE_DOWNSTREAM_DISTANCE=','')),
                                         DownstreamDistance=ifelse(is.na(DownstreamDistance),0,as.numeric(DownstreamDistance)))
  
  igKnownGenesDownstream = igKnownGenes %>% filter(DownstreamDistance>0) %>% 
    mutate(DownGeneStartNew=ifelse(DownStrand==1,DownGeneEnd,DownGeneStart-DownstreamDistance),
           DownGeneEndNew=ifelse(DownStrand==1,DownGeneEnd+DownstreamDistance,DownGeneStart),
           DownGeneStart=DownGeneStartNew,
           DownGeneEnd=DownGeneEndNew,
           DownStrand=ifelse(DownStrand==1,-1,1)) %>%
    select(-DownGeneStartNew,-DownGeneEndNew)
  
  igKnownGenes = rbind(igKnownGenes,igKnownGenesDownstream)
  
  igKnownGenesReverseOrient = igKnownGenes %>% mutate(UpStrand=ifelse(UpStrand==1,-1,1))
  
  igKnownGenes = rbind(igKnownGenes,igKnownGenesReverseOrient)
  # View(igKnownGenes)
  
  kpBedInfo = rbind(kpBedInfo,igKnownGenes %>% select(-DownstreamDistance))

  # View(kpBedInfo)

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
      ## REVERSE STRAND2 since for the downstream genes the orientation is opposite to upstream (ie +ve strand = -ve orientation and vice versa)
      Strand2=ifelse(DownStrand==1,'-','+'),
      Score=0)
  
  kpBedInfo = kpBedInfo %>% mutate(
    StartIsUp = UpChr < DownChr | (UpChr == DownChr & Start1 < Start2),
    StartChr = if_else(StartIsUp, UpChr, DownChr),
    StartPositionStart = ifelse(StartIsUp, Start1, Start2),
    StartPositionEnd = ifelse(StartIsUp, End1, End2),
    StartStrand = ifelse(StartIsUp, Strand1, Strand2),
    EndChr = if_else(!StartIsUp, UpChr, DownChr),
    EndPositionStart = ifelse(!StartIsUp, Start1, Start2),
    EndPositionEnd = ifelse(!StartIsUp, End1, End2),
    EndStrand = ifelse(!StartIsUp, Strand1, Strand2)) %>%
    arrange(StartChr, StartPositionStart)
  
  # pre-pend 'chr' for HG38
  if(isRef38)
  {
    kpBedInfo = kpBedInfo %>% mutate(UpChr=paste('chr',UpChr,sep=''),DownChr=paste('chr',DownChr,sep=''),
                                     StartChr=paste('chr',StartChr,sep=''),EndChr=paste('chr',EndChr,sep=''))
  }

  View(kpBedInfo)
  
  bedpe = kpBedInfo %>% select(StartChr, StartPositionStart, StartPositionEnd, EndChr, EndPositionStart, EndPositionEnd, Name, Score, StartStrand, EndStrand, UpChr, DownChr)
  write.table(bedpe,outputFile,row.names=F,col.names=F,sep="\t",quote=F)  
}

create_bedpe_file(rawFusionData %>% select(Type,FiveGene,ThreeGene,Overrides),ensemblGeneData37,F,
                  paste0(outputDir,'known_fusions.37_',fileVersion,'.bedpe'))

create_bedpe_file(rawFusionData %>% select(Type,FiveGene=FiveGeneRef38,ThreeGene=ThreeGeneRef38,Overrides=OverridesRef38),ensemblGeneData38,T,
                  paste0(outputDir,'known_fusions.38_',fileVersion,'.bedpe'))

# deployment
# cp known_fusion_data.37_v7.csv ~/hmf/repos/common-resources-public/fusions/37/known_fusion_data.37.csv
# cp known_fusion_data.38_v7.csv ~/hmf/repos/common-resources-public/fusions/38/known_fusion_data.38.csv
# cp known_fusions.37_v7.bedpe ~/hmf/repos/common-resources-public/fusions/37/known_fusions.37.bedpe
# cp known_fusions.38_v7.bedpe ~/hmf/repos/common-resources-public/fusions/38/known_fusions.38.bedpe


