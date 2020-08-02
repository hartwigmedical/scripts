geneFile = '~/hmf/resources/ensembl_gene_data.hg19.csv'
outputFile = '~/hmf/resources/KnownFusionPairs.hg19.bedpe'
#geneFile = '~/hmf/resources/ensembl_gene_data.hg38.csv'
#outputFile = '~/hmf/resources/KnownFusionPairs.hg38.bedpe'

ensemblGeneData = read.csv(geneFile)
knownPairs = read.csv('~/hmf/resources/known_fusion_data.csv') %>% filter(Type == 'KNOWN_PAIR')
#View(knownPairs)

kpBedInfo = merge(knownPairs,ensemblGeneData %>% select(GeneName,UpChr=Chromosome,UpStrand=Strand,UpGeneStart=GeneStart,UpGeneEnd=GeneEnd),by.x='FiveGene',by.y='GeneName',all.x=T)
kpBedInfo = merge(kpBedInfo,ensemblGeneData %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),by.x='ThreeGene',by.y='GeneName',all.x=T)

str(kpBedInfo)

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

write.table(bedpe, file = outputFile, row.names = F, col.names = F, sep = "\t", quote = F)
