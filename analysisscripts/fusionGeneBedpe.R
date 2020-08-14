# HG37
geneFile = '~/data/sv/ensembl_gene_data.csv'
outputFile = '~/data/sv/known_fusion_pairs.hg37.bedpe'
knownPairs = read.csv('~/data/sv/known_fusion_data.csv') %>% filter(Type == 'KNOWN_PAIR')
View(knownPairs)

# HG38 
geneFile = '~/data/sv/ensembl_hg38/ensembl_gene_data.csv'
outputFile = '~/data/sv/known_fusion_pairs.hg38.bedpe'
knownPairs = read.csv('~/data/sv/hg38//known_fusion_data.csv') %>% filter(Type == 'KNOWN_PAIR')

#View(knownPairs)

ensemblGeneInfo = read.csv(geneFile)
#ensemblGeneInfo = ensemblGeneData
#ensemblGeneInfo = ensemblGeneData38New

kpBedInfo = merge(knownPairs,ensemblGeneInfo %>% select(GeneName,UpChr=Chromosome,UpStrand=Strand,UpGeneStart=GeneStart,UpGeneEnd=GeneEnd),by.x='FiveGene',by.y='GeneName',all.x=T)
kpBedInfo = merge(kpBedInfo,ensemblGeneInfo %>% select(GeneName,DownChr=Chromosome,DownStrand=Strand,DownGeneStart=GeneStart,DownGeneEnd=GeneEnd),by.x='ThreeGene',by.y='GeneName',all.x=T)

rm(ensemblGeneInfo)

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

# run for HG38 only
kpBedInfo = kpBedInfo %>% mutate(UpChr=paste('chr',UpChr,sep=''),DownChr=paste('chr',DownChr,sep=''))


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
write.table(bedpe, file = outputFile, row.names = F, col.names = F, sep = "\t", quote = F)
