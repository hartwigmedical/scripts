

ensemblGeneData = read.csv('~/hmf/resources/ensembl_gene_data.csv')
knownPairs = read.csv('~/hmf/resources/known_fusion_data.csv') %>% filter(Type == 'KNOWN_PAIR')
View(knownPairs)

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

write.table(bedpe, file = "~/hmf/resources/gridss_hotspot_breakpoint.bedpe", row.names = F, col.names = F, sep = "\t", quote = F)
str(bedpe)

#write.csv(kpBedInfo %>% select(UpChr,Start1,End1,DownChr,Start2,End2,Name,Strand1,Strand2,Info),'~/data/sv/known_pairs_bedpe',row.names = F,quote = F)


#hotspots = read.csv(file = "~/hmf/resources/StructuralVariantHotspots.hg19.bedpe")

bedpe = hotspots %>%
  mutate(Start1 = Start1 - 1, Start2 = Start2 - 1, Score = 0) %>%
  select(Chr1, Start1, End1, Chr2, Start2, End2, Name, Score, Strand1, Strand2) %>%
  arrange(Chr1, Start1, End2, Chr2, Start2, End2)




