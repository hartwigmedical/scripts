library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)

wideGenes = read.csv('~/dev/bachelor/WIDE_genes.csv')
View(wideGenes)

ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)
ensemblTransData = read.csv('~/data/sv/ensembl_trans_data.csv')
View(ensemblTransData)

wideGeneData = merge(wideGenes,ensemblGeneData %>% select(GeneName,GeneId),by='GeneName',all.x=T)
View(wideGeneData)
View(wideGeneData %>% filter(is.na(GeneId)))

wideGeneData = wideGeneData %>% filter(!is.na(GeneId))

wideGeneData = merge(wideGeneData,ensemblTransData %>% filter(IsCanonical) %>% select(StableId,GeneId),by='GeneId',all.x=T)
write.csv(wideGeneData %>% select(GeneId,GeneName,TransId=StableId),'~/dev/bachelor/wide_gene_data.csv',row.names = F,quote = F)

# data for XML
write.csv(wideGeneData %>% select(GeneName,StableId),'~/dev/bachelor/bachelor_wide_genes.xml',row.names = F,quote = F)

# unmatched genes - FLT3, GSK3b, HPRT2, MSI, PTCH1, RAF-1, VEGFR, VHL 
