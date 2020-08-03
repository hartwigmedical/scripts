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





######
## Clinvar Filters
cvFilters = read.csv('~/logs/bachelor_clinvar_filters.csv')
nrow(cvFilters)
View(cvFilters)
View(cvFilters %>% group_by(CodingEffect) %>% count)
View(cvFilters %>% filter(CodingEffect==''))
View(cvFilters %>% filter(CodingEffect=='NONE'|CodingEffect=='SYNONYMOUS'))
View(cvFilters %>% group_by(ClinvarSignificance) %>% count)

View(cvFilters %>% filter(grepl('Conflicting',ClinvarSignificance)) %>% 
       filter(!grepl('Pathogenic',ClinvarSigInfo)&!grepl('Likely_pathogenic',ClinvarSigInfo)&!grepl('Benign',ClinvarSigInfo)&!grepl('Likely_benign',ClinvarSigInfo)))

sampleRecords = read.csv('~/logs/DRUP01090022T.bachelor.germline_variant.tsv',sep='\t')
View(sampleRecords)
View(sampleRecords %>% group_by(filter,reported, pathogenic,clinvarSignificance,effects) %>% count)

