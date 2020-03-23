

#####
## Arriba vs StarFusion RNA fusions


## StarFusion calls
starRnaFusions = read.csv('~/data/sv/rna/LNX_RNA_DATA.csv')
dnaRnaStarFusionResults = read.csv('~/data/sv/rna/LINX_dna_rna_combined_data.csv')
View(dnaRnaStarFusionResults)

sampleIdMapping = read.csv('~/data/sv/sample_id_mapping.csv')
dnaRnaStarFusionResults = merge(dnaRnaStarFusionResults,sampleIdMapping,by='HmfId',all.x=T)


# good candidates for comparison with Arriba
View(dnaRnaStarFusionResults %>% filter(MatchType=='DNA & RNA'&KnownCategory=='Known'&RegionTypeUp=='Intronic'&RegionTypeDown=='Intronic'
                                        &ExonsSkippedUp==0&ExonsSkippedDown==0&ChainLinks==0) %>%
       select(SampleId,GeneNameUp,GeneNameDown,everything()))

View(dnaRnaStarFusionResults %>% filter(GeneNameUp=='AHR'))

View(dnaRnaStarFusionResults)

sampleSubset = c('CPCT02020378T','CPCT02020618T','CPCT02020723T','CPCT02040216T','CPCT02140007T','DRUP01010047T')

View(dnaRnaStarFusionResults %>% filter(SampleId %in% sampleSubset) %>% 
       filter(MatchType=='DNA & RNA'&(KnownCategory=='Known'|GeneNameUp=='AHR')&RegionTypeUp=='Intronic'&RegionTypeDown=='Intronic'&ExonsSkippedUp==0&ExonsSkippedDown==0&ChainLinks==0))

View(dnaRnaStarFusionResults %>% filter(SampleId %in% sampleSubset) %>% group_by(GeneNameUp,GeneNameDown) %>% count)

compKnownFusions = dnaRnaStarFusionResults %>% filter(SampleId %in% sampleSubset) %>% 
  filter(MatchType=='DNA & RNA'&KnownCategory=='Known'&RegionTypeUp=='Intronic'&RegionTypeDown=='Intronic'&ExonsSkippedUp==0&ExonsSkippedDown==0&ChainLinks==0) %>%
  select(SampleId,GeneNameUp,GeneNameDown,everything())
View(compKnownFusions)
write.csv(compKnownFusions,'~/data/sv/rna/arriba/sf_comp_known_fusions.csv',row.names = F, quote = F)






## Arriba fusions

arribaFusions = read.csv('~/data/sv/rna/arriba/CPCT02020378T.fusions.tsv',sep='\t')
arribaDisFusions = read.csv('~/data/sv/rna/arriba/CPCT02020378T.fusions.discarded.tsv',sep='\t')
View(arribaFusions)
nrow(arribaDisFusions)



# Matching with StarFusion









# Field value examples

colnames(arribaFusions)
View(arribaFusions %>% filter(X.gene1=='TMPRSS2'&gene2=='ERG'))
View(arribaDisFusions %>% filter(X.gene1=='TMPRSS2'&gene2=='ERG'))

View(arribaDisFusions %>% filter(closest_genomic_breakpoint2!='.'))

View(arribaFusions %>% group_by(strand1.gene.fusion.) %>% count)
View(arribaFusions %>% group_by(strand2.gene.fusion.) %>% count)

View(arribaFusions %>% group_by(site1,site2) %>% count)
View(arribaFusions %>% group_by(type) %>% count)
View(arribaFusions %>% group_by(direction1,direction2) %>% count)

View(arribaFusions %>% group_by(reading_frame) %>% count)




# DEBUG






