
# HGNC which will become our reference gene set for HMF
hgncGenes = read.csv('~/data/ensembl_db/hgnc_ensembl_gene_data.csv')
nrow(hgncGenes) # 39647
View(hgncGenes)

# full data
hgncGeneDataEnsembl = hgncGeneData %>% filter(ensembl_gene_id!='')
View(hgncGeneDataEnsembl)


#####
# Ensembl GRCh38 vs HGNC
#####

# Ensembl gene query restricted by join to external DB = 1100 (HGNC)
geneTest38 = read.csv('~/data/ensembl_db/test/ensembl_gene_hgnc_38.csv',sep='\t')
nrow(geneTest38) # 39730
View(geneTest38)
View(geneTest38 %>% filter(as.character(GeneName)!=as.character(display_label))) # only 5 are different

geneData38 = merge(geneTest38 %>% select(HgncId=dbprimary_acc,GeneIdEns=GeneId,GeneNameEns=GeneName),
                   hgncGeneData %>% select(HgncId=hgnc_id,GeneNameHgnc=symbol,GeneIdHgnc=ensembl_gene_id),by='HgncId',all = T)
View(geneData38)

geneData38 = geneData38 %>% mutate(MapType=ifelse(is.na(GeneIdEns),'HGNC Only',ifelse(is.na(GeneIdHgnc),'Ensembl Only','Both')),
                                   GeneIdMap=ifelse(MapType!='Both','N/A',ifelse(as.character(GeneIdEns)==as.character(GeneIdHgnc),'Match','Diff')),
                                   GeneNameMap=ifelse(MapType!='Both','N/A',ifelse(as.character(GeneNameEns)==as.character(GeneNameHgnc),'Match','Diff')))

View(geneData38 %>% group_by(MapType,GeneIdMap,GeneNameMap) %>% count)

# FINAL: take the set of genes with HGNC IDs, using GeneId from Ensembl and GeneName (Symbol) from HGNC
geneData38Final = geneData38 %>% filter(MapType=='Both') %>% select(HgncId,GeneId=GeneIdEns,GeneName=GeneNameHgnc)
nrow(geneData38Final) # 39718
View(geneData38Final)
View(geneData38Final %>% group_by(GeneId) %>% count)
nrow(geneData38Final %>% group_by(GeneId) %>% count) # 39243

write.csv(geneData38Final,'~/data/ensembl_db/ensembl_hgnc_gene_set.csv',row.names = F,quote = F)

ensemblGeneData38New = read.csv('~/data/ensembl_db/ensembl_38_104/ensembl_gene_data.csv')
nrow(ensemblGeneData38New) # 39243

View(geneData38Final %>% filter(!(GeneId %in% ensemblGeneData38New$GeneId)))


# unmapped: 12 in Ensembl but not found in HGNC
View(geneData38 %>% filter(is.na(GeneIdEns)|is.na(GeneIdHgnc)))

# diff Ensembl IDs
View(geneData38 %>% filter(!is.na(GeneId.x)!is.na(GeneId.y)!as.character(GeneId.x)!=as.character(GeneId.y)|as.character(GeneName.x)!=as.character(GeneName.y)))

# diff names
View(geneData38 %>% filter(is.na(GeneId.x)|is.na(GeneId.y)|as.character(GeneId.x)!=as.character(GeneId.y)|as.character(GeneName.x)!=as.character(GeneName.y)))


# Ensembl 37 mapping to the HGNC / Ensembl 38 set of genes
nrow(ensemblGeneData37 %>% filter(GeneId %in% geneData38Final$GeneId)) # 37058
nrow(ensemblGeneData37 %>% filter(GeneId %in% geneData38Final$GeneId | GeneName %in% geneData38Final$GeneName)) # 37832


# FINAL: take the set in HGNC from v38
# v37 genes must match on 


nrow(ensemblGeneData37 %>% filter(GeneId %in% ensemblGeneData38$GeneId)) # 37186

nrow(ensemblGeneData37 %>% filter(GeneId %in% ensemblGeneData38$GeneId)) # 37186


nrow(ensemblGeneData37 %>% filter(GeneId %in% ensemblGeneData38$GeneId | GeneName %in% ensemblGeneData38$GeneName)) # 38366
View(ensemblGeneData37 %>% filter(!(GeneId %in% ensemblGeneData38$GeneId) & GeneName %in% ensemblGeneData38$GeneName)) # 

# FINAL: take the set in HGNC
nrow(ensemblGeneData37 %>% filter(GeneId %in% hgncGeneData$ensembl_gene_id | GeneName %in% hgncGeneData$symbol)) # 38178
nrow(ensemblGeneData37 %>% filter(GeneId %in% hgncGeneData$ensembl_gene_id)) # 37404
nrow(ensemblGeneData37 %>% filter(GeneName %in% hgncGeneData$symbol)) # 31403


ensemblGeneData37Final = ensemblGeneData37 %>% filter(GeneId %in% hgncGeneData$ensembl_gene_id | GeneName %in% hgncGeneData$symbol)
nrow(ensemblGeneData37Final) # 38178

# which are now in 38 and HGNC but not 37?
View(geneData38Final %>% filter(!(GeneId %in% ensemblGeneData37$GeneId) & !(GeneName %in% ensemblGeneData37$GeneName)))






# Ensembl Gene Data GRCH37 v89
ensemblGeneData37 = read.csv('~/data/ensembl_hg37/ensembl_gene_data.csv')

# mis-matches
nrow(hgncGenes %>% filter(!(GeneId %in% ensemblGeneData37$GeneId))) # 2241 not matched on GeneId
nrow(hgncGenes %>% filter(!(GeneId %in% ensemblGeneData37$GeneId) & !(GeneName %in% ensemblGeneData37$GeneName))) # 1486 not matched on GeneId or GeneName

nrow(hgncGenes %>% filter(!(GeneId %in% ensemblGeneData38$GeneId) & !(GeneName %in% ensemblGeneData38$GeneName))) # 411 not matched on GeneId or GeneName


# name differences
hgnvVsEnsById = merge(hgncGeneDataEnsembl,ensemblGeneData37,by.x='ensembl_gene_id',by='GeneId',all.x=T)
View(hgnvVsEnsById %>% filter(is.na(GeneName)))

# 7024 have different names
View(hgnvVsEnsById %>% filter(!is.na(GeneName) & as.character(GeneName)!=as.character(symbol)) %>% 
       select(ensembl_gene_id,symbol,name,GeneName,Chromosome,GeneStart,GeneEnd,location,alias_symbol,prev_symbol))

# 1257 match on a previous name, and many of the rest match amongst the collection of previous names
nrow(hgnvVsEnsById %>% filter(!is.na(GeneName) & as.character(GeneName)!=as.character(symbol)) %>% 
       filter(as.character(GeneName)==as.character(prev_symbol) | grepl(as.character(GeneName),prev_symbol)))

View(hgncGenes %>% filter(!(GeneId %in% ensemblGeneData37$GeneId) & !(GeneName %in% ensemblGeneData37$GeneName)))

genePanel37 = read.csv('~/data/ensembl_db/ensembl_37_89/all_genes.37.tsv',sep = '\t')
nrow(genePanel37) # 25965
View(genePanel37) # 25965

View(genePanel37 %>% filter(!(GeneId %in% ensemblGeneData37$GeneId))) # 19 without a match
View(genePanel37 %>% filter(!(GeneName %in% ensemblGeneData37$GeneName))) # 34 without a match
View(genePanel37 %>% filter(!(GeneName %in% hgncGenes$GeneName))) # 4K without a match
View(genePanel37 %>% filter(!(GeneId %in% hgncGenes$GeneId))) # 2424 without a match

# Driver Genes
driverGenes37 = read.csv('~/data/DriverGenePanel.hg19.tsv',sep='\t')
driverGenes37 = merge(driverGenes37,ensemblGeneData37 %>% select(gene=GeneName,GeneId,Chromosome,GeneStart,GeneEnd,Synonyms),by='gene',all.x=T)
nrow(driverGenes37) # 448
View(driverGenes37) # 448
nrow(driverGenes37 %>% filter(is.na(GeneId))) # all are in Ensembl

# driver genes in HGNC

# 5 genes are not matched on geneId, but 4 of these can be matched by GeneName
nrow(driverGenes37 %>% filter(!(GeneId %in% hgncGenes$GeneId))) 
nrow(driverGenes37 %>% filter(!(GeneId %in% hgncGenes$GeneId) & gene %in% hgncGenes$GeneName)) # 4 of these

# 10 genes have changed name
View(driverGenes37 %>% filter(!(gene %in% hgncGenes$GeneName))) 
nrow(driverGenes37 %>% filter(!(gene %in% hgncGenes$GeneName))) # 

driverDiffHgncName = merge(driverGenes37,hgncGenes,by='GeneId',all.x=T) %>% filter(!is.na(GeneName)&as.character(gene)!=as.character(GeneName))
View(driverDiffHgncName %>% select(GeneId,DriverGeneName=gene,HgncGeneName=GeneName,everything()))

# only HIST1H3B/ENSG00000124693 is not matched with HGNC by either GeneId or GeneName, but is matched on prev symbol with HGNC:4776 GeneName=H3C2 and GeneId = ENSG00000286522
View(driverGenes37 %>% filter(!(GeneId %in% hgncGenes$GeneId) & !(gene %in% hgncGenes$GeneName))) # 



View(driverGenes37 %>% filter(!(gene %in% hgncGenes$GeneName) & gene %in% hgncGeneDataEnsembl$name)) # 

knownFusionData = read.csv('~/data/fusion_ref/known_fusion_data.37_v3.csv')
View(knownFusionData)


fusionGeneData = rbind(knownFusionData %>% filter(FiveGene!=''&!(FiveGene %in% c('IGH','IGL','IGK'))) %>% mutate(Stream='Five') %>% select(GeneName=FiveGene,Stream,KnownType=Type),
                       knownFusionData %>% filter(ThreeGene!='') %>% mutate(Stream='Three') %>% select(GeneName=ThreeGene,Stream,KnownType=Type))

fusionGeneData = merge(fusionGeneData %>% group_by(GeneName,Stream) %>% summarise(KnownType=first(KnownType)) %>% ungroup(),ensemblGeneData37,by='GeneName',all.x=T)
View(fusionGeneData)
View(fusionGeneData %>% filter(!(GeneName %in% hgncGeneDataEnsembl$symbol)))
View(fusionGeneData %>% filter(!(GeneId %in% hgncGeneDataEnsembl$ensembl_gene_id)))

# 8 genes have different GeneIds but names match and map to Ensembl 38
fusionSameNameDiffId = fusionGeneData %>% filter(!(GeneId %in% hgncGeneDataEnsembl$ensembl_gene_id) & GeneName %in% hgncGeneDataEnsembl$symbol)
View(fusionSameNameDiffId)
fusionSameNameDiffId = merge(fusionSameNameDiffId,ensemblGeneData38 %>% select(GeneName,GeneId38=GeneId),by='GeneName',all.x=T)
View(fusionSameNameDiffId)
View(fusionSameNameDiffId %>% filter(!(GeneId38 %in% hgncGeneDataEnsembl$ensembl_gene_id)))

fusionDiffNameSameId = fusionGeneData %>% filter(GeneId %in% hgncGeneDataEnsembl$ensembl_gene_id & !(GeneName %in% hgncGeneDataEnsembl$symbol))
View(fusionDiffNameSameId)
fusionDiffNameSameId = merge(fusionDiffNameSameId,ensemblGeneData38 %>% select(GeneName38=GeneName,GeneId),by='GeneId',all.x=T)
View(fusionDiffNameSameId)

# 7 have new names: DUX4 -> DUX4L1, 
# only one gene cannot be mapped by GeneId or GeneName: RP11-356O9.1

missingFusionGenes = rbind(knownFusionData %>% filter((FiveGene!='' & !(FiveGene %in% c('IGH','IGL','IGK')) & !(FiveGene %in% hgncGenes$GeneName))) %>% mutate(Type='Five'),
                           knownFusionData %>% filter((ThreeGene!='' & !(ThreeGene %in% hgncGenes$GeneName))) %>% mutate(Type='Three'))
                  
View(missingFusionGenes)
         
View(knownFusionData %>% filter((FiveGene!='' & !(FiveGene %in% c('IGH','IGL','IGK')) & !(FiveGene %in% hgncGenes$GeneName))
                                |(ThreeGene!='' & !(ThreeGene %in% hgncGenes$GeneName))))








# Gene Mapping - v37 and v38

genes37 = read.csv('~/data/ensembl_db/ensembl_37_89/all_genes.37.tsv',sep = '\t')
allGenes38 = read.csv('~/data/ensembl_db/ensembl_38_89/all_genes.38.tsv',sep = '\t')
nrow(genes37) # 25965
nrow(allGenes38) # 25446

matchedGenes = merge(genes37 %>% select(GeneId,GeneName,TransId37=TranscriptId),
                     allGenes38 %>% select(GeneId,GeneName,TransId38=TranscriptId),by=c('GeneId','GeneName'),all=F)

nrow(matchedGenes) # 20868 are unch

genesDiff37 = genes37 %>% filter(!(GeneId %in% matchedGenes$GeneId))
genesDiff38 = allGenes38 %>% filter(!(GeneId %in% matchedGenes$GeneId))

mergedGenesById = merge(genesDiff37 %>% select(GeneId,GeneName37=GeneName,TransId37=TranscriptId),
                        genesDiff38 %>% select(GeneId,GeneName38=GeneName,TransId38=TranscriptId),by='GeneId',all=F)

nrow(mergedGenesById) # 2119
nrow(mergedGenesById %>% filter(as.character(TransId37)==as.character(TransId38))) # 1921 keep the same canonical transcript
View(mergedGenesById)

mergedGenesByName = merge(genesDiff37 %>% filter(!(GeneId %in% mergedGenesById$GeneId)) %>% select(GeneId37=GeneId,GeneName,TransId37=TranscriptId),
                          genesDiff38 %>% filter(!(GeneId %in% mergedGenesById$GeneId)) %>% select(GeneId38=GeneId,GeneName,TransId38=TranscriptId),by='GeneName',all=F)

nrow(mergedGenesByName) # 488

unmatchedGenes37 = genesDiff37 %>% filter(!(GeneId %in% mergedGenesById$GeneId) & !(GeneName %in% mergedGenesByName$GeneName))
unmatchedGenes38 = genesDiff38 %>% filter(!(GeneId %in% mergedGenesById$GeneId) & !(GeneName %in% mergedGenesByName$GeneName))
nrow(unmatchedGenes37)
nrow(unmatchedGenes38)

driverGenes37 = read.csv('~/data/DriverGenePanel.hg19.tsv',sep='\t')
nrow(driverGenes37)
driverGenes38 = read.csv('~/data/DriverGenePanel.38.tsv',sep='\t')
nrow(driverGenes38)

View(driverGenes37 %>% filter(gene %in% unmatchedGenes37$GeneName))
View(driverGenes38 %>% filter(gene %in% unmatchedGenes38$GeneName))

write.csv(unmatchedGenes37 %>% select(GeneId),'~/data/ensembl_db/unmatched_gene_ids.37.csv',row.names = F,quote = F)
write.csv(unmatchedGenes38 %>% select(GeneId),'~/data/ensembl_db/unmatched_gene_ids.38.csv',row.names = F,quote = F)

# BED of 37 genes for lift-over
write.table(unmatchedGenes37 %>% select(Chromosome,GeneStart,GeneEnd),'~/data/ensembl_db/unmatched_genes.37.bed',row.names = F,quote = F,sep='\t')

unmatchedLiftOverTo38 = read.csv('~/data/ensembl_db/unmatched_genes.37_lifted_over.tsv',sep='\t')
View(unmatchedLiftOverTo38)

View(genes38)
mergedLiftOver38 = merge(unmatchedGenes38,unmatchedLiftOverTo38 %>% select(Chromosome,GeneStart,GeneEnd),by=c('Chromosome','GeneStart','GeneEnd'),all=F)
View(mergedLiftOver38)

library(GenomicRanges)

# comparison using HGNC
extraGene37 = read.csv('~/data/ensembl_db/test/extra_gene_info.37.tsv',sep = '\t')
extraGene38 = read.csv('~/data/ensembl_db/test/extra_gene_info.38.tsv',sep = '\t')

View(extraGene37)
View(extraGene38)

extraGene38 = extraGene38 %>% mutate(Hgnc=stri_replace_all_fixed(entrezId_1100,'HGNC:',''))

mergedGenesByHgnc = merge(extraGene37 %>% filter(entrezId_1100!='') %>% select(GeneId37=gene_id,GeneName37=gene_name,entrezId_1100),
                          extraGene38 %>% filter(entrezId_1100!='') %>% select(GeneId38=gene_id,GeneName38=gene_name,,entrezId_1100=Hgnc),by='entrezId_1100',all=F)

View(mergedGenesByHgnc)

mergedGenesByEntrez = merge(extraGene37 %>% filter(entrezId_1300!='') %>% select(GeneId37=gene_id,GeneName37=gene_name,entrezId_1300),
                            extraGene38 %>% filter(entrezId_1300!='') %>% select(GeneId38=gene_id,GeneName38=gene_name,,entrezId_1300),by='entrezId_1300',all=F)

View(mergedGenesByEntrez)


# full Ensembl cache - 37 v89 vs 38 v104
ensemblGeneData = read.csv('~/data/ensembl_hg37/ensembl_gene_data.csv')
nrow(ensemblGeneData) # 55736
ensemblGeneData38 = read.csv('~/data/ensembl_db/ensembl_38_104/remote/ensembl_gene_data.csv')
nrow(ensemblGeneData38) # 39357

ensMatchedGenes = merge(ensemblGeneData %>% select(GeneId,GeneName),ensemblGeneData38 %>% select(GeneId,GeneName),by=c('GeneId','GeneName'),all=F)

nrow(ensMatchedGenes) # 30465 are unch

ensGenesDiff37 = ensemblGeneData %>% filter(!(GeneId %in% ensMatchedGenes$GeneId))
ensGenesDiff38 = ensemblGeneData38 %>% filter(!(GeneId %in% ensMatchedGenes$GeneId))

ensMergedGenesById = merge(ensGenesDiff37 %>% select(GeneId,GeneName37=GeneName),
                           ensGenesDiff38 %>% select(GeneId,GeneName38=GeneName),by='GeneId',all=F)

nrow(ensMergedGenesById) # 6669
View(ensMergedGenesById)

ensMergedGenesByName = merge(ensGenesDiff37 %>% filter(!(GeneId %in% ensMergedGenesById$GeneId)) %>% select(GeneId37=GeneId,GeneName),
                             ensGenesDiff38 %>% filter(!(GeneId %in% ensMergedGenesById$GeneId)) %>% select(GeneId38=GeneId,GeneName),by='GeneName',all=F)

nrow(ensMergedGenesByName) # 708

ensUnmatchedGenes37 = ensGenesDiff37 %>% filter(!(GeneId %in% ensMergedGenesById$GeneId) & !(GeneName %in% ensMergedGenesByName$GeneName))
ensUnmatchedGenes38 = ensGenesDiff38 %>% filter(!(GeneId %in% ensMergedGenesById$GeneId) & !(GeneName %in% ensMergedGenesByName$GeneName))
nrow(ensUnmatchedGenes37) # 17894
nrow(ensUnmatchedGenes38) # 1516

# BED of 37 genes for lift-over
write.table(ensUnmatchedGenes37 %>% mutate(Chromosome=paste0('chr',Chromosome)) %>% select(Chromosome,GeneStart,GeneEnd,GeneId),
            '~/data/ensembl_db/ens_unmatched_genes.37.bed',row.names = F,quote = F,sep='\t')


View(ensUnmatchedGenes37)
write.csv(ensUnmatchedGenes37,'~/data/ensembl_db/unmatched_ensembl_gene_data.37.csv',row.names = F,quote = F)

ensMapping = read.csv('~/data/ensembl_db/ensembl_gene_mapping_89_vs_104.csv')
nrow(ensMapping)
View(ensMapping)
View(ensMapping %>% group_by(MappingType) %>% count)
View(ensMapping %>% group_by(MappingType,OtherInfo) %>% count)

View(ensMapping %>% filter(MappingType!='NONE') %>% filter(as.character(GeneName37)!=as.character(GeneName38)))

View(ensMapping %>% filter(MappingType!='NONE') %>% filter(as.character(GeneName37)!=as.character(GeneName38)) %>%
       group_by(MappingType) %>% count)


View(ensMapping %>% filter(MappingType!='NONE') %>% filter(as.character(GeneName37)!=as.character(GeneName38)) %>% filter(GeneName37 %in% genes37$GeneName))


View(ensMapping %>% filter(MappingType!='NONE') %>% filter(as.character(GeneName37)!=as.character(GeneName38)) %>% filter(GeneName37 %in% driverGenes37$gene))


# fusion genes
View(ensMapping %>% filter(GeneName37 %in% c('MLLT4','RP11-356O9.1','PPAP2B','KIAA1524','LHFP','MKL2','WHSC1L1',
                                             'CXorf67','MGEA5','CARS','SEPT14','FGFR1OP','SEPT2','SEPT5','SEPT6','SEPT9','WHSC1')))

# driver genes
driverGenes37 = read.csv('~/data/DriverGenePanel.hg19.tsv',sep='\t')
nrow(driverGenes37) # 448
driverGenes38 = read.csv('~/data/DriverGenePanel.38.tsv',sep='\t')
nrow(driverGenes38) # 460

View(ensMapping %>% filter(GeneName37 %in% driverGenes37$gene))
View(ensMapping %>% filter(GeneName37 %in% driverGenes37$gene) %>% group_by(MappingType) %>% count)
View(ensMapping %>% filter(GeneName38 %in% driverGenes38$gene) %>% group_by(MappingType) %>% count)
View(driverGenes38 %>% filter(!(gene %in% ensMapping$GeneName38)))



# gene panel
View(ensMapping %>% filter(GeneName37 %in% genes37$GeneName) %>% group_by(MappingType) %>% count)


# mapping between Ensembl version 89
ensMapping89 = read.csv('~/data/ensembl_db/ensembl_38_89/ensembl_gene_mapping.csv')
nrow(ensMapping89)
View(ensMapping89 %>% group_by(MappingType) %>% count)
View(ensMapping89 %>% group_by(MappingType,OtherInfo) %>% count)

# driver genes
View(ensMapping89 %>% filter(GeneName37 %in% driverGenes37$gene))
View(ensMapping89 %>% filter(GeneName37 %in% driverGenes37$gene) %>% group_by(MappingType) %>% count)

# gene panel
View(ensMapping89 %>% filter(GeneName37 %in% genes37$GeneName) %>% group_by(MappingType) %>% count)


#####
## Ensembl transcript mapping v89 vs v104
ensTransMapping = read.csv('~/data/ensembl_db/ensembl_transcript_mapping.csv')
View(ensTransMapping)
View(ensTransMapping %>% group_by(MatchInfo) %>% count)

# changed trans
View(ensTransMapping %>% filter(MatchInfo=='MATCHED') %>% group_by(ContainsCanonical,ChangedCanonical) %>% count)

unmatchedTrans = ensTransMapping %>% filter(MatchInfo=='NO_MATCH')

allGenes38 = read.csv('~/data/ensembl_db/ensembl_38_89/all_genes.38.tsv',sep = '\t')
View(allGenes38 %>% filter(GeneName %in% unmatchedTrans$GeneName))

View(driverGenes38 %>% filter(gene %in% unmatchedTrans$GeneName))



#####
## Ensembl v102 vs v104

ensemblGeneData38v102 = read.csv('~/data/ensembl_hg38/ensembl_gene_data.csv') # based on 102
ensemblGeneData38v104 = read.csv('~/data/ensembl_db/ensembl_38_104/cache/ensembl_gene_data.csv')
View(ensemblGeneDataHg38)

View(ensemblGeneData38v102 %>% filter(!(GeneId %in% ensemblGeneData38v104$GeneId))) # about 20K dropped from 102 to 104
View(ensemblGeneData38v104 %>% filter(!(GeneId %in% ensemblGeneData38v102$GeneId))) # only 33 'new' genes

View(ensemblGeneData38v104 %>% filter(!(GeneId %in% ensemblGeneData38v102$GeneId) 
                                      & !(GeneName %in% ensemblGeneData38v102$GeneName))) # only 2 genes if GeneName used

View(ensemblGeneData38v104 %>% filter(!(GeneId %in% ensemblGeneData38v102$GeneId) 
                                      & GeneName %in% ensemblGeneData38v102$GeneName)) # 31 changed GeneId







ensemblGeneData37_104 = read.csv('~/data/ensembl_db/ensembl_37_104/ensembl_gene_data.csv')
nrow(ensemblGeneData37_104)

ensembl37Cmp = merge(ensemblGeneData37,ensemblGeneData37_104,by='GeneId',all=T)
nrow(ensembl37Cmp)
View(ensembl37Cmp %>% filter(as.character(GeneName.x)!=as.character(GeneName.y)))
View(ensembl37Cmp %>% filter(GeneStart.x!=GeneStart.y))
rm(ensembl37Cmp)
rm(ensemblGeneData37_104)


# gene name duplicates in Ensembl 37 104 - all explained by the link to external databases - no duplicates by name with different coords

geneData37_104 = read.csv('~/data/ensembl_db/test/ensembl_gene_query.csv',sep='\t')
nrow(geneData37_104)
View(geneData37_104)  
View(geneData37_104 %>% group_by(GeneId) %>% count %>% filter(n>1))  
View(geneData37_104 %>% group_by(GeneName,GeneId) %>% count %>% filter(n>1) %>% group_by(GeneName) %>% count %>% filter(n>1))  


