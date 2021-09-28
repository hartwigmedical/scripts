
#####
# Transcript comparisons - Ensembl 37 vs 38 vs HGNC
#####

# limit comparison to the genes in HGNC


ensTransData37 = read.csv('~/data/ensembl_hg37/ensembl_trans_exon_data.csv')
ensTransData37 = ensTransExonData37 %>% filter(ExonRank==1)
nrow(ensTransData37)
View(ensTransData37)
ensTransData37 = merge(ensTransData37,ensemblGeneData37 %>% select(GeneId,GeneName),by='GeneId',all.x=T)
ensCanTransData37 = ensTransData37 %>% filter(CanonicalTranscriptId==TransId)
nrow(ensCanTransData37) # 55736 as expected

View(ensTransData37 %>% filter(GeneName %in% c('BRCA1','CDKN2A','CHEK2','MITF','MET','MUTYH','NF1','PTCH1','RAD51D','WT1')))

ensCanTransData37 = merge(ensCanTransData37,ensemblGeneData37 %>% select(GeneId,GeneName,Chromosome),by='GeneId',all.x=T)

allGenes37 = read.csv('~/data/ensembl_db/ensembl_37_89/all_genes.37.tsv',sep='\t')
allGenes37 = read.csv('~/data/ensembl_db/ensembl_37_104/all_genes.37.tsv',sep='\t')

View(allGenes37)
View(allGenes37 %>% filter(!(GeneId %in% ensemblGeneData37$GeneId)))
View(allGenes37 %>% filter(Chromosome=='MT'))

View(allGenes37 %>% filter(GeneId %in% ensemblGeneData37$GeneId & !(TranscriptId %in% ensCanTransData37$TransName)))



ensTransExonData38 = read.csv('~/data/ensembl_db/ensembl_38_104/ensembl_trans_exon_data.csv')
nrow(ensTransExonData38) # 
ensTransData38 = ensTransExonData38 %>% filter(ExonRank==1)
nrow(ensTransData38) # 205K
View(ensTransData38 %>% filter(TransName %in% c('ENST00000508376','ENST00000357654','ENST00000304494','ENST00000328354','ENST00000397752',
                                                'ENST00000394351','ENST00000450313','ENST00000356175','ENST00000437951','ENST00000345365','ENST00000452863')))

View(ensTransData38 %>% filter(GeneId=='ENSG00000183765'))
ensCanTransData38 = ensTransData38 %>% filter(CanonicalTranscriptId==TransId)
nrow(ensCanTransData38) # 39357


hgncGeneDataEnsembl = hgncGeneDataEnsembl %>% mutate(TransName=ifelse(mane_select!='',substr(mane_select,1,15),''))
View(hgncGeneDataEnsembl)
nrow(hgncGeneDataEnsembl %>% filter(TransName=='')) # 21892 have no transcript specified
nrow(hgncGeneDataEnsembl %>% filter(TransName!='')) # 17755 have transcript specified

# but all the ones which do have the same canonical transcript as marked in Ensembl
hgncGeneDataEns = merge(hgncGeneDataEnsembl %>% 
                          filter(ensembl_gene_id %in% ensemblGeneData38$GeneId) %>% 
                          filter(TransName!='') %>% 
                          select(GeneId=ensembl_gene_id,GeneName=symbol,HgncId=hgnc_id,TransName),
                        ensTransData38 %>% mutate(IsCanonical=CanonicalTranscriptId==TransId) %>% select(TransName,IsCanonical,BioType),by='TransName',all.x=T)

View(hgncGeneDataEns)

# ENST00000373997.8|NM_014576.4


# compare canonical transcripts between 37 and 38
ensCanTransData38Final = ensCanTransData38 %>% filter(GeneId %in% geneData38Final$GeneId)
nrow(ensCanTransData38Final) # 39225
View(ensCanTransData38Final) 

nrow(ensTransData37) # 194K
View(ensTransData37 %>% filter(GeneId=='ENSG00000004139'))

ensTransData37Final = ensTransData37 %>% filter(GeneId %in% ensemblGeneData37Final$GeneId)
ensTransData37Final = ensTransData37Final %>% mutate(Canonical37=(CanonicalTranscriptId==TransId))
nrow(ensTransData37Final) # 172K
View(ensTransData37Final) 

ensTransData37Final = merge(ensTransData37Final,ensCanTransData38Final %>% mutate(Canonical38=T) %>% select(TransName,Canonical38),by='TransName',all.x=T)
View(ensTransData37Final)

ensCanTransData37Final = ensTransData37Final %>% filter(Canonical37|!is.na(Canonical38))
nrow(ensCanTransData37Final)
View(ensCanTransData37Final)

ens37vs38Changes = read.csv('~/logs/ens37_38_trans_mapping.csv')
ens37vs38Changes = ens37vs38Changes %>% select(-Time)
View(ens37vs38Changes)
View(ens37vs38Changes %>% filter(MatchType!='NoRefByGeneId') %>% group_by(MatchType) %>% count)
ens37vs38Changes = merge(ens37vs38Changes,ensemblGeneData37 %>% select(GeneId,GeneName),by='GeneId',all.x=T)

ensNmMappings = read.csv('~/data/ensembl_db/test/ensembl_trans_nm_query.csv',sep='\t')
View(ensNmMappings)

ens37vs38Changes = merge(ens37vs38Changes,ensNmMappings %>% select(SelectedTrans=stable_id,NmId=display_label),by='SelectedTrans',all.x=T)
ens37vs38Changes = merge(ens37vs38Changes,ensNmMappings %>% select(CanonicalTrans=stable_id,NmIdCanonical=display_label),by='CanonicalTrans',all.x=T)
View(ens37vs38Changes)

View(ens37vs38Changes %>% group_by(CanonicalTrans) %>% count)

ens37vs38Changes = ens37vs38Changes %>% mutate(NmIdSelected=ifelse(is.na(NmId),'',as.character(NmId)),
                                               NmIdCanonical=ifelse(is.na(NmIdCanonical),'',as.character(NmIdCanonical))) %>% select(-NmId)


write.csv(ens37vs38Changes,'~/data/ensembl_db/ens37_38_trans_mapping.csv',row.names = F,quote = F)

write.csv(ens37vs38Changes %>% group_by(SelectedTrans,CanonicalTrans,Ref38Trans,MatchType,CodingBasesSelected,CodingBasesCanonical38,CodingBasesCanonical37) %>%
            summarise(NmIdSelected=first(NmIdSelected),NmIdCanonical=first(NmIdCanonical)),'~/data/ensembl_db/ens37_38_trans_mapping.csv',row.names = F,quote = F)


# focus on the differences
ensCanTransData37Diffs = ensCanTransData37Final %>% filter(is.na(Canonical38)|Canonical38!=Canonical37)
View(ensCanTransData37Diffs)
View(ensCanTransData37Diffs %>% group_by(GeneId) %>% count)

View(rbind(ensTransExonData38 %>% filter(TransName=='ENST00000585482'),
           ensTransExonData37 %>% filter(TransName=='ENST00000457710')))

View(rbind(ensTransExonData37 %>% filter(GeneId=='ENSG00000141622') %>% mutate(Version='V37') %>% select(Version,everything()),
           ensTransExonData38 %>% filter(GeneId=='ENSG00000141622') %>% mutate(Version='V38') %>% select(Version,everything())))

# 10:41:15.560 [main] [DEBUG] geneId(ENSG00000204130) choosing trans(ENST00000602465 matching ref name as canonical over designated canonical(ENST00000388768)
View(rbind(ensTransExonData37 %>% filter(TransName %in% c('ENST00000602465','ENST00000388768')) %>% mutate(Version='V37',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything()),
           ensTransExonData38 %>% filter(TransName=='ENST00000602465') %>% mutate(Version='V38',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything())))

# 10:41:07.501 [main] [WARN ] geneId(ENSG00000168386) non-canonical trans(ENST00000477258) name match to ref has different definition, vs designated canonical(ENST00000354552)
View(rbind(ensTransExonData37 %>% filter(TransName %in% c('ENST00000477258','ENST00000354552')) %>% mutate(Version='V37',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything()),
           ensTransExonData38 %>% filter(TransName=='ENST00000477258') %>% mutate(Version='V38',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything())))

# 10:57:05.916 [main] [WARN ] geneId(ENSG00000204387) canonical trans(ENST00000375640) name match but has different definitions
View(rbind(ensTransExonData37 %>% filter(GeneId=='ENSG00000204387') %>% mutate(Version='V37',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything()),
           ensTransExonData38 %>% filter(TransName=='ENST00000375640') %>% mutate(Version='V38',ExonLength=ExonEnd-ExonStart) %>% 
             select(Version,GeneId,TransId,CanonicalTranscriptId,TransName,ExonRank,ExonLength,ExonPhase,ExonEndPhase,everything())))

View(rbind(ensTransExonData37 %>% filter(TransName %in% c('ENST00000372067','ENST00000372055')) %>% mutate(Version='V37') %>% select(Version,everything()),
           ensTransExonData38 %>% filter(TransName %in% c('ENST00000672860')) %>% mutate(Version='V38') %>% select(Version,everything())))

# trans(ENST00000372067 matches ref=ENST00000672860) as canonical over designated canonical(ENST00000372055) 

tmpTransData = read.csv('~/data/ensembl_db/test/ensembl_trans_exon_data.csv')
View(tmpTransData)




