######################
## Gene Panel Creation

# GRCh38

# HGNC gene definitions
hgncGeneData = read.csv('~/data/ensembl_db/hgnc_complete_set_20210705.csv',sep='\t')
nrow(hgncGeneData) # 42630
View(hgncGeneData)
nrow(hgncGeneData %>% group_by(symbol) %>% count %>% filter(n>1)) # no duplicates by symbol
nrow(hgncGeneData %>% group_by(hgnc_id) %>% count %>% filter(n>1)) # no duplicates by hgnc ID

# only use those with Ensembl GeneIds
# select and rename columns for consistency
hgncGenes = hgncGeneData %>% filter(ensembl_gene_id!='') %>% select(HgncId=hgnc_id,GeneName=symbol,GeneId=ensembl_gene_id)
View(hgncGenes)
nrow(hgncGenes) # 39647

write.csv(hgncGenes,'~/data/ensembl_db/hgnc_ensembl_gene_data.csv',row.names = F,quote = F)

hgncGenes = read.csv('~/data/ensembl_db/hgnc_ensembl_gene_data.csv')


#####
# Ensembl data cache
#####

# first run GeneUtils for v38 to generate the reference gene list - intersection of Ensembl v38 (v104) with HGNC

# Ensembl gene data cache
ensemblGeneData38 = read.csv('~/data/ensembl_db/ensembl_38_104/ensembl_gene_data.csv')
nrow(ensemblGeneData38) # 39225
View(ensemblGeneData38) # 39243

# then generate v37 using v38 as the reference set
ensemblGeneData37 = read.csv('~/data/ensembl_db/ensembl_37_104/ensembl_gene_data.csv')
nrow(ensemblGeneData37) # 37873
View(ensemblGeneData37)


# no duplicates by geneName in either set
nrow(ensemblGeneData38 %>% group_by(GeneName) %>% count %>% filter(n>1)) # none
nrow(ensemblGeneData37 %>% group_by(GeneName) %>% count %>% filter(n>1)) # none

# compare to existing

# Ensembl Gene Data GRCH37 v89
ensemblGeneData37Old = read.csv('~/data/ensembl_hg37/ensembl_gene_data.csv')

# Ensembl Gene Data GRCH38 v104
ensemblGeneData38Old = read.csv('~/data/ensembl_db/ensembl_38_89/ensembl_gene_data.csv')
View(ensemblGeneData38)

# compare to each other
ensemblGeneDataCmp = merge(ensemblGeneData37 %>% select(GeneId37=GeneId,GeneName,Synonyms37=Synonyms),ensemblGeneData38 %>% select(GeneId38=GeneId,GeneName,Synonyms38=Synonyms),
                           by='GeneName',all=T)

ens37to38MapById = merge(ensemblGeneData37 %>% select(GeneId,GeneName37=GeneName,Chr37=Chromosome,KarBand37=KaryotypeBand),
                         ensemblGeneData38 %>% select(GeneId,GeneName38=GeneName,Chr38=Chromosome,KarBand38=KaryotypeBand),by='GeneId',all.x=T)

View(ens37to38MapById)
nrow(ens37to38MapById %>% filter(!is.na(GeneName38)&as.character(GeneName37)==as.character(GeneName38))) # 37182
View(ens37to38MapById %>% filter(!is.na(GeneName38)&as.character(GeneName37)!=as.character(GeneName38))) # 4 genes, all in same location
View(ens37to38MapById %>% filter(!is.na(GeneName38)&as.character(Chr37)!=as.character(Chr38))) # 16 genes, all have same name
# View(ens37to38Mapping %>% filter(is.na(GeneName38)))

ens37to38MapByNameLoc = merge(ensemblGeneData37 %>% select(GeneId37=GeneId,GeneName,Chromosome,KarBand37=KaryotypeBand),
                              ensemblGeneData38 %>% select(GeneId38=GeneId,GeneName,Chromosome,KarBand38=KaryotypeBand),by=c('GeneName','Chromosome'),all.x=T)

nrow(ens37to38MapByNameLoc)
View(ens37to38MapByNameLoc %>% filter(is.na(GeneId38)))

# 876 where name and chromosome are the same
View(ens37to38MapByNameLoc %>% filter(!is.na(GeneId38)&as.character(GeneId37)!=as.character(GeneId38)))

nrow(ens37to38MapByNameLoc %>% filter(!is.na(GeneId38)&as.character(GeneId37)!=as.character(GeneId38)&as.character(KarBand37)==as.character(KarBand38))) # 783 KB match, 93 don't
View(ens37to38MapByNameLoc %>% filter(!is.na(GeneId38)&as.character(GeneId37)!=as.character(GeneId38)&as.character(KarBand37)!=as.character(KarBand38))) # 93 don't but are typically close

# map by name only where geneId and chromosome aren't the same
ens37to38MapByName = merge(ensemblGeneData37 %>% select(GeneId37=GeneId,GeneName,Chr37=Chromosome,KarBand37=KaryotypeBand),
                           ensemblGeneData38 %>% select(GeneId38=GeneId,GeneName,Chr38=Chromosome,KarBand38=KaryotypeBand),by=c('GeneName'),all.x=T)

View(ens37to38MapByName %>% filter(!is.na(GeneId38)&as.character(Chr37)!=as.character(Chr38)))




#####
# Driver genes
#####

# driver genes will now use the same gene names for both v37 ad v38 as well

driverGenes37 = read.csv('~/hmf/repos/common-resources-public/gene_panel/37/DriverGenePanel.37.tsv',sep='\t')
nrow(driverGenes37)
View(driverGenes37)

View(driverGenes37 %>% filter(!(gene %in% ensemblGeneData37$GeneName)))


driverGene37NameChanges = merge(driverGenes37 %>% select(GeneName=gene),ensemblGeneData37Old %>% select(GeneId,GeneName),by='GeneName',all.x=T)
driverGene37NameChanges = merge(driverGene37NameChanges,ensemblGeneData37 %>% select(GeneId,GeneNameNew=GeneName),by='GeneId',all.x=T)
View(driverGene37NameChanges)
View(driverGene37NameChanges %>% filter(as.character(GeneName)!=as.character(GeneNameNew)))

driverGenes38 = read.csv('~//hmf/repos/common-resources-public/gene_panel/38/DriverGenePanel.38.tsv',sep='\t')
nrow(driverGenes38)

View(driverGenes38 %>% filter(!(gene %in% ensemblGeneData38$GeneName)))

ensemblGeneData38Old = read.csv('~/data/ensembl_hg38/ensembl_gene_data.csv')
View(driverGenes38 %>% filter(!(gene %in% ensemblGeneData38Old$GeneName)))

driverGene38NameChanges = merge(driverGenes38 %>% select(GeneName=gene),ensemblGeneData38Old %>% select(GeneId,GeneName),by='GeneName',all.x=T)
driverGene38NameChanges = merge(driverGene38NameChanges,ensemblGeneData38 %>% select(GeneId,GeneNameNew=GeneName),by='GeneId',all.x=T)
View(driverGene38NameChanges)
View(driverGene38NameChanges %>% filter(as.character(GeneName)!=as.character(GeneNameNew)))
View(driverGenes38 %>% filter(!(gene %in% ensemblGeneData38Old$GeneName)))


# manually name-changed gene: 
View(ensemblGeneData38Old %>% filter(GeneName=='H3C2'))
View(ensemblGeneData38Old %>% filter(GeneName=='HIST1H3B'))

# Charless-MBP:ensembl_db charlesshale$ grep ENSG00000286522 ./ensembl_38_104/ensembl_gene_data.csv
#ENSG00000286522,H3C2,6,-1,26031589,26032099,p22.2,HGNC:4776
# ENSG00000274267,HIST1H3B,6,-1,26031650,26032060,p22.2,HGNC:4776


#####
# Fusion Genes
#####


# known fusion data will now use the same gene names for both v37 ad v38 as well

knownFusionData37 = read.csv('~/data/fusion_ref/known_fusion_data.37_v3.csv')
View(knownFusionData37)

fusionGeneData37 = rbind(knownFusionData37 %>% filter(FiveGene!=''&!(FiveGene %in% c('IGH','IGL','IGK'))) %>% mutate(Stream='Five') %>% select(GeneName=FiveGene,Stream,KnownType=Type),
                         knownFusionData37 %>% filter(ThreeGene!='') %>% mutate(Stream='Three') %>% select(GeneName=ThreeGene,Stream,KnownType=Type))

fusionGeneNameChanges37 = fusionGeneData37 %>% filter(!(GeneName %in% ensemblGeneData37$GeneName))
View(fusionGeneNameChanges37)
fusionGeneNameChanges37 = merge(fusionGeneNameChanges37,ensemblGeneData37Old %>% select(GeneId,GeneName),by='GeneName',all.x=T)
fusionGeneNameChanges37 = merge(fusionGeneNameChanges37,ensemblGeneData37 %>% select(GeneId,GeneNameNew=GeneName),by='GeneId',all.x=T)
View(fusionGeneNameChanges37)
View(fusionGeneData37 %>% filter(!(GeneName %in% ensemblGeneData37$GeneName)))

knownFusionData38 = read.csv('~/data/fusion_ref/known_fusion_data.38_v3.csv')
View(knownFusionData38)

fusionGeneData38 = rbind(knownFusionData38 %>% filter(FiveGene!=''&!(FiveGene %in% c('IGH','IGL','IGK'))) %>% mutate(Stream='Five') %>% select(GeneName=FiveGene,Stream,KnownType=Type),
                         knownFusionData38 %>% filter(ThreeGene!='') %>% mutate(Stream='Three') %>% select(GeneName=ThreeGene,Stream,KnownType=Type))

fusionGeneNameChanges38 = fusionGeneData38 %>% filter(!(GeneName %in% ensemblGeneData38$GeneName))
View(fusionGeneNameChanges38)
fusionGeneNameChanges38 = merge(fusionGeneNameChanges38,ensemblGeneData38Old %>% select(GeneId,GeneName),by='GeneName',all.x=T)
fusionGeneNameChanges38 = merge(fusionGeneNameChanges38,ensemblGeneData38 %>% select(GeneId,GeneNameNew=GeneName),by='GeneId',all.x=T)
View(fusionGeneNameChanges38)

View(ensemblGeneData37Old)
View(ensemblGeneData38Old)

# RP11-356O9.1
