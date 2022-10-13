
driverGeneFile <- '/data/resources/public/gene_panel/37/DriverGenePanel.37.tsv'
ensemblGeneDataFile <- '/data/resources/public/ensembl_data_cache/37/ensembl_gene_data.csv'
ensemblTransExonDataFile <- '/data/resources/public/ensembl_data_cache/37/ensembl_trans_exon_data.csv'

ensemblGeneData <- read.csv(ensemblGeneDataFile)
ensemblTransExonData <- read.csv(ensemblTransExonDataFile)
driverGenePanel <- read.csv(driverGeneFile,sep='\t')
driverGeneData <- ensemblGeneData %>% filter(GeneName %in% driverGenePanel$gene) %>% select(GeneId,GeneName)
driverTranscripts <- ensemblTransExonData %>% filter(ExonRank==1) %>% filter(CanonicalTranscriptId==TransId | TransName=='ENST00000579755') # manually include the CDKN2A alt transcript
driverGeneTransData <- merge(driverGeneData,driverTranscripts %>% select(GeneId,TransName),by='GeneId',all.x=T)

# outputFile <-  paste0('hmf_reportable_gene_transcripts.tsv')
# write.table(driverGeneTransData, file = "outputFile", row.names=TRUE, sep="\t")