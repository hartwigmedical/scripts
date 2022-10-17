library(dplyr)

ensemblGeneData <- read.csv('/data/resources/public/ensembl_data_cache/37/ensembl_gene_data.csv')
ensemblTransExonData <- read.csv('/data/resources/public/ensembl_data_cache/37/ensembl_trans_exon_data.csv')
driverGenePanel <- read.csv('/data/resources/public/gene_panel/37/DriverGenePanel.37.tsv',sep='\t')
driverGeneData <- ensemblGeneData %>% filter(GeneName %in% driverGenePanel$gene) %>% select(GeneId, GeneName)
driverTranscripts <- ensemblTransExonData %>% filter(ExonRank==1) %>% filter(CanonicalTranscriptId==TransId | TransName=='ENST00000579755') # manually include the CDKN2A alt transcript
driverGeneTransData <- merge(driverGeneData,driverTranscripts %>% select(GeneId, TransName),by='TransName',all.x=T)

outputFile <-  paste0('/home/lschoenmaker/tmphmf_reportable_gene_transcripts.tsv')
write.table(driverGeneTransData, file = outputFile, row.names=TRUE, sep="\t")