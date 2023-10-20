#!/usr/bin/Rscript

library(dplyr)

ensemblGeneData <- read.csv('/data/resources/report_resources/37/ensembl_gene_data.csv')
ensemblTransExonData <- read.csv('/data/resources/report_resources/37/ensembl_trans_exon_data.csv')
driverGenePanel <- read.csv('/data/resources/report_resources/37/DriverGenePanel.37.tsv',sep='\t')
driverGeneData <- ensemblGeneData %>% filter(GeneName %in% driverGenePanel$gene) %>% select(GeneId,GeneName)
driverTranscripts <- ensemblTransExonData %>% filter(ExonRank==1) %>% filter(CanonicalTranscriptId==TransId | TransName=='ENST00000579755') # manually include the CDKN2A alt transcript
driverGeneTransData <- merge(driverGeneData,driverTranscripts %>% select(GeneId,TransName),by='GeneId',all.x=T)

outputFile <-  paste0('/data/resources/report_resources/oncoact_wgs_gene_transcripts.tsv')
write.table(driverGeneTransData, file = outputFile, row.names=TRUE, sep="\t")