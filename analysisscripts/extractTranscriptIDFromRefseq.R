#!/usr/bin/Rscript

library("biomaRt")

grch37 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
refSeqNKI = c("XXXX")

getBM(attributes = c('external_gene_name', 'refseq_mrna', 'ensembl_transcript_id', 'chromosome_name', 'band'),
filters = 'refseq_mrna',
values = refSeqNKI,
mart = grch37)

