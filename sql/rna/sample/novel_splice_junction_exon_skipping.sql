SELECT geneId, transcriptId, primaryTumorLocation, a.*, driversInGene, reportedFusions FROM
(SELECT * FROM novelSpliceJunction
WHERE gene in ('EGFR','ALK','BRAF','KIT','MET','FLT3','BRAF','CTNNB1','EGFR','AHR','PDGFRA','KMT2A') AND type IN ('SKIPPED_EXONS') AND fragmentCount>5 AND cohortFrequency<30) AS a
LEFT JOIN (SELECT sampleId, gene, group_concat(likelihoodMethod) AS driversInGene FROM driverCatalog d GROUP BY 1,2) AS b
ON a.sampleId=b.sampleId AND a.gene=b.gene
LEFT JOIN (SELECT sampleId, group_concat(name) AS reportedFusions FROM svFusion WHERE reported OR name IN ('EGFR_EGFR','ALK_ALK','BRAF_BRAF','KIT_KIT','MET_MET','FLT3_FLT3','BRAF_BRAF','CTNNB1_CTNNB1','EGFR_EGFR','AHR_AHR','PDGFRA_PDFGRA','KMT2A_KMT2A') GROUP BY 1) AS c
ON a.sampleId=c.sampleId
LEFT JOIN clinical
ON a.sampleId=clinical.sampleId
LEFT JOIN canonicalTranscript
ON a.gene=canonicalTranscript.gene
WHERE a.sampleId in ('XXX')
GROUP BY a.sampleId, a.gene
ORDER BY a.gene, junctionStart;