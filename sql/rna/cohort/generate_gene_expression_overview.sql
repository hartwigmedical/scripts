SELECT purity AS "Tumor purity", r.qcStatus AS "RNA QC status",
primaryTumorLocation AS "Primary tumor location", g.sampleId AS "Sample ID", gene AS "Gene", round(tpm,3) AS "TPM",
round(medianTpmCancer,2) AS "Median TPM for cancer type", round(percentileCancer,2) AS "Percentile for cancer type",
round(medianTpmCohort,2) AS "Median TPM for database", round(percentileCohort,2) AS "Percentile for database",
IF((percentileCancer>=0.8 OR percentileCohort>=0.8),1,0) AS "Top 20% tumor type or database?",
IF((percentileCancer>=0.75 OR percentileCohort>=0.75),1,0) AS "Top 25% tumor type or database?"
FROM geneExpression g
INNER JOIN clinical c ON g.sampleId=c.sampleId
INNER JOIN rnaStatistics r ON g.sampleId=r.sampleId
INNER JOIN purity p ON g.sampleId=p.sampleId
WHERE g.sampleId LIKE 'XXX%' AND gene IN ('XXX') ORDER BY g.sampleId;