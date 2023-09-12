SELECT g.sampleId AS "Sample ID", gene AS "Gene", primaryTumorLocation AS "Primary tumor location",
purity AS "Tumor purity", p.qcStatus AS "WGS QC status", r.qcStatus AS "RNA QC status", round(tpm,3) AS "TPM",
round(medianTpmCancer,2) AS "Median TPM for cancer type", round(percentileCancer,2) AS "Percentile for cancer type",
round(medianTpmCohort,2) AS "Median TPM for database", round(percentileCohort,2) AS "Percentile for database",
IF((percentileCancer>=0.8 OR percentileCohort>=0.8),1,0) AS "Top 20% tumor type or database?",
IF((percentileCancer>=0.75 OR percentileCohort>=0.75),1,0) AS "Top 25% tumor type or database?"
FROM geneExpression g
INNER JOIN sample s ON s.sampleId=g.sampleId
LEFT JOIN baseline b ON s.patientId=b.patientId
INNER JOIN rnaStatistics r ON g.sampleId=r.sampleId
INNER JOIN purity p ON g.sampleId=p.sampleId
WHERE g.sampleId LIKE 'XXX%' AND gene IN ('XXX') AND p.qcStatus NOT LIKE '%FAIL_NO_TUMOR%' ORDER BY g.sampleId;