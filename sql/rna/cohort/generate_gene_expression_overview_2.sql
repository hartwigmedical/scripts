SELECT g.sampleId AS "Sample ID", gene AS "Gene",
primaryTumorLocation AS "Primary tumor location", IF(primaryTumorSubType = "", primaryTumorType, primaryTumorSubType) AS "Primary tumor type",
purity AS "Tumor purity", p.qcStatus AS "WGS QC status", r.qcStatus AS "RNA QC status", round(tpm,3) AS "TPM",
round(medianTpmCancer,2) AS "Median TPM for cancer type", round(percentileCancer,2) AS "Percentile for cancer type",
round(medianTpmCohort,2) AS "Median TPM for database", round(percentileCohort,2) AS "Percentile for database",
IF(percentileCancer>=0.5,1,0) AS "Top 50% for cancer type?"
FROM geneExpression g
INNER JOIN sample s ON s.sampleID=g.sampleId
LEFT JOIN baseline b ON s.patientId=b.patientId
INNER JOIN rnaStatistics r ON g.sampleId=r.sampleId
INNER JOIN purity p ON g.sampleId=p.sampleId
WHERE g.sampleId LIKE 'XXX%' AND gene IN ('XXX') AND primaryTumorLocation IN ('XXX') AND p.qcStatus NOT LIKE '%FAIL_NO_TUMOR%' ORDER BY g.sampleId;