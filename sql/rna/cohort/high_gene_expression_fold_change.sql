SELECT r.qcStatus, primaryTumorLocation, purity, g.sampleId, g.gene, round(gn.minCopyNumber,1) AS minCopyNr, ploidy, medianTpmCancer, medianTpmCohort, tpm, percentileCancer, percentileCohort, round(tpm/medianTpmCancer,1) AS foldChangeCancer, round(tpm/medianTpmCohort,1) AS foldChangeCohort, group_concat(distinct d.driver) AS driversInGene, group_concat(distinct dr.gene), group_concat(distinct name) AS reportedFusions, hrStatus
FROM geneExpression g
LEFT JOIN driverCatalog d ON g.sampleId=d.sampleId AND g.gene=d.gene
LEFT JOIN (SELECT * FROM driverCatalog WHERE driverLikelihood=1) AS dr ON g.sampleId=dr.sampleId
LEFT JOIN (SELECT sampleId, name FROM svFusion WHERE reported) AS b ON b.sampleId=g.sampleId
LEFT JOIN chord ch ON g.sampleId=ch.sampleId
INNER JOIN geneCopyNumber gn ON g.sampleId=gn.sampleId AND g.gene=gn.gene
INNER JOIN purity pu ON g.sampleId=pu.sampleId
INNER JOIN clinical c ON g.sampleId=c.sampleId
INNER JOIN rnaStatistics r ON r.sampleId=g.sampleId
WHERE r.qcStatus = 'PASS' AND pu.qcStatus='PASS' AND g.gene IN ('XXX')
GROUP BY 4,5
HAVING (foldChangeCohort > 10 OR foldChangeCancer > 10)
ORDER BY gene, tpm DESC;
