SELECT a.*, IF(driver IS NOT NULL, group_concat(driver," (",round(driverLikelihood,2),")"), "") AS hasDriverInGene, round(g.minCopyNumber,2) AS minCopyNumber, round(g.maxCopyNumber,2) AS maxCopyNumber, purity
FROM
(SELECT primaryTumorLocation, primaryTumorType, g.* FROM geneExpression g INNER JOIN datarequest dr ON g.sampleId=dr.sampleId
WHERE (primaryTumorLocation, g.gene) IN
	(SELECT primaryTumorLocation, gene FROM geneExpression g INNER JOIN clinical c ON g.sampleId=c.sampleId
	WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='TSG')
	AND percentileCancer < 0.05
	AND g.sampleId IN ('XXX'))
) AS a
LEFT JOIN driverCatalog d ON d.sampleId=a.sampleId AND d.gene=a.gene
INNER JOIN geneCopyNumber g ON g.sampleId=a.sampleId AND g.gene=a.gene
INNER JOIN purity p ON a.sampleId=p.sampleId
WHERE percentileCancer < 0.05 AND a.sampleId IN (SELECT sampleId FROM rnaStatistics WHERE qcStatus='PASS')
GROUP BY gene, sampleId
ORDER BY gene, percentileCancer;