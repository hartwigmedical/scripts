SELECT a.*, IF(driver IS NOT NULL, group_concat(driver," (",round(driverLikelihood,2),")"), "") AS hasDriverInGene, round(g.maxCopyNumber,2) as maxCopyNumber, purity
FROM
(SELECT primaryTumorLocation, primaryTumorSubLocation, primaryTumorType, primaryTumorSubType, g.* FROM geneExpression g INNER JOIN datarequest dr ON g.sampleId=dr.sampleId
WHERE g.gene IN
	(SELECT gene FROM geneExpression g
	WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='ONCO')
	AND percentileCohort > 0.95
	AND g.sampleId IN ('XXX')
    GROUP BY g.sampleId, gene)
) AS a
LEFT JOIN driverCatalog d ON d.sampleId=a.sampleId AND d.gene=a.gene
INNER JOIN geneCopyNumber g ON g.sampleId=a.sampleId AND g.gene=a.gene
INNER JOIN purity p ON a.sampleId=p.sampleId
WHERE percentileCohort > 0.95
GROUP BY gene, sampleId
ORDER BY gene, tpm DESC;