SELECT * FROM geneExpression g
WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='TSG')
AND (percentileCohort < 0.05 AND percentileCancer < 0.05)
AND sampleId in ('XXX')
ORDER BY percentileCohort;