SELECT * FROM geneExpression g
WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='TSG')
AND (percentileCohort < 0.05 OR percentileCancer < 0.1)
AND sampleId in ('XXX')
ORDER BY percentileCohort DESC;