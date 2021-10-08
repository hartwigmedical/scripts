SELECT * FROM geneExpression g
WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='ONCO')
AND (percentileCohort > 0.95 OR percentileCancer > 0.9)
AND sampleId in ('XXX')
ORDER BY percentileCohort DESC;