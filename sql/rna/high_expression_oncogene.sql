SELECT * FROM geneExpression g
WHERE gene IN (SELECT gene FROM driverGenePanel WHERE likelihoodType='ONCO')
AND (percentileCohort > 0.90 AND percentileCancer > 0.90)
AND sampleId in ('XXX')
ORDER BY percentileCohort DESC;