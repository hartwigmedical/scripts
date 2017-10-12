SELECT c.sampleId, 
round(pow(SUM(IF(start<>1 AND structuralVariantSupport='NONE',1,0)),3)/pow(SUM(IF(start<>1,1,0)),2)/AVG(ploidy),0) AS QCscore
FROM copyNumber c INNER JOIN purity p ON c.sampleId = p.sampleId
WHERE c.sampleId IN ('xxx')
GROUP BY c.sampleId