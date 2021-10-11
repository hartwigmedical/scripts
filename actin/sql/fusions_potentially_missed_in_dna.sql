SELECT primaryTumorLocation, cuppaTumorLocation, b.reportedType, dnaFusionCount, a.*
FROM
(SELECT * FROM rnaFusion WHERE name IN (SELECT name FROM svFusion WHERE reportedType <> 'NONE'))
	AS a
LEFT JOIN (SELECT distinct name, reportedType, count(*) as dnaFusionCount FROM svFusion WHERE reportedType <> 'NONE' GROUP BY 1,2) AS b
    ON a.name=b.name
LEFT JOIN svFusion
    ON a.sampleId=svFusion.sampleId and b.name=svFusion.name
LEFT JOIN clinical
    on a.sampleId=clinical.sampleId
LEFT JOIN cuppa
    on a.sampleId=cuppa.sampleId
WHERE svFusion.sampleId IS NULL AND a.sampleId IN ('XXX')
ORDER BY b.reportedType, a.name, a.sampleId;