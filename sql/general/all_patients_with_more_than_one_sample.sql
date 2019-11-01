SELECT count(sampleId), if(isnull(targetId), patientId, targetId) AS effectivePatientId FROM
(SELECT sampleId, patientId, targetId FROM datarequestFiltered
LEFT JOIN patientMapping on datarequestFiltered.patientId = patientMapping.sourceId) AS x
GROUP BY 2 HAVING count(sampleId) > 1 ORDER BY 1 DESC;