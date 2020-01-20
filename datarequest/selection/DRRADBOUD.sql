SELECT
	DISTINCT patientId AS '#patientId'
FROM
	datarequest
WHERE sampleId LIKE '%CPCT0207%' OR sampleId LIKE '%DRUP0107%'
ORDER BY 1;