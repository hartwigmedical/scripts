SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE sampleId LIKE '%CPCT0202%'
ORDER BY 1;