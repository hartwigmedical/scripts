SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE cancerSubtype LIKE '%glioblastoma%' AND treatmentName LIKE '%temozolomide%'
ORDER BY 1;
