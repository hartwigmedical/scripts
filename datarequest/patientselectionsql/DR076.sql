SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE cancerSubtype LIKE '%glioblastoma%' AND treatmentName LIKE '%temozolomide%'
ORDER BY 1;
