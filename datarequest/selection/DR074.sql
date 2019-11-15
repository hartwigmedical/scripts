SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE sampleId IN (SELECT sampleId FROM rna)
ORDER BY 1;
