SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE sampleId IN (SELECT sampleId FROM rna)
ORDER BY 1;
