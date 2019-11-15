SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation = 'Breast' AND sampleId IN (SELECT sampleId FROM rna)
ORDER BY 1;
