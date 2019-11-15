SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Breast' AND sampleId IN (SELECT sampleId FROM rna)
ORDER BY 1;
