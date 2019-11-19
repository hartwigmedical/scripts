SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation IN ('Breast', 'Prostate') AND sampleId LIKE '%CPCT%'
ORDER BY 1;
