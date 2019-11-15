SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation IN ('Breast', 'Prostate') AND sampleId LIKE '%CPCT%'
ORDER BY 1;
