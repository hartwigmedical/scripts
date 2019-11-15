SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation = 'Prostate'
ORDER BY 1;
