SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Colon/Rectum' or primaryTumorLocation = 'Lung' or primaryTumorLocation = 'Breast'
ORDER BY 1;