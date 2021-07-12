SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Colorectum' or primaryTumorLocation = 'Lung' or primaryTumorLocation = 'Breast'
ORDER BY 1;