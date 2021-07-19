SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE (primaryTumorLocation = 'Colorectum' or primaryTumorLocation = 'Lung' or primaryTumorLocation = 'Breast') AND purpleVersion < 2.54
ORDER BY 1;