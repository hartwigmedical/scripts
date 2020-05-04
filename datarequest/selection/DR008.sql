SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE
    primaryTumorLocation='Lung'
    AND cancerSubtype='Non-Small Cell'
    AND (concatenatedTreatmentType LIKE '%Targeted%' OR concatenatedTreatmentType LIKE '%Immuno%')
ORDER BY 1;
