SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Nervous system' AND
( preTreatments like '%Temozolomide%' OR treatment like '%Temozolomide%')
ORDER BY 1;