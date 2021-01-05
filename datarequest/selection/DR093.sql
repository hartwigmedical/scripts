SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Nervous system' AND primaryTumorSubType = 'Glioblastoma Multiforme' AND (treatment like '%Temozolomide%')
ORDER BY 1;