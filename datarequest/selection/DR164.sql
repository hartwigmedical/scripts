SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation='Nervous system' and cancerSubtype = 'Glioblastoma Multiforme'
ORDER BY 1;