SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation = 'Lung'
AND primaryTumorSubType like '%Non-Small Cell%'
AND concatenatedTreatmentType LIKE '%Immunotherapy%'
ORDER BY 1;