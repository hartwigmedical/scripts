SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation = 'Lung'
AND cancerSubtype = 'Non-Small Cell'
AND concatenatedTreatmentType LIKE '%Immunotherapy%'
ORDER BY 1;