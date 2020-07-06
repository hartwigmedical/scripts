SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation='Lung'
ORDER BY 1;