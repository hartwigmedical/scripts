SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation REGEXP 'Esophagus|Pancreas|Stomach'
ORDER BY 1;