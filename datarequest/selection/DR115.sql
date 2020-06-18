SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation REGEXP 'Breast|Colon/Rectum|Esophagus|Stomach|Pancreas|Small intestine|Liver|Biliary'
ORDER BY 1;