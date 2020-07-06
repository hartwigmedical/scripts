SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE
primaryTumorLocation='Skin'
AND cancerSubtype='Melanoma'
ORDER BY 1;