SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE
primaryTumorLocation='Skin'
AND primaryTumorType='Melanoma'
ORDER BY 1;