SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation = 'Prostate'
AND (preTreatments REGEXP 'Abirateron|Enzalutamide' OR treatment REGEXP 'Abirateron|Enzalutamide')
ORDER BY 1;