SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Prostate' AND
( treatment like '%Enzalutamide%' OR  treatment like '%Abirateron%')
ORDER BY 1;