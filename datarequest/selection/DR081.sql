SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Biliary' AND cancerSubtype = 'Gall bladder'
ORDER BY 1;
