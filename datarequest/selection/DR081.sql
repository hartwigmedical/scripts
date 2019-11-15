SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation = 'Biliary' AND cancerSubtype = 'Gall bladder'
ORDER BY 1;
