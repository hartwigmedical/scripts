SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Urinary tract' and cancerSubtype = 'Bladder'
ORDER BY 1;