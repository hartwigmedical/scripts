SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Kidney' AND cancerSubtype = 'Renal cell'
ORDER BY 1;
