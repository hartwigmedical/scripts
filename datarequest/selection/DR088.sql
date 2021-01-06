SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Kidney' AND primaryTumorSubType like '%renal cell%'
ORDER BY 1;