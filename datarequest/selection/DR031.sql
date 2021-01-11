SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE
    primaryTumorLocation = "Urothelial tract"
ORDER BY 1;