SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = "Urothelial tract" AND (treatment like "%pembro%")
ORDER BY 1;
