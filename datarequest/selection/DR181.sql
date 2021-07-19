SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE purpleVersion < 2.54
ORDER BY 1;
