SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE patientId IN (# See DR archive)
ORDER BY 1;