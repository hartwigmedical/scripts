SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE sampleId IN ( # See DR archive)
ORDER BY 1
