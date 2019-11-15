SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE sampleId IN ( # See DR archive)
ORDER BY 1
