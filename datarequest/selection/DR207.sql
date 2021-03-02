SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest where hospital = 'LUMC'
ORDER BY 1;