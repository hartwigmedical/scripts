SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE concatenatedTreatmentType LIKE "%Targeted therapy%"
ORDER BY 1;