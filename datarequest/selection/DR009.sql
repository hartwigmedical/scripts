SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE concatenatedTreatmentType LIKE "%Immuno%"
ORDER BY 1;