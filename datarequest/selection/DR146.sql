SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE concatenatedTreatmentType like "%Targeted therapy%"
ORDER BY 1;