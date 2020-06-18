SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE
primaryTumorLocation = 'Colon/Rectum'
AND (concatenatedTreatmentType like "%chemo%")
ORDER BY 1;
