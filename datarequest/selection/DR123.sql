SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE hasRNA = 1
ORDER BY 1;