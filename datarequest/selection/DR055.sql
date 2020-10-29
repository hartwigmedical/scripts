SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest
WHERE primaryTumorLocation in ('Bile duct',  'Hepatobiliary system', 'Gallbladder')
ORDER BY 1;