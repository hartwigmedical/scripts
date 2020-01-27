SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Ovary' AND DATE(informedConsentDate) > '2015-12-31'
ORDER BY 1;
