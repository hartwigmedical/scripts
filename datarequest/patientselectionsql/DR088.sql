SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation = 'Kidney' AND cancerSubtype = 'Renal cell'
ORDER BY 1;
