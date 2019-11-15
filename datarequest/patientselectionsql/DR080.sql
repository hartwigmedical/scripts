SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation IN ('Stomach', 'Colon/Rectum', 'Esophagus') AND
(biopsyPostDrugMechanisms LIKE '%anti-pd-1%' OR biopsyPostDrugMechanisms LIKE '%anti-pd-l1%' OR biopsyPostDrugMechanisms LIKE '%anti-ctla-4%')
ORDER BY 1;
