SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE primaryTumorLocation = 'Skin' AND cancerSubType = 'Melanoma' AND (biopsyPostDrugTypes LIKE '%Immuno%' OR biopsyPostDrugTypes LIKE '%Target%')
ORDER BY 1;
