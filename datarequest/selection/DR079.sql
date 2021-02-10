SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Skin' AND primaryTumorType = 'Melanoma' AND (concatenatedTreatmentType LIKE '%Immuno%' OR concatenatedTreatmentType LIKE '%Target%')
ORDER BY 1;
