SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE (primaryTumorLocation = 'Skin' and cancerSubtype = 'Melanoma') or (primaryTumorLocation = 'Lung' and cancerSubType = 'Non-small Cell')
ORDER BY 1;
