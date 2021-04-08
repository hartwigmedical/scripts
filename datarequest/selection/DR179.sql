SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation in ('Breast',  'Colorectum', 'Uterus', 'Ovary', 'Pancreas')
ORDER BY 1;
