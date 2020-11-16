SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE (left(sampleId, 4) = 'WIDE' or hospital in ('ANTONI VAN LEEUWENHOEK ZIEKENHUIS', 'NKI-AVL, Amsterdam')) and concatenatedTreatmentType like '%Chemo%'
ORDER BY 1;