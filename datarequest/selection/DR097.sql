SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Ovary' AND
( preTreatments like '%Carboplatin%' OR  preTreatments like  '%Cisplatin%' OR preTreatments like '%Olaparib%'
OR treatment like '%Carboplatin%' OR treatment like  '%Cisplatin%' OR treatment like '%Olaparib%')
ORDER BY 1;