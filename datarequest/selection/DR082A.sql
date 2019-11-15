SELECT    
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE biopsyPostDrugMechanisms LIKE '%taxane%' AND sampleId LIKE '%CPCT%'
ORDER BY 1;
