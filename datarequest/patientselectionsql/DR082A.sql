SELECT    
    DISTINCT patientId AS '#patientId'
FROM
    datarequestFiltered
WHERE biopsyPostDrugMechanisms LIKE '%taxane%' AND sampleId LIKE '%CPCT%'
ORDER BY 1;
