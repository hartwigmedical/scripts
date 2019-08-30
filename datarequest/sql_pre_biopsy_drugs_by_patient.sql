SELECT
    patient.patientIdentifier AS '#patientId',
    preTreatmentDrug.startDate,
    preTreatmentDrug.endDate,
    preTreatmentDrug.name,
    preTreatmentDrug.type,
    preTreatmentDrug.mechanism
FROM preTreatmentDrug
    INNER JOIN patient ON preTreatmentDrug.patientId = patient.id
WHERE NOT isnull(preTreatmentDrug.startDate)
ORDER BY 1,2;