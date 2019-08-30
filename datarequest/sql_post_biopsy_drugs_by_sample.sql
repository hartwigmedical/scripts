SELECT
    patient.patientIdentifier AS '#patientId',
    biopsy.sampleId,
    drug.startDate,
    drug.endDate,
    drug.name,
    drug.type,
    drug.mechanism
FROM drug
    INNER JOIN treatment ON drug.treatmentId = treatment.id
    INNER JOIN biopsy ON treatment.biopsyId = biopsy.id
    INNER JOIN patient ON biopsy.patientId = patient.id
WHERE treatmentGiven = "Yes" AND NOT isnull(drug.startDate)
ORDER BY 1,2,3;