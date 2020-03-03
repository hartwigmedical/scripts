SELECT
    patient.patientIdentifier AS '#patientId',
    biopsy.sampleId,
    treatmentResponse.responseDate,
    treatmentResponse.response,
    treatment.startDate,
    treatment.endDate,
    treatment.name,
    treatment.type,
    treatment.mechanism
FROM treatmentResponse
    INNER JOIN treatment ON treatmentResponse.treatmentId = treatment.id
    INNER JOIN biopsy ON treatment.biopsyId = biopsy.id
    INNER JOIN patient ON biopsy.patientId = patient.id
    INNER JOIN datarequest ON biopsy.sampleId = datarequest.sampleId
WHERE treatment.treatmentGiven = "Yes" AND treatmentResponse.measurementDone = "Yes"
ORDER BY 1,2,3;
