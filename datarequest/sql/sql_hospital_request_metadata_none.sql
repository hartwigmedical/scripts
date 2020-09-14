SELECT
    clinical.patientId AS '#patientId',
    clinical.sampleId,
    clinical.setName,
    IFNULL(datarequest.tumorPurity, 'na') as tumorPurity,
    clinical.hmfPatientId,
    clinical.hmfSampleId
FROM clinical
LEFT JOIN datarequest
ON clinical.sampleId = datarequest.sampleId;