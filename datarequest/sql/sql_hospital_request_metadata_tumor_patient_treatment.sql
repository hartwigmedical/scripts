SELECT
    clinical.patientId AS '#patientId',
    clinical.sampleId,
    clinical.setName,
    IFNULL(datarequest.tumorPurity, 'na') as tumorPurity,
    clinical.hmfPatientId,
    clinical.hmfSampleId,
    clinical.primaryTumorLocation,
    clinical.cancerSubtype,
    clinical.biopsyDate,
    clinical.biopsySite,
    clinical.biopsyLocation,
    clinical.gender,
    clinical.birthYear,
    clinical.deathDate,
    clinical.hasSystemicPreTreatment,
    clinical.hasRadiotherapyPreTreatment,
    clinical.treatmentGiven,
    clinical.treatmentStartDate,
    clinical.treatmentEndDate,
    clinical.treatment,
    clinical.consolidatedTreatmentType as treatmentType,
    clinical.responseDate,
    clinical.responseMeasured,
    clinical.firstResponse
FROM clinical
LEFT JOIN datarequest
ON clinical.sampleId = datarequest.sampleId;