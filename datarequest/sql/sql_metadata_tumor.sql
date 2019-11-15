SELECT
    patientId AS '#patientId',
    sampleId,
    setName,
    tumorPurity,
    hmfPatientId,
    hmfSampleId,
    primaryTumorLocation,
    cancerSubtype,
    biopsyDate,
    biopsySite,
    biopsyLocation
FROM datarequest;