SELECT
    patientId AS '#patientId',
    sampleId,
    setName,
    tumorPurity,
    hmfPatientId,
    hmfSampleId,
    primaryTumorLocation,
    cancerSubtype,
    primaryTumorLocation,
    primaryTumorSubLocation,
    primaryTumorType,
    primaryTumorSubType,
    primaryTumorExtraDetails,
    doids,
    biopsyDate,
    biopsySite,
    biopsyLocation
FROM datarequest;