SELECT
    patientId AS '#patientId',
    sampleId,
    setName,
    tumorPurity,
    hmfPatientId,
    hmfSampleId,
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