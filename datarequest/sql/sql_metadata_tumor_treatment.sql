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
    biopsyLocation,
    hasSystemicPreTreatment,
    hasRadiotherapyPreTreatment,
    treatmentGiven,
    treatmentStartDate,
    treatmentEndDate,
    treatment,
    consolidatedTreatmentType as treatmentType,
    responseDate,
    responseMeasured,
    firstResponse
FROM datarequest;