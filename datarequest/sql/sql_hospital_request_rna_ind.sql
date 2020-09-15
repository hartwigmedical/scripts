SELECT
    patientId AS '#patientId',
    sampleId,
    hmfPatientId,
    hmfSampleId
FROM clinical WHERE hasRNA=1;