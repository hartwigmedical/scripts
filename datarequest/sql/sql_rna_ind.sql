SELECT
    patientId AS '#patientId',
    sampleId,
    hmfPatientId,
    hmfSampleId
FROM datarequest WHERE hasRNA=1;