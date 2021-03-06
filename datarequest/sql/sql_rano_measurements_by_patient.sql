SELECT
	patient.patientIdentifier AS '#patientId',
    ranoMeasurement.responseDate,
    ranoMeasurement.therapyGiven,
    ranoMeasurement.targetLesionResponse,
    ranoMeasurement.noTargetLesionResponse,
    ranoMeasurement.overallResponse
FROM ranoMeasurement
	INNER JOIN patient ON ranoMeasurement.patientId = patient.id
ORDER BY 1,2;
