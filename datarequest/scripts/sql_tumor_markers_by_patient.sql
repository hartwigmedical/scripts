SELECT
	patient.patientIdentifier AS '#patientId',
    tumorMarker.date,
    tumorMarker.marker,
    tumorMarker.measurement,
    tumorMarker.unit
FROM tumorMarker
	INNER JOIN patient ON tumorMarker.patientId = patient.id
ORDER BY 1,2;
