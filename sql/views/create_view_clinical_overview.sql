DROP VIEW IF EXISTS clinical;

CREATE VIEW clinical AS

	SELECT sample.sampleId, sample.arrivalDate as sampleArrivalDate, patient.cpctId AS patientId, patient.hospital, patient.gender, patient.birthYear, patient.registrationDate, biopsy.biopsyDate, 
	treatment.startDate AS treatmentStartDate, treatment.endDate AS treatmentEndDate, firstMatchedTreatmentResponse.responseDate, patient.deathDate, patient.cancerType,
    patient.cancerSubtype, biopsy.biopsySite, biopsy.biopsyLocation, treatment.treatmentGiven, treatment.name AS treatment, treatment.type AS treatmentType,
    firstMatchedTreatmentResponse.measurementDone as responseMeasured, firstMatchedTreatmentResponse.response
	FROM sample 
		INNER JOIN patient ON sample.patientId = patient.id
		LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
		LEFT JOIN treatment on treatment.biopsyId = biopsy.id
		LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId
		