DROP VIEW IF EXISTS clinical;

CREATE VIEW clinical AS

	SELECT sample.sampleId, sample.arrivalDate as sampleArrivalDate, patient.patientIdentifier AS patientId, baseline.informedConsentDate, baseline.hospital, baseline.gender,
    baseline.birthYear, baseline.hasSystemicPreTreatment, baseline.hasRadiotherapyPreTreatment, baseline.preTreatments, biopsy.biopsyDate,
	treatment.startDate AS treatmentStartDate, treatment.endDate AS treatmentEndDate, firstMatchedTreatmentResponse.responseDate, baseline.deathDate,
	baseline.primaryTumorLocation, baseline.cancerSubtype, biopsy.biopsyType, biopsy.biopsySite, biopsy.biopsyLocation,
	treatment.treatmentGiven, treatment.radiotherapyGiven, treatment.name AS treatment, treatment.type AS treatmentType,
	firstMatchedTreatmentResponse.measurementDone as responseMeasured, firstMatchedTreatmentResponse.response as firstResponse
	FROM sample
		INNER JOIN patient ON sample.patientId = patient.id
        LEFT JOIN baseline ON patient.id = baseline.patientId
		LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
		LEFT JOIN treatment on treatment.biopsyId = biopsy.id
		LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId;