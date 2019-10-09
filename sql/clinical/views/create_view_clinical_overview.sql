CREATE OR REPLACE VIEW clinical AS

	SELECT sample.sampleId, sample.arrivalDate as sampleArrivalDate, patient.patientIdentifier AS patientId,
	baseline.registrationDate, baseline.informedConsentDate, baseline.hospital, baseline.gender,
    baseline.birthYear, baseline.hasSystemicPreTreatment, baseline.hasRadiotherapyPreTreatment, baseline.preTreatments, baseline.preTreatmentsType,
    baseline.preTreatmentsMechanism, biopsy.biopsyDate,
	treatment.startDate AS treatmentStartDate, treatment.endDate AS treatmentEndDate, firstMatchedTreatmentResponse.responseDate, baseline.deathDate,
	baseline.primaryTumorLocation, baseline.cancerSubtype, biopsy.biopsyType, biopsy.biopsySite, biopsy.biopsyLocation,
	treatment.treatmentGiven, treatment.radiotherapyGiven, treatment.name AS treatment,
	treatment.type AS consolidatedTreatmentType, GROUP_CONCAT(drug.type SEPARATOR '|') as concatenatedTreatmentType,
	treatment.mechanism AS consolidatedTreatmentMechanism, GROUP_CONCAT(drug.mechanism SEPARATOR '|') as concatenatedTreatmentMechanism,
	firstMatchedTreatmentResponse.measurementDone AS responseMeasured, firstMatchedTreatmentResponse.response AS firstResponse
	FROM sample
		INNER JOIN patient ON sample.patientId = patient.id
        LEFT JOIN baseline ON patient.id = baseline.patientId
		LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
		LEFT JOIN treatment ON treatment.biopsyId = biopsy.id
		LEFT JOIN drug on drug.treatmentId = treatment.id
		LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId;;