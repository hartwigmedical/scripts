CREATE OR REPLACE VIEW clinical AS

SELECT sample.sampleId, sample.arrivalDate as sampleArrivalDate, patient.patientIdentifier AS patientId,
	baseline.registrationDate, baseline.informedConsentDate, baseline.hospital, baseline.gender,
    baseline.birthYear, baseline.hasSystemicPreTreatment, baseline.hasRadiotherapyPreTreatment, baseline.preTreatments, baseline.preTreatmentsType,
    baseline.preTreatmentsMechanism, biopsy.biopsyDate,
	treatment.startDate AS treatmentStartDate, treatment.endDate AS treatmentEndDate, firstMatchedTreatmentResponse.responseDate, baseline.deathDate,
	baseline.primaryTumorLocation, baseline.cancerSubtype, biopsy.biopsyType, biopsy.biopsySite, biopsy.biopsyLocation,
	treatment.treatmentGiven, treatment.radiotherapyGiven, treatment.name AS treatment,
	treatment.type AS consolidatedTreatmentType, biopsyDrugs.treatmentType AS concatenatedTreatmentType,
	treatment.mechanism AS consolidatedTreatmentMechanism, biopsyDrugs.treatmentMechanism as concatenatedTreatmentMechanism,
	firstMatchedTreatmentResponse.measurementDone AS responseMeasured, firstMatchedTreatmentResponse.response AS firstResponse
FROM sample
    INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN baseline ON patient.id = baseline.patientId
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment ON treatment.biopsyId = biopsy.id
    LEFT JOIN
    (SELECT treatment.biopsyId,
        GROUP_CONCAT(drug.type SEPARATOR '/') AS treatmentType,
        GROUP_CONCAT(drug.mechanism SEPARATOR '/') AS treatmentMechanism
        FROM drug INNER JOIN treatment ON drug.treatmentId = treatment.id GROUP BY biopsyId)
        biopsyDrugs ON biopsy.id = biopsyDrugs.biopsyId
    LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId;