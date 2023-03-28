CREATE OR REPLACE VIEW clinical AS

SELECT amberAnonymous.hmfSampleId, left(amberAnonymous.hmfSampleId, 9) as hmfPatientId,
    sample.sampleId, patient.patientIdentifier AS patientId, not(isnull(rnaStatistics.sampleId)) as hasRNA, setName,
    sample.arrivalDate as sampleArrivalDate, patient.blacklisted, baseline.registrationDate, baseline.informedConsentDate,
    baseline.inHMFDatabase, baseline.outsideEU, baseline.deathDate,
    baseline.primaryTumorLocation, baseline.primaryTumorSubLocation, baseline.primaryTumorType, baseline.primaryTumorSubType,
    baseline.primaryTumorExtraDetails, doidView.doids,
    baseline.hospital, baseline.gender, baseline.birthYear,
    baseline.hasSystemicPreTreatment, baseline.hasRadiotherapyPreTreatment,
    baseline.preTreatments, baseline.preTreatmentsType, baseline.preTreatmentsMechanism,
    biopsy.biopsyDate, biopsy.biopsyType, biopsy.biopsySite, biopsy.biopsyLocation,
	treatment.treatmentGiven, treatment.radiotherapyGiven, treatment.startDate AS treatmentStartDate, treatment.endDate AS treatmentEndDate,
    treatment.name AS treatment, treatment.type AS consolidatedTreatmentType, biopsyDrugs.treatmentType AS concatenatedTreatmentType,
	treatment.mechanism AS consolidatedTreatmentMechanism, biopsyDrugs.treatmentMechanism as concatenatedTreatmentMechanism,
    firstMatchedTreatmentResponse.measurementDone AS responseMeasured, firstMatchedTreatmentResponse.responseDate,
    firstMatchedTreatmentResponse.response AS firstResponse
FROM sample
    INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN amberAnonymous on sample.sampleId = amberAnonymous.sampleId AND deleted = 0
    LEFT JOIN rnaStatistics on sample.sampleId = rnaStatistics.sampleId
    LEFT JOIN baseline ON patient.id = baseline.patientId
    LEFT JOIN (SELECT patientId, group_concat(doid separator ",") AS doids FROM doidNode GROUP BY 1) AS doidView
        ON patient.id = doidView.patientId
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment ON treatment.biopsyId = biopsy.id
    LEFT JOIN
    (SELECT treatment.biopsyId,
        GROUP_CONCAT(drug.type SEPARATOR '/') AS treatmentType,
        GROUP_CONCAT(drug.mechanism SEPARATOR '/') AS treatmentMechanism
        FROM drug INNER JOIN treatment ON drug.treatmentId = treatment.id GROUP BY biopsyId)
        biopsyDrugs ON biopsy.id = biopsyDrugs.biopsyId
    LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId
    INNER JOIN consentsLAMA ON SUBSTRING_INDEX(SUBSTRING_INDEX(sample.setName, '_', - 2), '_', 1)=consentsLAMA.barcode
 WHERE consentsLAMA.allowInternalUse = 'true'
