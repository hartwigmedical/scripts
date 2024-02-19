# TODO (KD) Could potentially reuse clinical view once consents are integrated into sample table.

CREATE OR REPLACE VIEW datarequest AS

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
    firstMatchedTreatmentResponse.response AS firstResponse, purity.purity AS tumorPurity, purity.qcStatus AS purpleQC,
    metric.sufficientCoverage AS sufficientCoverage
FROM sample
    INNER JOIN consentsLAMA ON sample.sampleId = consentsLAMA.sampleId AND allowInternalUse = 'true' AND allowExternalUseWithoutCheck = 'true'
    INNER JOIN purity ON sample.sampleId = purity.sampleId AND purity.qcStatus = 'PASS'
    INNER JOIN metric ON sample.sampleId = metric.sampleId AND metric.sufficientCoverage = 1
    INNER JOIN patient ON sample.patientId = patient.id AND blacklisted = 0
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