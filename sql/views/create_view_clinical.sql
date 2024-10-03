CREATE OR REPLACE VIEW clinical AS
SELECT amberAnonymous.hmfSampleId,
       left(amberAnonymous.hmfSampleId, 9)           as hmfPatientId,
       sample.sampleId,
       patient.donorId                               AS patientId,
       not (isnull(rnaStatistics.sampleId))          as hasRNA,
       specimen.arrivalDate                          AS sampleArrivalDate,
       patient.blacklisted,
       baseline.registrationDate,
       baseline.informedConsentDate,
       baseline.inHMFDatabase,
       baseline.outsideEU,
       baseline.deathDate,
       baseline.primaryTumorLocation,
       baseline.primaryTumorType,
       baseline.primaryTumorExtraDetails,
       doidView.doids,
       baseline.hospital,
       baseline.gender,
       baseline.birthYear,
       baseline.hasSystemicPreTreatment,
       baseline.hasRadiotherapyPreTreatment,
       baseline.preTreatments,
       baseline.preTreatmentsType,
       baseline.preTreatmentsMechanism,
       biopsy.biopsyDate,
       biopsy.biopsyType,
       biopsy.biopsySite,
       biopsy.biopsyLocation,
       treatment.treatmentGiven,
       treatment.radiotherapyGiven,
       treatment.startDate                           AS treatmentStartDate,
       treatment.endDate                             AS treatmentEndDate,
       treatment.name                                AS treatment,
       treatment.type                                AS consolidatedTreatmentType,
       biopsyDrugs.treatmentType                     AS concatenatedTreatmentType,
       treatment.mechanism                           AS consolidatedTreatmentMechanism,
       biopsyDrugs.treatmentMechanism                as concatenatedTreatmentMechanism,
       firstMatchedTreatmentResponse.measurementDone AS responseMeasured,
       firstMatchedTreatmentResponse.responseDate,
       firstMatchedTreatmentResponse.response        AS firstResponse
FROM sample
         INNER JOIN specimen ON specimen.specimenId = sample.specimenId
         INNER JOIN patient ON specimen.patientId = patient.id
         LEFT JOIN amberAnonymous on sample.sampleId = amberAnonymous.sampleId AND deleted = 0
         LEFT JOIN rnaStatistics on sample.sampleId = rnaStatistics.sampleId
         LEFT JOIN baseline ON patient.id = baseline.patientId
         LEFT JOIN (SELECT patientId, group_concat(doid separator ',') AS doids FROM doidNode GROUP BY 1) AS doidView
                   ON patient.id = doidView.patientId
         LEFT JOIN biopsy ON biopsy.specimenId = sample.specimenId
         LEFT JOIN first_treatment_after_biopsy as treatment ON treatment.patientId = specimen.patientId
         LEFT JOIN
     (SELECT treatment.biopsyId,
             GROUP_CONCAT(drug.type SEPARATOR '/')      AS treatmentType,
             GROUP_CONCAT(drug.mechanism SEPARATOR '/') AS treatmentMechanism
      FROM drug
               INNER JOIN treatment ON drug.treatmentId = treatment.id
      GROUP BY biopsyId) biopsyDrugs ON biopsy.id = biopsyDrugs.biopsyId
         LEFT JOIN firstMatchedTreatmentResponse ON treatment.id = firstMatchedTreatmentResponse.treatmentId
WHERE allowInternalUse = 1;