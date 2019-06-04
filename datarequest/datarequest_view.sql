use hmfpatients;
SELECT 
    `sample`.`sampleId` AS `sampleId`,
    `sample`.`arrivalDate` AS `sampleArrivalDate`,
    `purity`.`purity` AS `purplePurity`,
    `purity`.`status` AS `purpleStatus`,
    `purity`.`qcStatus` AS `purpleQC`,
    `baseline`.`primaryTumorLocation` AS `primaryTumorLocation`,
    `baseline`.`cancerSubtype` AS `cancerSubtype`,
    `biopsy`.`biopsySite` AS `biopsySite`,
    `biopsy`.`biopsyLocation` AS `biopsyLocation`,
    `biopsy`.`biopsyType` AS `biopsyType`,
    `biopsy`.`biopsyDate` AS `biopsyDate`,
    `treatment`.`radiotherapyGiven` AS `radiotherapyGiven`,
    `treatment`.`treatmentGiven` AS `treatmentGiven`,
    `treatment`.`name` AS `treatmentName`,
    `treatment`.`type` AS `treatmentType`,
    `treatment`.`mechanism` AS `treatmentMechanism`,
    `treatment`.`startDate` AS `treatmentStartDate`,
    `treatment`.`endDate` AS `treatmentEndDate`,
    `biopsyDrugs`.`biopsyPostDrugs` AS `biopsyPostDrugs`,
    `biopsyDrugs`.`biopsyPostDrugTypes` AS `biopsyPostDrugTypes`,
    `biopsyDrugs`.`biopsyPostDrugMechanisms` AS `biopsyPostDrugMechanisms`,
    `patientDrugs`.`patientPostDrugs` AS `patientPostDrugs`,
    `patientDrugs`.`patientPostDrugTypes` AS `patientPostDrugTypes`,
    `patientDrugs`.`patientPostDrugMechanisms` AS `patientPostDrugMechanisms`,
    `firstMatchedTreatmentResponse`.`responseDate` AS `responseDate`,
    `firstMatchedTreatmentResponse`.`measurementDone` AS `responseMeasured`,
    `firstMatchedTreatmentResponse`.`response` AS `firstResponse`,
    `baseline`.`hospital` AS `hospital`,
    `baseline`.`gender` AS `gender`,
    `baseline`.`birthYear` AS `birthYear`,
    `baseline`.`hasRadiotherapyPreTreatment` AS `hasRadiotherapyPreTreatment`,
    `baseline`.`hasSystemicPreTreatment` AS `hasSystemicPreTreatment`,
    `baseline`.`preTreatments` AS `preTreatments`,
    `baseline`.`preTreatmentsType` AS `preTreatmentsType`,
    `baseline`.`preTreatmentsMechanism` AS `preTreatmentsMechanism`,
    `baseline`.`registrationDate` AS `registrationDate`,
    `baseline`.`informedConsentDate` AS `informedConsentDate`,
    `patient`.`patientIdentifier` AS `patientId`
FROM
    ((((((`sample`
    JOIN `patient` ON ((`sample`.`patientId` = `patient`.`id`)))
    LEFT JOIN `baseline` ON ((`patient`.`id` = `baseline`.`patientId`)))
    LEFT JOIN `biopsy` ON ((`biopsy`.`sampleId` = `sample`.`sampleId`)))
    LEFT JOIN `purity` ON ((`purity`.`sampleId` = `sample`.`sampleId`)))
    LEFT JOIN `treatment` ON ((`treatment`.`biopsyId` = `biopsy`.`id`)))
    LEFT JOIN `firstMatchedTreatmentResponse` ON ((`treatment`.`id` = `firstMatchedTreatmentResponse`.`treatmentId`)))
        LEFT JOIN
    (SELECT 
        treatment.patientId,
            GROUP_CONCAT(drug.name
                SEPARATOR '|') AS patientPostDrugs,
            GROUP_CONCAT(drug.type
                SEPARATOR '|') AS patientPostDrugTypes,
            GROUP_CONCAT(drug.mechanism
                SEPARATOR '|') AS patientPostDrugMechanisms
    FROM
        drug
    LEFT JOIN treatment ON drug.treatmentId = treatment.id
    WHERE
        treatment.treatmentGiven = 'Yes'
    GROUP BY patientId) patientDrugs ON patient.id = patientDrugs.patientId
        LEFT JOIN
    (SELECT 
        treatment.biopsyId,
            GROUP_CONCAT(drug.name
                SEPARATOR '|') AS biopsyPostDrugs,
            GROUP_CONCAT(drug.type
                SEPARATOR '|') AS biopsyPostDrugTypes,
            GROUP_CONCAT(drug.mechanism
                SEPARATOR '|') AS biopsyPostDrugMechanisms
    FROM
        drug
    LEFT JOIN treatment ON drug.treatmentId = treatment.id
    WHERE
        treatment.treatmentGiven = 'Yes'
    GROUP BY biopsyId) biopsyDrugs ON biopsy.id = biopsyDrugs.biopsyId
WHERE
    (sample.sampleId LIKE 'CPCT%'
        OR sample.sampleId LIKE 'DRUP%'
        OR sample.sampleId LIKE 'WIDE%')
    AND informedConsentDate > '2016-04-20'
ORDER BY sampleId
;
