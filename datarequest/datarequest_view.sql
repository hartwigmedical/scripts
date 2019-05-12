/* Descr: Creates a starting point for data requests. */
use hmfpatients;
SELECT 
    `sample`.`sampleId` AS `sampleId`,
    `sample`.`arrivalDate` AS `sampleArrivalDate`,
    `purity`.`purity` AS `purplePurity`,
    `purity`.`status` AS `purpleStatus`,
    `purity`.`qcStatus` AS `purpleQC`,
    `patientTreatment`.`patientPostTreatments` AS `patientPostTreatments`,
    `patientTreatment`.`patientPostTreatmentTypes` AS `patientPostTreatmentTypes`,
    `patientTreatment`.`patientPostTreatmentMechanisms` AS `patientPostTreatmentMechanisms`,
    `biopsyTreatment`.`biopsyPostTreatments` AS `biopsyPostTreatments`,
    `biopsyTreatment`.`biopsyPostTreatmentTypes` AS `biopsyPostTreatmentTypes`,
    `biopsyTreatment`.`biopsyPostTreatmentMechanisms` AS `biopsyPostTreatmentMechanisms`,
    `patient`.`patientIdentifier` AS `patientId`,
    `baseline`.`registrationDate` AS `registrationDate`,
    `baseline`.`informedConsentDate` AS `informedConsentDate`,
    `baseline`.`hospital` AS `hospital`,
    `baseline`.`gender` AS `gender`,
    `baseline`.`birthYear` AS `birthYear`,
    `baseline`.`hasSystemicPreTreatment` AS `hasSystemicPreTreatment`,
    `baseline`.`hasRadiotherapyPreTreatment` AS `hasRadiotherapyPreTreatment`,
    `baseline`.`preTreatments` AS `preTreatments`,
    `baseline`.`preTreatmentsType` AS `preTreatmentsType`,
    `baseline`.`preTreatmentsMechanism` AS `preTreatmentsMechanism`,
    `biopsy`.`biopsyDate` AS `biopsyDate`,
    `treatment`.`startDate` AS `treatmentStartDate`,
    `treatment`.`endDate` AS `treatmentEndDate`,
    `firstMatchedTreatmentResponse`.`responseDate` AS `responseDate`,
    `baseline`.`deathDate` AS `deathDate`,
    `baseline`.`primaryTumorLocation` AS `primaryTumorLocation`,
    `baseline`.`cancerSubtype` AS `cancerSubtype`,
    `biopsy`.`biopsyType` AS `biopsyType`,
    `biopsy`.`biopsySite` AS `biopsySite`,
    `biopsy`.`biopsyLocation` AS `biopsyLocation`,
    `treatment`.`treatmentGiven` AS `treatmentGiven`,
    `treatment`.`radiotherapyGiven` AS `radiotherapyGiven`,
    `treatment`.`name` AS `treatment`,
    `treatment`.`type` AS `treatmentType`,
    `treatment`.`mechanism` AS `treatmentMechanism`,
    `firstMatchedTreatmentResponse`.`measurementDone` AS `responseMeasured`,
    `firstMatchedTreatmentResponse`.`response` AS `firstResponse`
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
        /* List all treatment info for patient (only explicitly given treatments). */
        patientId,
            GROUP_CONCAT(DISTINCT (treatment.name)
                SEPARATOR '/') AS patientPostTreatments,
            GROUP_CONCAT(DISTINCT (treatment.type)
                SEPARATOR '/') AS patientPostTreatmentTypes,
            GROUP_CONCAT(DISTINCT (treatment.mechanism)
                SEPARATOR '/') AS patientPostTreatmentMechanisms
    FROM
        treatment
    WHERE
        treatment.treatmentGiven = 'Yes'
    GROUP BY patientId) patientTreatment ON patient.id = patientTreatment.patientId
        LEFT JOIN
    (SELECT 
		/* List all treatment info for biopsy (only explicitly given treatments). */
        biopsyId,
            GROUP_CONCAT(DISTINCT (treatment.name)
                SEPARATOR '/') AS biopsyPostTreatments,
            GROUP_CONCAT(DISTINCT (treatment.type)
                SEPARATOR '/') AS biopsyPostTreatmentTypes,
            GROUP_CONCAT(DISTINCT (treatment.mechanism)
                SEPARATOR '/') AS biopsyPostTreatmentMechanisms
    FROM
        treatment
    WHERE
        treatment.treatmentGiven = 'Yes'
    GROUP BY biopsyId) biopsyTreatment ON biopsy.id = biopsyTreatment.biopsyId
WHERE
    (sample.sampleId LIKE 'CPCT%'
        OR sample.sampleId LIKE 'DRUP%'
        OR sample.sampleId LIKE 'WIDE%')
        AND informedConsentDate > '2016-04-20'
ORDER BY sampleId
;

