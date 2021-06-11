CREATE OR REPLACE VIEW osiResistance AS (
SELECT *,
    CASE
    WHEN (preTreatmentLine1 LIKE '%Osi%' AND preTreatmentLine1PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine1 > 0 AND daysBiopsyAfterPDLine1 < 90) THEN "INCLUDE"
    WHEN (preTreatmentLine2 LIKE '%Osi%' AND preTreatmentLine2PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine2 > 0 AND daysBiopsyAfterPDLine2 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine3 LIKE '%Osi%' AND preTreatmentLine3PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine3 > 0 AND daysBiopsyAfterPDLine3 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine4 LIKE '%Osi%' AND preTreatmentLine4PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine4 > 0 AND daysBiopsyAfterPDLine4 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine5 LIKE '%Osi%' AND preTreatmentLine5PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine5 > 0 AND daysBiopsyAfterPDLine5 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine6 LIKE '%Osi%' AND preTreatmentLine6PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine6 > 0 AND daysBiopsyAfterPDLine6 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine7 LIKE '%Osi%' AND preTreatmentLine7PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine7 > 0 AND daysBiopsyAfterPDLine7 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine8 LIKE '%Osi%' AND preTreatmentLine8PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine8 > 0 AND daysBiopsyAfterPDLine8 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine9 LIKE '%Osi%' AND preTreatmentLine9PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine9 > 0 AND daysBiopsyAfterPDLine9 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine10 LIKE '%Osi%' AND preTreatmentLine10PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine10 > 0 AND daysBiopsyAfterPDLine10 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine11 LIKE '%Osi%' AND preTreatmentLine11PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine11 > 0 AND daysBiopsyAfterPDLine11 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine12 LIKE '%Osi%' AND preTreatmentLine12PdDate <> '1000-10-10' AND daysBiopsyAfterStartLine12 > 0 AND daysBiopsyAfterPDLine12 < 90) THEN "INCLUDE"
	ELSE "EXCLUDE"
    END AS "include?"
 FROM
	(SELECT a.patientId, n.sampleId, dateBiopsy, driver,
	pretreatmentLine1, pretreatmentLine1startDate, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(dateBiopsy,pretreatmentLine1StartDate) AS daysBiopsyAfterStartLine1, dateDiff(dateBiopsy,pretreatmentLine1PdDate) AS daysBiopsyAfterPDLine1,
	pretreatmentLine2, pretreatmentLine2startDate, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(dateBiopsy,pretreatmentLine2StartDate) AS daysBiopsyAfterStartLine2, dateDiff(dateBiopsy,pretreatmentLine2PdDate) AS daysBiopsyAfterPDLine2,
	pretreatmentLine3, pretreatmentLine3startDate, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(dateBiopsy,pretreatmentLine3StartDate) AS daysBiopsyAfterStartLine3, dateDiff(dateBiopsy,pretreatmentLine3PdDate) AS daysBiopsyAfterPDLine3,
	pretreatmentLine4, pretreatmentLine4startDate, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(dateBiopsy,pretreatmentLine4StartDate) AS daysBiopsyAfterStartLine4, dateDiff(dateBiopsy,pretreatmentLine4PdDate) AS daysBiopsyAfterPDLine4,
	pretreatmentLine5, pretreatmentLine5startDate, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(dateBiopsy,pretreatmentLine5StartDate) AS daysBiopsyAfterStartLine5, dateDiff(dateBiopsy,pretreatmentLine5PdDate) AS daysBiopsyAfterPDLine5,
	pretreatmentLine6, pretreatmentLine6startDate, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(dateBiopsy,pretreatmentLine6StartDate) AS daysBiopsyAfterStartLine6, dateDiff(dateBiopsy,pretreatmentLine6PdDate) AS daysBiopsyAfterPDLine6,
	pretreatmentLine7, pretreatmentLine7startDate, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(dateBiopsy,pretreatmentLine7StartDate) AS daysBiopsyAfterStartLine7, dateDiff(dateBiopsy,pretreatmentLine7PdDate) AS daysBiopsyAfterPDLine7,
	pretreatmentLine8, pretreatmentLine8startDate, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(dateBiopsy,pretreatmentLine8StartDate) AS daysBiopsyAfterStartLine8, dateDiff(dateBiopsy,pretreatmentLine8PdDate) AS daysBiopsyAfterPDLine8,
	pretreatmentLine9, pretreatmentLine9startDate, pretreatmentLine9PdDate, pretreatmentLine9StopReason, dateDiff(dateBiopsy,pretreatmentLine9StartDate) AS daysBiopsyAfterStartLine9, dateDiff(dateBiopsy,pretreatmentLine9PdDate) AS daysBiopsyAfterPDLine9,
	pretreatmentLine10, pretreatmentLine10startDate, pretreatmentLine10PdDate, pretreatmentLine10StopReason, dateDiff(dateBiopsy,pretreatmentLine10StartDate) AS daysBiopsyAfterStartLine10, dateDiff(dateBiopsy,pretreatmentLine10PdDate) AS daysBiopsyAfterPDLine10,
	pretreatmentLine11, pretreatmentLine11startDate, pretreatmentLine11PdDate, pretreatmentLine11StopReason, dateDiff(dateBiopsy,pretreatmentLine11StartDate) AS daysBiopsyAfterStartLine11, dateDiff(dateBiopsy,pretreatmentLine11PdDate) AS daysBiopsyAfterPDLine11,
	pretreatmentLine12, pretreatmentLine12startDate, pretreatmentLine12PdDate, pretreatmentLine12StopReason, dateDiff(dateBiopsy,pretreatmentLine12StartDate) AS daysBiopsyAfterStartLine12, dateDiff(dateBiopsy,pretreatmentLine12PdDate) AS daysBiopsyAfterPDLine12
    FROM
	nsclc.clinical n
    INNER JOIN hmfpatients.purity p ON n.sampleId=p.sampleId
    INNER JOIN hmfpatients.amberPatient a ON n.sampleId=a.sampleId
	WHERE
	(qcStatus NOT LIKE '%WARN_LOW_PURITY%' AND qcStatus <> 'FAIL_NO_TUMOR' AND n.sampleId NOT LIKE 'CORE%')
    AND
    (preTreatmentLine1 LIKE '%Osi%' OR
    preTreatmentLine2 LIKE '%Osi%' OR
    preTreatmentLine3 LIKE '%Osi%' OR
    preTreatmentLine4 LIKE '%Osi%' OR
    preTreatmentLine5 LIKE '%Osi%' OR
    preTreatmentLine6 LIKE '%Osi%' OR
    preTreatmentLine7 LIKE '%Osi%' OR
    preTreatmentLine8 LIKE '%Osi%' OR
    preTreatmentLine9 LIKE '%Osi%' OR
    preTreatmentLine10 LIKE '%Osi%' OR
    preTreatmentLine11 LIKE '%Osi%' OR
    preTreatmentLine12 LIKE '%Osi%'))
		AS potentialSelection
);