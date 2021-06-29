CREATE OR REPLACE VIEW nsclc.osiResistance AS (
SELECT *,
    CASE
    WHEN driver NOT LIKE '%EGFR%'
		THEN "EXCLUDE_0_NO_EGFR_DRIVER"
    WHEN (preTreatmentLine1 NOT LIKE '%Osi%' AND preTreatmentLine2 NOT LIKE '%Osi%' AND preTreatmentLine3 NOT LIKE '%Osi%' AND preTreatmentLine4 NOT LIKE '%Osi%' AND preTreatmentLine5 NOT LIKE '%Osi%' AND preTreatmentLine6 NOT LIKE '%Osi%' AND preTreatmentLine7 NOT LIKE '%Osi%' AND preTreatmentLine8 NOT LIKE '%Osi%' AND preTreatmentLine9 NOT LIKE '%Osi%' AND preTreatmentLine10 NOT LIKE '%Osi%' AND preTreatmentLine11 NOT LIKE '%Osi%' AND preTreatmentLine12 NOT LIKE '%Osi%')
		THEN "EXCLUDE_1_NO_OSI_BEFORE_BIOPSY"
    WHEN sampleId LIKE 'CORE%'
		THEN "EXCLUDE_2_COREDB_SAMPLE"
	WHEN (qcStatus LIKE '%WARN_LOW_PURITY%' OR qcStatus = 'FAIL_NO_TUMOR')
		THEN "EXCLUDE_3_QC_CHECK"
    WHEN (preTreatmentLine1 LIKE '%Osi%' AND preTreatmentLine1PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine1 < 100) THEN "INCLUDE"
    WHEN (preTreatmentLine2 LIKE '%Osi%' AND preTreatmentLine2PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine2 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine3 LIKE '%Osi%' AND preTreatmentLine3PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine3 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine4 LIKE '%Osi%' AND preTreatmentLine4PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine4 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine5 LIKE '%Osi%' AND preTreatmentLine5PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine5 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine6 LIKE '%Osi%' AND preTreatmentLine6PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine6 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine7 LIKE '%Osi%' AND preTreatmentLine7PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine7 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine8 LIKE '%Osi%' AND preTreatmentLine8PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine8 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine9 LIKE '%Osi%' AND preTreatmentLine9PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine9 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine10 LIKE '%Osi%' AND preTreatmentLine10PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine10 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine11 LIKE '%Osi%' AND preTreatmentLine11PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine11 < 100) THEN "INCLUDE"
	WHEN (preTreatmentLine12 LIKE '%Osi%' AND preTreatmentLine12PdDate <> '1000-10-10' AND daysBiopsyAfterPDLine12 < 100) THEN "INCLUDE"
	ELSE "EXCLUDE_4_NO_PD_OR_TIME_RULES"
    END AS "include"
 FROM
	(SELECT n.sampleId, qcStatus, note, driver, dateDiagnosis, dateStageIV, dateDeathOrLastContact, alive, dateBiopsy,
	pretreatmentLine1, pretreatmentLine1startDate, pretreatmentLine1endDate, pretreatmentLine1PdDate, pretreatmentLine1StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine1StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine1StartDate)) AS daysBiopsyAfterStartLine1, IF(dateDiff(dateBiopsy,pretreatmentLine1PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine1PdDate)) AS daysBiopsyAfterPDLine1,
	pretreatmentLine2, pretreatmentLine2startDate, pretreatmentLine2endDate, pretreatmentLine2PdDate, pretreatmentLine2StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine2StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine2StartDate)) AS daysBiopsyAfterStartLine2, IF(dateDiff(dateBiopsy,pretreatmentLine2PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine2PdDate)) AS daysBiopsyAfterPDLine2,
	pretreatmentLine3, pretreatmentLine3startDate, pretreatmentLine3endDate, pretreatmentLine3PdDate, pretreatmentLine3StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine3StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine3StartDate)) AS daysBiopsyAfterStartLine3, IF(dateDiff(dateBiopsy,pretreatmentLine3PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine3PdDate)) AS daysBiopsyAfterPDLine3,
	pretreatmentLine4, pretreatmentLine4startDate, pretreatmentLine4endDate, pretreatmentLine4PdDate, pretreatmentLine4StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine4StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine4StartDate)) AS daysBiopsyAfterStartLine4, IF(dateDiff(dateBiopsy,pretreatmentLine4PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine4PdDate)) AS daysBiopsyAfterPDLine4,
	pretreatmentLine5, pretreatmentLine5startDate, pretreatmentLine5endDate, pretreatmentLine5PdDate, pretreatmentLine5StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine5StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine5StartDate)) AS daysBiopsyAfterStartLine5, IF(dateDiff(dateBiopsy,pretreatmentLine5PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine5PdDate)) AS daysBiopsyAfterPDLine5,
	pretreatmentLine6, pretreatmentLine6startDate, pretreatmentLine6endDate, pretreatmentLine6PdDate, pretreatmentLine6StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine6StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine6StartDate)) AS daysBiopsyAfterStartLine6, IF(dateDiff(dateBiopsy,pretreatmentLine6PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine6PdDate)) AS daysBiopsyAfterPDLine6,
	pretreatmentLine7, pretreatmentLine7startDate, pretreatmentLine7endDate, pretreatmentLine7PdDate, pretreatmentLine7StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine7StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine7StartDate)) AS daysBiopsyAfterStartLine7, IF(dateDiff(dateBiopsy,pretreatmentLine7PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine7PdDate)) AS daysBiopsyAfterPDLine7,
	pretreatmentLine8, pretreatmentLine8startDate, pretreatmentLine8endDate, pretreatmentLine8PdDate, pretreatmentLine8StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine8StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine8StartDate)) AS daysBiopsyAfterStartLine8, IF(dateDiff(dateBiopsy,pretreatmentLine8PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine8PdDate)) AS daysBiopsyAfterPDLine8,
	pretreatmentLine9, pretreatmentLine9startDate, pretreatmentLine9endDate, pretreatmentLine9PdDate, pretreatmentLine9StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine9StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine9StartDate)) AS daysBiopsyAfterStartLine9, IF(dateDiff(dateBiopsy,pretreatmentLine9PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine9PdDate)) AS daysBiopsyAfterPDLine9,
	pretreatmentLine10, pretreatmentLine10startDate, pretreatmentLine10endDate, pretreatmentLine10PdDate, pretreatmentLine10StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine10StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine10StartDate)) AS daysBiopsyAfterStartLine10, IF(dateDiff(dateBiopsy,pretreatmentLine10PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine10PdDate)) AS daysBiopsyAfterPDLine10,
	pretreatmentLine11, pretreatmentLine11startDate, pretreatmentLine11endDate, pretreatmentLine11PdDate, pretreatmentLine11StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine11StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine11StartDate)) AS daysBiopsyAfterStartLine11, IF(dateDiff(dateBiopsy,pretreatmentLine11PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine11PdDate)) AS daysBiopsyAfterPDLine11,
	pretreatmentLine12, pretreatmentLine12startDate, pretreatmentLine12endDate, pretreatmentLine12PdDate, pretreatmentLine12StopReason, IF(dateDiff(dateBiopsy,pretreatmentLine12StartDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine12StartDate)) AS daysBiopsyAfterStartLine12, IF(dateDiff(dateBiopsy,pretreatmentLine12PdDate)>300000,"NA",dateDiff(dateBiopsy,pretreatmentLine12PdDate)) AS daysBiopsyAfterPDLine12
    FROM
	nsclc.clinical n
    INNER JOIN hmfpatients.purity p ON n.sampleId=p.sampleId
    INNER JOIN hmfpatients.amberPatient a ON n.sampleId=a.sampleId)
		AS potentialSelection
);