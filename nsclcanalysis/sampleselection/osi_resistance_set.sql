SELECT *,
	CASE
    WHEN (preTreatmentLine1 like '%Osi%' and preTreatmentLine1PdDate <> '1000-10-10' and daysBiopsyAfterStartLine1 > 0 and daysBiopsyAfterPDLine1 < 90) THEN "INCLUDE"
    WHEN (preTreatmentLine2 like '%Osi%' and preTreatmentLine2PdDate <> '1000-10-10' and daysBiopsyAfterStartLine2 > 0 and daysBiopsyAfterPDLine2 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine3 like '%Osi%' and preTreatmentLine3PdDate <> '1000-10-10' and daysBiopsyAfterStartLine3 > 0 and daysBiopsyAfterPDLine3 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine4 like '%Osi%' and preTreatmentLine4PdDate <> '1000-10-10' and daysBiopsyAfterStartLine4 > 0 and daysBiopsyAfterPDLine4 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine5 like '%Osi%' and preTreatmentLine5PdDate <> '1000-10-10' and daysBiopsyAfterStartLine5 > 0 and daysBiopsyAfterPDLine5 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine6 like '%Osi%' and preTreatmentLine6PdDate <> '1000-10-10' and daysBiopsyAfterStartLine6 > 0 and daysBiopsyAfterPDLine6 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine7 like '%Osi%' and preTreatmentLine7PdDate <> '1000-10-10' and daysBiopsyAfterStartLine7 > 0 and daysBiopsyAfterPDLine7 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine8 like '%Osi%' and preTreatmentLine8PdDate <> '1000-10-10' and daysBiopsyAfterStartLine8 > 0 and daysBiopsyAfterPDLine8 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine9 like '%Osi%' and preTreatmentLine9PdDate <> '1000-10-10' and daysBiopsyAfterStartLine9 > 0 and daysBiopsyAfterPDLine9 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine10 like '%Osi%' and preTreatmentLine10PdDate <> '1000-10-10' and daysBiopsyAfterStartLine10 > 0 and daysBiopsyAfterPDLine10 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine11 like '%Osi%' and preTreatmentLine11PdDate <> '1000-10-10' and daysBiopsyAfterStartLine11 > 0 and daysBiopsyAfterPDLine11 < 90) THEN "INCLUDE"
	WHEN (preTreatmentLine12 like '%Osi%' and preTreatmentLine12PdDate <> '1000-10-10' and daysBiopsyAfterStartLine12 > 0 and daysBiopsyAfterPDLine12 < 90) THEN "INCLUDE"
	ELSE "EXCLUDE"
    END AS "include?"
 FROM
	(SELECT a.patientId, n.sampleId, dateBiopsy, driver,
	pretreatmentLine1, pretreatmentLine1startDate, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(dateBiopsy,pretreatmentLine1StartDate) as daysBiopsyAfterStartLine1, dateDiff(dateBiopsy,pretreatmentLine1PdDate) as daysBiopsyAfterPDLine1,
	pretreatmentLine2, pretreatmentLine2startDate, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(dateBiopsy,pretreatmentLine2StartDate) as daysBiopsyAfterStartLine2, dateDiff(dateBiopsy,pretreatmentLine2PdDate) as daysBiopsyAfterPDLine2,
	pretreatmentLine3, pretreatmentLine3startDate, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(dateBiopsy,pretreatmentLine3StartDate) as daysBiopsyAfterStartLine3, dateDiff(dateBiopsy,pretreatmentLine3PdDate) as daysBiopsyAfterPDLine3,
	pretreatmentLine4, pretreatmentLine4startDate, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(dateBiopsy,pretreatmentLine4StartDate) as daysBiopsyAfterStartLine4, dateDiff(dateBiopsy,pretreatmentLine4PdDate) as daysBiopsyAfterPDLine4,
	pretreatmentLine5, pretreatmentLine5startDate, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(dateBiopsy,pretreatmentLine5StartDate) as daysBiopsyAfterStartLine5, dateDiff(dateBiopsy,pretreatmentLine5PdDate) as daysBiopsyAfterPDLine5,
	pretreatmentLine6, pretreatmentLine6startDate, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(dateBiopsy,pretreatmentLine6StartDate) as daysBiopsyAfterStartLine6, dateDiff(dateBiopsy,pretreatmentLine6PdDate) as daysBiopsyAfterPDLine6,
	pretreatmentLine7, pretreatmentLine7startDate, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(dateBiopsy,pretreatmentLine7StartDate) as daysBiopsyAfterStartLine7, dateDiff(dateBiopsy,pretreatmentLine7PdDate) as daysBiopsyAfterPDLine7,
	pretreatmentLine8, pretreatmentLine8startDate, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(dateBiopsy,pretreatmentLine8StartDate) as daysBiopsyAfterStartLine8, dateDiff(dateBiopsy,pretreatmentLine8PdDate) as daysBiopsyAfterPDLine8,
	pretreatmentLine9, pretreatmentLine9startDate, pretreatmentLine9PdDate, pretreatmentLine9StopReason, dateDiff(dateBiopsy,pretreatmentLine9StartDate) as daysBiopsyAfterStartLine9, dateDiff(dateBiopsy,pretreatmentLine9PdDate) as daysBiopsyAfterPDLine9,
	pretreatmentLine10, pretreatmentLine10startDate, pretreatmentLine10PdDate, pretreatmentLine10StopReason, dateDiff(dateBiopsy,pretreatmentLine10StartDate) as daysBiopsyAfterStartLine10, dateDiff(dateBiopsy,pretreatmentLine10PdDate) as daysBiopsyAfterPDLine10,
	pretreatmentLine11, pretreatmentLine11startDate, pretreatmentLine11PdDate, pretreatmentLine11StopReason, dateDiff(dateBiopsy,pretreatmentLine11StartDate) as daysBiopsyAfterStartLine11, dateDiff(dateBiopsy,pretreatmentLine11PdDate) as daysBiopsyAfterPDLine11,
	pretreatmentLine12, pretreatmentLine12startDate, pretreatmentLine12PdDate, pretreatmentLine12StopReason, dateDiff(dateBiopsy,pretreatmentLine12StartDate) as daysBiopsyAfterStartLine12, dateDiff(dateBiopsy,pretreatmentLine12PdDate) as daysBiopsyAfterPDLine12
    FROM
	nsclc.clinical n
    inner join hmfpatients.purity p on n.sampleId=p.sampleId
    inner join hmfpatients.amberPatient a on n.sampleId=a.sampleId
	WHERE
	(qcStatus NOT LIKE '%WARN_LOW_PURITY%' AND qcStatus <> 'FAIL_NO_TUMOR' AND n.sampleId not like 'CORE%')
    AND
    (preTreatmentLine1 like '%Osi%' OR
    preTreatmentLine2 like '%Osi%' OR
    preTreatmentLine3 like '%Osi%' OR
    preTreatmentLine4 like '%Osi%' OR
    preTreatmentLine5 like '%Osi%' OR
    preTreatmentLine6 like '%Osi%' OR
    preTreatmentLine7 like '%Osi%' OR
    preTreatmentLine8 like '%Osi%' OR
    preTreatmentLine9 like '%Osi%' OR
    preTreatmentLine10 like '%Osi%' OR
    preTreatmentLine11 like '%Osi%' OR
    preTreatmentLine12 like '%Osi%'))
		as potentialSelection;