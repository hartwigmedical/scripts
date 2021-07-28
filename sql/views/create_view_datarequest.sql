CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical.*, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical
	INNER JOIN metric ON metric.sampleId = clinical.sampleId
    INNER JOIN snpcheck ON snpcheck.sampleId = clinical.sampleId
WHERE
	blacklisted = 0 AND purpleQC = 'PASS' AND metric.sufficientCoverage = 1 AND snpcheck.isPass = 1 AND
	(clinical.sampleId LIKE '%CPCT%' OR clinical.sampleId LIKE '%WIDE%' OR clinical.sampleId LIKE '%DRUP%' OR clinical.sampleId LIKE '%ACTN%')