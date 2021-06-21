CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical.*, purity AS tumorPurity, qcStatus AS purpleQC, version AS purpleVersion, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical
	INNER JOIN purity ON purity.sampleId = clinical.sampleId
	INNER JOIN metric ON metric.sampleId = clinical.sampleId
    INNER JOIN snpcheck ON snpcheck.sampleId = clinical.sampleId
WHERE
	blacklisted = 0 AND purity.qcStatus = 'PASS' AND metric.sufficientCoverage = 1 AND snpcheck.isPass = 1 AND
	(clinical.sampleId LIKE '%CPCT%' OR clinical.sampleId LIKE '%WIDE%' OR clinical.sampleId LIKE '%DRUP%' OR clinical.sampleId LIKE '%ACTN%');
