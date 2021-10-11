CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical.*, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical
	INNER JOIN metric ON metric.sampleId = clinical.sampleId
    INNER JOIN snpcheck ON snpcheck.sampleId = clinical.sampleId
    INNER JOIN purity on purity.sampleId = clinical.sampleId
WHERE
	blacklisted = 0 AND purity.qcStatus = 'PASS' AND metric.sufficientCoverage = 1 AND snpcheck.isPass = 1 AND inHMFDatabase =1 AND outsideEU =1 AND
	(clinical.sampleId LIKE '%CPCT%' OR clinical.sampleId LIKE '%WIDE%' OR clinical.sampleId LIKE '%DRUP%' OR clinical.sampleId LIKE '%ACTN%');