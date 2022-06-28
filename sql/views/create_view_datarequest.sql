CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical_purity.*, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical_purity
	INNER JOIN metric ON metric.sampleId = clinical_purity.sampleId
    INNER JOIN snpcheck ON snpcheck.sampleId = clinical_purity.sampleId
WHERE
	blacklisted = 0 AND purpleQC = 'PASS' AND metric.sufficientCoverage = 1 AND snpcheck.isPass = 1 AND inHMFDatabase =1 AND outsideEU =1 AND
	(clinical_purity.sampleId LIKE '%CPCT%' OR clinical_purity.sampleId LIKE '%ACTN%' OR clinical_purity.sampleId LIKE '%WIDE%');