CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical_purity.*, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical_purity
	INNER JOIN metric ON metric.sampleId = clinical_purity.sampleId
    INNER JOIN snpcheck ON snpcheck.sampleId = clinical_purity.sampleId
    INNER JOIN consentsLAMA ON SUBSTRING_INDEX(SUBSTRING_INDEX(clinical_purity.setName, '_', - 2), '_', 1)=consentsLAMA.barcode
WHERE
	blacklisted = 0
	AND purpleQC = 'PASS'
	AND metric.sufficientCoverage = 1
	AND snpcheck.isPass = 1
	AND inHMFDatabase =1
	AND outsideEU =1
	AND consentsLAMA.allowExternalUseWithoutCheck = 'true'
	;