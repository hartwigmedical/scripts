use hmfpatients;
CREATE OR REPLACE VIEW datarequest_CB01 AS
SELECT
	clinical_purity.*, metric.sufficientCoverage AS sufficientCoverage
FROM
	clinical_purity
	INNER JOIN metric ON metric.sampleId = clinical_purity.sampleId
    INNER JOIN consentsLAMA ON clinical_purity.sampleId=consentsLAMA.sampleId
WHERE
	blacklisted = 0
	AND purpleQC = 'PASS'
	AND metric.sufficientCoverage = 1
	AND consentsLAMA.allowExternalUseWithCheck = 'true'
	;