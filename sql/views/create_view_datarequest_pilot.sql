CREATE OR REPLACE VIEW datarequest AS

SELECT
	clinical.*, purity AS tumorPurity, qcStatus AS purpleQC, status AS purpleStatus, version AS purpleVersion, metric.metricQC AS metricQC
FROM
	clinical
	INNER JOIN purity ON purity.sampleId = clinical.sampleId
	INNER JOIN metric ON metricQC.sampleId = clinical.sampleId
WHERE
	blacklisted = 0 AND qcStatus = 'PASS' AND metricQC = 1 AND status <> 'NO_TUMOR' AND purity > 0.195 AND
	((clinical.sampleId LIKE '%CPCT%' AND informedConsentDate > '2016-04-20')
	    OR clinical.sampleId LIKE '%WIDE%'
	    OR clinical.sampleId LIKE '%DRUP%');
