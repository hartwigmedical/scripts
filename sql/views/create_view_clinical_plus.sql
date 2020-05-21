CREATE OR REPLACE VIEW clinicalPlus AS

SELECT
	clinical.*, purity AS tumorPurity, qcStatus AS purpleQC, status AS purpleStatus, version AS purpleVersion
FROM
	clinical INNER JOIN purity ON purity.sampleId = clinical.sampleId;