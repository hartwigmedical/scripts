CREATE OR REPLACE VIEW clinical_purity AS

SELECT
	clinical.*,
	purity AS tumorPurity,
    qcStatus AS purpleQC
FROM
	clinical
	INNER JOIN purity ON purity.sampleId = clinical.sampleId