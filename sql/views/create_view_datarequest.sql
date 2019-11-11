CREATE OR REPLACE VIEW datarequest AS

SELECT 
	clinical.*, purity AS purplePurity, qcStatus AS purpleQC, status AS purpleStatus, version AS purpleVersion 
FROM 
	clinical INNER JOIN purity ON purity.sampleId = clinical.sampleId
WHERE 
	qcStatus = 'PASS' AND status <> 'NO_TUMOR' AND purity > 0.195 AND 
	((clinical.sampleId LIKE '%CPCT%' AND informedConsentDate > '2016-04-20')
	    OR clinical.sampleId LIKE '%WIDE%'
	    OR clinical.sampleId LIKE '%DRUP%');
