CREATE OR REPLACE VIEW clinical_purity AS

SELECT clinical.*,
       lower(purity.gender) as gender,
       purity               AS tumorPurity,
       qcStatus             AS purpleQC
FROM clinical
         INNER JOIN purity ON purity.sampleId = clinical.sampleId