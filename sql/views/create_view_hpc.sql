#for patients with multiples samples, this view (as opposed to the datarequest view) contains only the sample with the highest purity
CREATE OR REPLACE VIEW hpc AS

SELECT SUBSTRING_INDEX(group_concat(purity.sampleId ORDER BY purity DESC, id ASC SEPARATOR ','), ",", 1) AS sampleId
FROM purity
    INNER JOIN amberPatient ON purity.sampleId = amberPatient.sampleId
    INNER JOIN metric ON purity.sampleId = metric.sampleId
WHERE qcStatus = 'PASS' and sufficientCoverage = 1
GROUP BY patientId;