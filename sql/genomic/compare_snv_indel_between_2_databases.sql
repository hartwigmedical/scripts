SELECT min(pipeline) as pipeline, chromosome, position, ref, alt, count(*) FROM
(SELECT id, 'Truth' AS pipeline, chromosome, position, ref, alt FROM reference_validation_sets.somaticVariant WHERE sampleId = 'XXX' AND filter = "PASS"
UNION
SELECT id, 'New' AS pipeline, chromosome, position, ref, alt FROM pipeline_v5_validation.somaticVariant WHERE sampleId = 'XXX' AND filter = "PASS") AS a
GROUP BY 2,3,4,5 HAVING count(*) = 1;