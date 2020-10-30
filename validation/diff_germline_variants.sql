SELECT
    MIN(pipeline) AS pipeline,
    chromosome,
    position,
    ref,
    alt,
    COUNT(*)
FROM
    (SELECT
        id,
            'Truth' AS pipeline,
            chromosome,
            position,
            ref,
            alt
    FROM
        reference_validation_sets.germlineVariant
    WHERE
        sampleId = 'COLO829v003T'
            AND filter = 'PASS' UNION SELECT
        id,
            'New' AS pipeline,
            chromosome,
            position,
            ref,
            alt
    FROM
        pipeline_v5_validation.germlineVariant
    WHERE
        sampleId = 'COLO829v003T'
            AND filter = 'PASS') AS a
GROUP BY 2 , 3 , 4 , 5
HAVING COUNT(*) != 2;