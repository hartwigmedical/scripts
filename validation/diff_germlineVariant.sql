SELECT
    MIN(pipeline) AS pipeline,
    chromosome,
    position,
    ref,
    alt,
    adjustedVaf,
    reported,
    COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
            chromosome,
            position,
            ref,
            alt,
            ROUND(adjustedVaf,2) AS adjustedVaf,
            reported
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.germlineVariant
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
            AND filter = 'PASS' UNION SELECT
        'OnlyInNew' AS pipeline,
            chromosome,
            position,
            ref,
            alt,
            ROUND(adjustedVaf,2) AS adjustedVaf,
            reported
    FROM
        VARIABLE_NEW_DB_SCHEMA.germlineVariant
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID'
            AND filter = 'PASS') AS a
GROUP BY chromosome , position , ref , alt , adjustedVaf, reported
HAVING COUNT(*) != 2;