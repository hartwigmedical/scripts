SELECT
    MIN(pipeline) AS pipeline,
    chromosome,
    position,
    ref,
    alt,
    reported,
    COUNT(*),
    MIN(adjustedVaf)
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
        VARIABLE_TRUTH_DB_SCHEMA.somaticVariant
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
        VARIABLE_NEW_DB_SCHEMA.somaticVariant
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID'
            AND filter = 'PASS') AS a
GROUP BY chromosome , position , ref , alt , reported
HAVING COUNT(*) != 2;
