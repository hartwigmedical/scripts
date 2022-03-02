SELECT
    MIN(pipeline) AS pipeline,
    gene,
    disruptive,
    reportedDisruption,
    COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
            gene,
            disruptive,
            reportedDisruption
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.svBreakend
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline,
            gene,
            disruptive,
            reportedDisruption
    FROM
        VARIABLE_NEW_DB_SCHEMA.svBreakend
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY gene , disruptive , reportedDisruption
HAVING COUNT(*) != 2;