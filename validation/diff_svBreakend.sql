SELECT
    MIN(pipeline) AS Source, gene, disruptive, reportedDisruption, COUNT(*)
FROM
    (SELECT
        id, 'OnlyInTruth' AS pipeline, gene, disruptive, reportedDisruption
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.svBreakend
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
    UNION
    SELECT
        id, 'OnlyInNew' AS pipeline, gene, disruptive, reportedDisruption
    FROM
        VARIABLE_NEW_DB_SCHEMA.svBreakend
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID')
    AS a
GROUP BY 2 , 3 , 4
HAVING COUNT(*) = 1
ORDER BY Gene;