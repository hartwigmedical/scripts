SELECT
    MIN(pipeline) AS pipeline, chromosome, position, ref, alt, reported, COUNT(*)
FROM
        (SELECT
            id, 'OnlyInTruth' AS pipeline, chromosome, position, ref, alt, reported
        FROM
            VARIABLE_TRUTH_DB_SCHEMA.germlineVariant
        WHERE
            sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' AND filter = 'PASS'
    UNION
        SELECT
            id, 'OnlyInNew' AS pipeline, chromosome, position, ref, alt, reported
        FROM
            VARIABLE_NEW_DB_SCHEMA.germlineVariant
        WHERE
            sampleId = 'VARIABLE_NEW_SAMPLE_ID' AND filter = 'PASS')
    AS a
GROUP BY 2 , 3 , 4 , 5, 6
HAVING COUNT(*) != 2;