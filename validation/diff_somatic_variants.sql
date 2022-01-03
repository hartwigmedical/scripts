SELECT 
    MIN(pipeline) AS pipeline, chromosome, position, ref, alt, COUNT(*), MIN(adjustedVaf)
FROM
        (SELECT
            id, 'OnlyInTruth' AS pipeline, chromosome, position, ref, alt, adjustedVaf
        FROM
            VARIABLE_TRUTH_DB_SCHEMA.somaticVariant
        WHERE
            sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' AND filter = 'PASS'
    UNION
        SELECT
            id, 'OnlyInNew' AS pipeline, chromosome, position, ref, alt, adjustedVaf
        FROM
            VARIABLE_NEW_DB_SCHEMA.somaticVariant
        WHERE
            sampleId = 'VARIABLE_NEW_SAMPLE_ID' AND filter = 'PASS')
    AS a
GROUP BY 2 , 3 , 4 , 5
HAVING COUNT(*) != 2;
