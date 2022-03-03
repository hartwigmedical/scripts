SELECT
    MIN(pipeline) AS pipeline, gene, haplotype, function, linkedDrugs, COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline, gene, haplotype, function, linkedDrugs
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.peachGenotype
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline, gene, haplotype, function, linkedDrugs
    FROM
        VARIABLE_NEW_DB_SCHEMA.peachGenotype
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY gene, haplotype, function, linkedDrugs
HAVING COUNT(*) != 2;