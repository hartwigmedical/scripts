SELECT
    MIN(pipeline),
    startChromosome,
    startPosition,
    endChromosome,
    endPosition,
    TYPE,
    endPosition - startPosition AS length,
    COUNT(*),
    MIN(junctionCopyNumber),
    MIN(qualScore)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
            startChromosome,
            startPosition,
            endChromosome,
            endPosition,
            type,
            junctionCopyNumber,
            qualScore
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.structuralVariant
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
            AND filter <> 'PON' UNION SELECT
        'OnlyInNew' AS pipeline,
            startChromosome,
            startPosition,
            endChromosome,
            endPosition,
            type,
            junctionCopyNumber,
            qualScore
    FROM
        VARIABLE_NEW_DB_SCHEMA.structuralVariant
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID'
            AND filter <> 'PON') AS a
GROUP BY startChromosome , startPosition , endChromosome , endPosition , TYPE , length
HAVING COUNT(*) != 2;