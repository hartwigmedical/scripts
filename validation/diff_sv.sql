SELECT 
    MIN(pipeline), startChromosome, startPosition, endChromosome, endPosition, TYPE, endPosition - startPosition AS len,
    COUNT(*), MAX(junctionCopyNumber), MAX(qualScore)
FROM
        (SELECT
            id, 'OnlyInTruth' AS pipeline, startChromosome, startPosition, endChromosome, endPosition, type, junctionCopyNumber, qualScore
        FROM
            VARIABLE_TRUTH_DB_SCHEMA.structuralVariant
        WHERE
            sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' AND filter <> 'PON'
    UNION
        SELECT
            id, 'OnlyInNew' AS pipeline, startChromosome, startPosition, endChromosome, endPosition, type, junctionCopyNumber, qualScore
        FROM
            VARIABLE_NEW_DB_SCHEMA.structuralVariant
        WHERE
            sampleId = 'VARIABLE_NEW_SAMPLE_ID' AND filter <> 'PON')
    AS a
GROUP BY 2 , 3 , 4 , 5 , 6 , 7
HAVING COUNT(*) = 1;
