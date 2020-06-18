SELECT 
    MIN(pipeline),
    startChromosome,
    startPosition,
    endChromosome,
    endPosition,
    TYPE,
    endPosition - startPosition AS len,
    COUNT(*),
    MAX(ploidy),
    MAX(qualScore)
FROM
    (SELECT 
        id,
            'Truth' AS pipeline,
            startChromosome,
            startPosition,
            endChromosome,
            endPosition,
            TYPE,
            ploidy,
            qualScore
    FROM
        reference_validation_sets.structuralVariant
    WHERE
        sampleId = 'COLO829v003T' AND filter <> 'PON' UNION SELECT 
        id,
            'New' AS pipeline,
            startChromosome,
            startPosition,
            endChromosome,
            endPosition,
            TYPE,
            ploidy,
            qualScore
    FROM
        pipeline_v5_validation.structuralVariant
    WHERE
        sampleId = 'COLO829v003T' AND filter <> 'PON') AS a
GROUP BY 2 , 3 , 4 , 5 , 6 , 7
HAVING COUNT(*) = 1;
