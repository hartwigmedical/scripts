SELECT
    MIN(pipeline) AS pipeline,
    gene,
    chromosome,
    positionV37,
    positionV38,
    refV37,
    refV38,
    allele1,
    allele2,
    rsid,
    variantAnnotationV37,
    filterV37,
    variantAnnotationV38,
    filterV38,
    COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
            gene,
            chromosome,
            positionV37,
            positionV38,
            refV37,
            refV38,
            allele1,
            allele2,
            rsid,
            variantAnnotationV37,
            filterV37,
            variantAnnotationV38,
            filterV38
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.peachCalls
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline,
            gene,
            chromosome,
            positionV37,
            positionV38,
            refV37,
            refV38,
            allele1,
            allele2,
            rsid,
            variantAnnotationV37,
            filterV37,
            variantAnnotationV38,
            filterV38
    FROM
        VARIABLE_NEW_DB_SCHEMA.peachCalls
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY gene , chromosome , positionV37 , positionV38 , refV37 , refV38 , allele1 , allele2 , rsid , variantAnnotationV37 , filterV37 , variantAnnotationV38 , filterV38
HAVING COUNT(*) != 2;