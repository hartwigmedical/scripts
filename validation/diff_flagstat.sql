SELECT
    MIN(pipeline) AS pipeline,
    sampleId, refUniqueReadCount, refSecondaryCount, refSupplementaryCount, refDuplicateProportion, refMappedProportion, refPairedInSequencingProportion, refProperlyPairedProportion, refWithItselfAndMateMappedProportion, refSingletonProportion, tumorUniqueReadCount, tumorSecondaryCount, tumorSupplementaryCount, tumorDuplicateProportion, tumorMappedProportion, tumorPairedInSequencingProportion, tumorProperlyPairedProportion, tumorWithItselfAndMateMappedProportion, tumorSingletonProportion, passQC
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
        sampleId, refUniqueReadCount, refSecondaryCount, refSupplementaryCount, refDuplicateProportion, refMappedProportion, refPairedInSequencingProportion, refProperlyPairedProportion, refWithItselfAndMateMappedProportion, refSingletonProportion, tumorUniqueReadCount, tumorSecondaryCount, tumorSupplementaryCount, tumorDuplicateProportion, tumorMappedProportion, tumorPairedInSequencingProportion, tumorProperlyPairedProportion, tumorWithItselfAndMateMappedProportion, tumorSingletonProportion, passQC
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.flagstat
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline,
        sampleId, refUniqueReadCount, refSecondaryCount, refSupplementaryCount, refDuplicateProportion, refMappedProportion, refPairedInSequencingProportion, refProperlyPairedProportion, refWithItselfAndMateMappedProportion, refSingletonProportion, tumorUniqueReadCount, tumorSecondaryCount, tumorSupplementaryCount, tumorDuplicateProportion, tumorMappedProportion, tumorPairedInSequencingProportion, tumorProperlyPairedProportion, tumorWithItselfAndMateMappedProportion, tumorSingletonProportion, passQC
    FROM
        VARIABLE_NEW_DB_SCHEMA.flagstat
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY sampleId, refUniqueReadCount, refSecondaryCount, refSupplementaryCount, refDuplicateProportion, refMappedProportion, refPairedInSequencingProportion, refProperlyPairedProportion, refWithItselfAndMateMappedProportion, refSingletonProportion, tumorUniqueReadCount, tumorSecondaryCount, tumorSupplementaryCount, tumorDuplicateProportion, tumorMappedProportion, tumorPairedInSequencingProportion, tumorProperlyPairedProportion, tumorWithItselfAndMateMappedProportion, tumorSingletonProportion, passQC
HAVING COUNT(*) != 2;