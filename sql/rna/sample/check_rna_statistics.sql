SELECT *, duplicates/totalFragments AS duplicateRate, totalFragments-duplicates AS nonDuplicates
FROM rnaStatistics
WHERE sampleId IN ('XXX');