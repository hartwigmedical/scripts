SELECT sampleId, readLength, qcStatus, totalFragments, duplicates,
totalFragments-duplicates AS nonDuplicates,
duplicates/totalFragments AS duplicateRate,
splicedPercent, unsplicedPercent, alternateSplicePercent, chimericPercent,
round(totalFragments*splicedPercent,0) AS splicedFragments,
round(totalFragments*alternateSplicePercent,0) AS alternateSpliceFragments,
round(totalFragments*chimericPercent,0) AS chimericFragments,
fragmentLengthPct05, fragmentLengthPct50, fragmentLengthPct95, enrichedGenePercent, medianGCRatio
FROM rnaStatistics
WHERE sampleId IN ('XXX');