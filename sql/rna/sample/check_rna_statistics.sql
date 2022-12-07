SELECT a.sampleId, readLength, qcStatus, IF(genesWithSplicedFragments>17000,TRUE,FALSE) AS sufficientGenesWithSplicedFragments,
totalFragments, duplicates, totalFragments-duplicates AS nonDuplicates, duplicates/totalFragments AS duplicateRate,
splicedPercent, unsplicedPercent, alternateSplicePercent, chimericPercent,
round(totalFragments*splicedPercent,0) AS splicedFragments,
round(totalFragments*unsplicedPercent,0) AS unsplicedFragments,
round(totalFragments*alternateSplicePercent,0) AS alternateSpliceFragments,
round(totalFragments*chimericPercent,0) AS chimericFragments,
enrichedGenePercent,
fragmentLengthPct05, fragmentLengthPct50, fragmentLengthPct95, medianGCRatio
FROM rnaStatistics r
INNER JOIN (SELECT sampleId, count(*) AS genesWithSplicedFragments FROM geneExpression WHERE splicedFragments>0 AND sampleId IN ('XXX') GROUP BY 1) AS a
ON r.sampleId=a.sampleId;