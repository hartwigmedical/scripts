SELECT structuralVariant.id, sampleId, gene, strand, isStartEnd, exonRankUpstream, exonRankDownstream, exonMax, startOrientation, endOrientation, type, ploidy FROM structuralVariantDisruption
INNER JOIN structuralVariantBreakend ON structuralVariantBreakend.id = structuralVariantDisruption.breakendId
INNER JOIN structuralVariant ON structuralVariant.id = structuralVariantBreakend.structuralVariantId
WHERE isReported = 1
AND sampleId IN ('XXX');