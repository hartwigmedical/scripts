SELECT fiveVariant.sampleId, fiveVariant.id, threeVariant.id, fiveBreakend.gene, fiveBreakend.exonRankUpstream, fiveBreakend.exonRankDownstream, 
threeBreakend.gene, threeBreakend.exonRankUpstream, threeBreakend.exonRankDownstream, fiveVariant.ploidy FROM structuralVariantFusion
INNER JOIN structuralVariantBreakend AS fiveBreakend ON fiveBreakend.id = structuralVariantFusion.fivePrimeBreakendId 
INNER JOIN structuralVariant AS fiveVariant ON fiveVariant.id = fiveBreakend.structuralVariantId
INNER JOIN structuralVariantBreakend AS threeBreakend ON threeBreakend.id = structuralVariantFusion.threePrimeBreakendId 
INNER JOIN structuralVariant AS threeVariant ON threeVariant.id = threeBreakend.structuralVariantId
WHERE isReported = 1
AND fiveVariant.sampleId IN ('XXX');