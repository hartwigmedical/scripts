SELECT fiveVariant.sampleId, fiveVariant.startOrientation, fiveVariant.startChromosome, fiveVariant.startPosition, fiveVariant.endChromosome, fiveVariant.endOrientation, fiveVariant.endPosition,
fiveBreakend.gene as fiveGene, fiveBreakend.strand as fiveStrand, "1" as firstFiveExon, fiveBreakend.exonRankUpstream as finalFiveExon, 
threeBreakend.gene as threeGene, threeBreakend.strand as threeStrand, threeBreakend.exonRankDownstream as firstThreeExon, threeBreakend.exonMax as finalThreeExon, fiveVariant.ploidy FROM structuralVariantFusion
INNER JOIN structuralVariantBreakend AS fiveBreakend ON fiveBreakend.id = structuralVariantFusion.fivePrimeBreakendId 
INNER JOIN structuralVariant AS fiveVariant ON fiveVariant.id = fiveBreakend.structuralVariantId
INNER JOIN structuralVariantBreakend AS threeBreakend ON threeBreakend.id = structuralVariantFusion.threePrimeBreakendId 
INNER JOIN structuralVariant AS threeVariant ON threeVariant.id = threeBreakend.structuralVariantId
WHERE isReported = 1
AND fiveVariant.sampleId IN  ('XXX');