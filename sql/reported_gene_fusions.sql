SELECT sampleId, startChromosome, startPosition, startOrientation, endChromosome, endPosition, endOrientation, ploidy,
fiveBreakend.gene as fiveGene, fiveBreakend.transcriptId as fiveTrranscript, fiveBreakend.strand as fiveStrand, "1" as firstFiveExon, fiveBreakend.exonRankUpstream as finalFiveExon, 
threeBreakend.gene as threeGene, threeBreakend.transcriptId as threeTranscript, threeBreakend.strand as threeStrand, threeBreakend.exonRankDownstream as firstThreeExon, threeBreakend.exonMax as finalThreeExon 
FROM structuralVariantFusion
INNER JOIN structuralVariantBreakend AS fiveBreakend ON fiveBreakend.id = structuralVariantFusion.fivePrimeBreakendId 
INNER JOIN structuralVariantBreakend AS threeBreakend ON threeBreakend.id = structuralVariantFusion.threePrimeBreakendId 
INNER JOIN structuralVariant  ON structuralVariant.id = fiveBreakend.structuralVariantId
WHERE isReported = 1
AND sampleId IN  ('XXX');