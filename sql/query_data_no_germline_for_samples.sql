select sampleId, chromosome, start, end, segmentStartSupport, segmentEndSupport, bafCount, observedBaf, actualBaf, copyNumber, copyNumberMethod from copyNumber where sampleId in ('XXX');

select sampleId, chromosome, start, end, svCluster, ratioSupport, segmentStartSupport, bafCount, observedBaf, observedTumorRatio, 
observedNormalRatio, observedTumorRatioCount, gcContent, modelPloidy, modelBaf, modelTumorRatio, actualTumorBaf, actualTumorCopyNumber, 
refNormalisedTumorCopyNumber, cnvDeviation, bafDeviation, totalDeviation, ploidyPenalty, fittedBaf, fittedCopyNumber from copyNumberRegion where sampleId in ('XXX');

select sampleId, chromosome, start, end, gene, chromosomeBand, transcriptId, transcriptVersion, minCopyNumber, maxCopyNumber, 
somaticRegions,  minRegions, minRegionStart, minRegionEnd, minRegionStartSupport, minRegionEndSupport,
minRegionMethod, nonsenseBiallelicVariants, nonsenseNonBiallelicVariants, nonsenseNonBiallelicPloidy, spliceBiallelicVariants, spliceNonBiallelicVariants,
spliceNonBiallelicPloidy, missenseBiallelicVariants, missenseNonBiallelicVariants, missenseNonBiallelicPloidy, minMinorAllelePloidy
 from geneCopyNumber where sampleId in ('XXX');

select sampleId, purity, normFactor, score, ploidy, diploidProportion from purityRange where sampleId in ('XXX');

select sampleId, chromosome, position, filter, type, ref, alt, gene, cosmicId, dbsnpId, effect, codingEffect, microhomology,
repeatSequence, repeatCount, alleleReadCount, totalReadCount, adjustedVaf, adjustedCopyNumber, highConfidence, trinucleotideContext,
clonality, biallelic, hotspot, mappability, minorAllelePloidy from somaticVariant where sampleId in ('XXX');

select sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, startHomologySequence,
endHomologySequence, startAF, endAF, ploidy, adjustedStartAF, adjustedEndAF, adjustedStartCopyNumber, adjustedEndCopyNumber, adjustedStartCopyNumberChange, adjustedEndCopyNumberChange
insertSequence, type from structuralVariant where sampleId in ('XXX');

select sampleId, gender, status, qcStatus, purity, normFactor, score, ploidy, diploidProportion, polyclonalProportion,
minPurity, maxPurity, minPloidy, maxPloidy, minDiploidProportion, maxDiploidProportion
 from purity where sampleId in ('XXX');







