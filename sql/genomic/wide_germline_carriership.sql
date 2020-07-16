SELECT sampleId, chromosome, position, filter, type, ref, alt, gene, transcript, codingEffect, hgvsCoding, hgvsProtein, refStatus
FROM germlineVariant
WHERE sampleId LIKE '%WIDE%' and reported = 1
ORDER by sampleId;