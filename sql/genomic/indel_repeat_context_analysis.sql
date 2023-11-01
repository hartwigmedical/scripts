SELECT sampleId,
least(5,greatest(-5,length(alt)-length(ref))) AS indelLength,
count(*) AS totalCount,
sum(IF(repeatCount>=4,1,0)) AS repeatCount,
sum(IF(repeatCount<4 AND microhomology <>'' AND microhomology <> '.',1,0)) AS MHCount,
sum(IF(repeatCount<4 AND (microhomology ='' OR microhomology='.'),1,0)) AS nonRepeatNonMHCount
FROM somaticVariant WHERE sampleId = 'XXX' AND type = 'INDEL' AND filter = 'PASS' GROUP BY 1,2;

