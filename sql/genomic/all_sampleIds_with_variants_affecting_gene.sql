SELECT count(distinct datarequestFiltered.sampleId) FROM hmfpatients.datarequestFiltered
INNER JOIN
((SELECT sampleId FROM hmfpatients.somaticVariant
WHERE gene = 'XXX' AND (worstCodingEffect = 'MISSENSE' OR worstCodingEffect = 'NONSENSE_OR_FRAMESHIFT' OR worstCodingEffect = 'SPLICE'))
UNION
(SELECT geneCopyNumber.sampleId FROM hmfpatients.geneCopyNumber
INNER JOIN
(SELECT sampleId, ploidy FROM hmfpatients.purity) as b on geneCopyNumber.sampleId = b.sampleId
WHERE gene = 'XXX' AND minCopyNumber > (ploidy*3))
UNION
(SELECT sampleId FROM geneCopyNumber WHERE gene = 'XXX' AND minCopynumber < 0.5)
UNION
(SELECT sampleId FROM svBreakend WHERE gene = 'XXX' AND disruptive = 1)
UNION
(SELECT sampleId from svFusion WHERE name like '%XXX%' AND phased = 1) )
as a on datarequestFiltered.sampleId = a.sampleId;