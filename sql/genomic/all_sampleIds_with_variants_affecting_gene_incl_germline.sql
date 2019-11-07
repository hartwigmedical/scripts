SELECT count(distinct sampleId) FROM datarequestFiltered
WHERE sampleId IN
(
		SELECT sampleId FROM somaticVariant
		WHERE gene = 'XXX' AND worstCodingEffect IN ('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE')
    UNION
        SELECT geneCopyNumber.sampleId FROM geneCopyNumber INNER JOIN purity ON geneCopyNumber.sampleId = purity.sampleId
        WHERE gene = 'XXX' AND minCopyNumber > (ploidy*3)
    UNION
        SELECT sampleId FROM geneCopyNumber WHERE gene = 'XXX' AND minCopynumber < 0.5
    UNION
        SELECT sampleId FROM svBreakend WHERE gene = 'XXX' AND disruptive = 1
    UNION
        SELECT sampleId from germlineVariant
        WHERE gene = 'XXX' AND codingEffect IN ('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE')
)