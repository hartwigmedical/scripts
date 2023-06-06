SELECT
(SELECT COUNT(*) FROM germlineVariant
WHERE sampleId='XXX' AND filter="PASS" AND type="SNP" AND gene != '' AND rnaTotalReadCount>=10 AND rnaAlleleReadCount>=1)
/
(SELECT COUNT(*) FROM germlineVariant
WHERE sampleId='XXX' AND filter="PASS" AND type="SNP" AND gene != '' AND rnaTotalReadCount>=10)
AS fraction;