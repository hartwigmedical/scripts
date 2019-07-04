SELECT sampleId, count(*)/2859 AS indelsPerMb, if(count(*)/2859 > 4, "MSI", "MSS" ) AS status FROM somaticVariant
WHERE filter = 'PASS' AND TYPE = 'INDEL' AND repeatCount >= 4 AND length(alt) <= 50 AND length(ref) <= 50
AND (
	(length(repeatSequence) BETWEEN 2 AND 4 ) OR
	(length(repeatSequence) = 1 AND repeatCount >= 5)
)
AND sampleId IN ('XXX')
GROUP BY sampleId
ORDER BY 2 DESC;