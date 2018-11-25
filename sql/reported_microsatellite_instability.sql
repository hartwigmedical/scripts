select sampleId, count(*)/2859 as indelsPerMb, if(count(*)/2859 > 4, "MSI", "MSS" ) as status from somaticVariant
where filter = 'PASS' and type = 'INDEL' and repeatCount >= 4 and length(alt) <= 50 and length(ref) <= 50
and (
	(length(repeatSequence) between 2 and 4 ) OR
	(length(repeatSequence) = 1 and repeatCount >= 5)
)
and sampleId IN ('XXX')
group by sampleId
order by 2 desc;