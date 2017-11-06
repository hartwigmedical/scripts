select s.sampleId as Sample, count(*)/3095 as SIndel_PerMb, if(count(*)/3095 > 0.909, "MSI", "MSS" ) as Status from somaticVariant s left join purity on s.sampleId = purity.sampleId
where filter = 'PASS' and repeatCount > 0 and length(alt) <> length(ref) and length(alt) <= 50 and length(ref) <= 50
  and (
(length(repeatSequence) between 2 and 4 ) OR
(length(repeatSequence) = 1 and repeatCount >= 5)
)
and s.sampleId IN ('XXX')
group by s.sampleId
order by 2 desc;