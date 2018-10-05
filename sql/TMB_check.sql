select sampleId, count(*)/2859 as TMB, if (count(*)/2859 > 10, "High", "Low") as status from somaticVariant where  filter in ('PASS')
group by sampleId
order by 2 desc;
