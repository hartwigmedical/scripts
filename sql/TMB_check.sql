select sampleId, count(*)/3095 as TMB, if (count(*)/3095 > 10, "High", "Low") as status from somaticVariant where  filter in ('PASS')
group by sampleId
order by 2 desc;
