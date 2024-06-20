use hmfpatients_pilot;

## Lung cancer Sig 4 >15% allocation v2 -> 410
select count(distinct s.sampleId) as samplecount from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'Sig4' and c.primaryTumorLocation = 'Lung' and percent >=0.15;

## Lung cancer Sig 4 >15% allocation v3 -> 414
select count(distinct s.sampleId) from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS4_P' and c.primaryTumorLocation = 'Lung' and percent >=0.15;

## Lung cancer Sig 4 v2 among v3 patients -> 404 (98.5% of samples in v2 overlap with v3)
select count(distinct v2.sampleId) as common_samplecount
from (
    select distinct s.sampleId 
    from hmfpatients.signature s
    inner join hmfpatients.hpc h on h.sampleId = s.sampleId
    inner join hmfpatients.clinical c on c.sampleId = s.sampleId
    where s.signature = 'Sig4' 
      and c.primaryTumorLocation = 'Lung' 
      and percent >= 0.15
) v2
inner join (
    select distinct s.sampleId 
    from hmfpatients_pilot.signature s
    inner join hmfpatients.hpc h on h.sampleId = s.sampleId
    inner join hmfpatients.clinical c on c.sampleId = s.sampleId
    where s.signature = 'BI_COMPOSITE_SNV_SBS4_P' 
      and c.primaryTumorLocation = 'Lung' 
      and percent >= 0.15
) v3 on v2.sampleId = v3.sampleId;


## Misallocation per sample v2
select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'MISALLOC';

## Signature 7 >15% in v2 - skin cancer -> 293
select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'Sig7' and c.primaryTumorLocation = 'Skin' and percent >= 0.15;

## Signature 7a/b/c >15% in v3 - skin cancer
select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS7a_S', 'BI_COMPOSITE_SNV_SBS7b_S','BI_COMPOSITE_SNV_SBS7c_S') and c.primaryTumorLocation = 'Skin';

## SUM OF Signature 7a/b/c >15% in v3 - skin cancer - count -> 290
select s.sampleId, SUM(s.percent) as sig7_total_percent from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS7a_S', 'BI_COMPOSITE_SNV_SBS7b_S','BI_COMPOSITE_SNV_SBS7c_S') and c.primaryTumorLocation = 'Skin'
group by s.sampleId
having SUM(s.percent) >0.15;

## Overlap of v2 Skin cancer Sig7 >15% with v3 7a/b/c >15% -> 290 (98.98% of v2 samples overlap with v3)
select count(distinct v2.sampleId) as common_samplecount
from (
    select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.clinical c on c.sampleId = s.sampleId
	where s.signature = 'Sig7' and c.primaryTumorLocation = 'Skin' and percent >= 0.15
) v2
inner join (
    select s.sampleId, SUM(s.percent) as sig7_total_percent from hmfpatients_pilot.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.clinical c on c.sampleId = s.sampleId
	where s.signature in ('BI_COMPOSITE_SNV_SBS7a_S', 'BI_COMPOSITE_SNV_SBS7b_S','BI_COMPOSITE_SNV_SBS7c_S') and c.primaryTumorLocation = 'Skin'
	group by s.sampleId
	having SUM(s.percent) >0.15
) v3 on v2.sampleId = v3.sampleId;

## SUM OF signature 7a/b/c, SBS38, SBS55, SBS65, SBS67, SBS75 >15% in v3 - skin cancer - count -> 300
select s.sampleId, SUM(s.percent) as UV_sig_total_percent from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS7a_S', 'BI_COMPOSITE_SNV_SBS7b_S','BI_COMPOSITE_SNV_SBS7c_S', 'BI_COMPOSITE_SNV_SBS38_S', 'BI_COMPOSITE_SNV_SBS55_S', 'BI_COMPOSITE_SNV_SBS65_S', 'BI_COMPOSITE_SNV_SBS67_S', 'BI_COMPOSITE_SNV_SBS75_S') 
and c.primaryTumorLocation = 'Skin'
group by s.sampleId
having SUM(s.percent) >0.15;

## Overlap of v2 Skin cancer Sig7 >15% with v3 7a/b/c, SBS38, SBS55, SBS65, SBS67, SBS75 >15% -> 293 (100% of v2 samples overlap with v3)
select count(distinct v2.sampleId) as common_samplecount
from (
    select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.clinical c on c.sampleId = s.sampleId
	where s.signature = 'Sig7' and c.primaryTumorLocation = 'Skin' and percent >= 0.15
) v2
inner join (
    select s.sampleId, SUM(s.percent) as sig7_total_percent from hmfpatients_pilot.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.clinical c on c.sampleId = s.sampleId
	where s.signature in ('BI_COMPOSITE_SNV_SBS7a_S', 'BI_COMPOSITE_SNV_SBS7b_S','BI_COMPOSITE_SNV_SBS7c_S', 'BI_COMPOSITE_SNV_SBS38_S', 'BI_COMPOSITE_SNV_SBS55_S', 'BI_COMPOSITE_SNV_SBS65_S', 'BI_COMPOSITE_SNV_SBS67_S', 'BI_COMPOSITE_SNV_SBS75_S')
    and c.primaryTumorLocation = 'Skin'
	group by s.sampleId
	having SUM(s.percent) >0.15
) v3 on v2.sampleId = v3.sampleId;

## Total MSI samples -> 202
select distinct s.sampleId from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where msStatus = 'MSI';

## Total MS-stable samples -> 6353
select distinct s.sampleId from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where msStatus = 'MSS';

## Sum of signature 6/15/21/26 >15% in MSI samples v2 -> 186
select s.sampleId, SUM(s.percent) from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where s.signature in ('Sig6', 'Sig15', 'Sig21', 'Sig26') and msStatus = 'MSI'
group by s.sampleId
having SUM(s.percent) >0.15;

## Sum of signature 6/15/21/26 >15% in MSI samples v3 -> 176
select s.sampleId, SUM(s.percent) from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS6_S', 'BI_COMPOSITE_SNV_SBS15_S', 'BI_COMPOSITE_SNV_SBS21_S', 'BI_COMPOSITE_SNV_SBS26_S') and msStatus = 'MSI'
group by s.sampleId
having SUM(s.percent) >0.15;

## Overlap of 6/15/21/26 v2 among v3 samples -> 173 v2 samples are in the v3 set
select count(distinct v2.sampleId) as common_samplecount
from (
    select s.sampleId, SUM(s.percent) from hmfpatients.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.purity p on p.sampleId = s.sampleId
	where s.signature in ('Sig6', 'Sig15', 'Sig21', 'Sig26') and msStatus = 'MSI'
	group by s.sampleId
	having SUM(s.percent) >0.15
) v2
inner join (
    select s.sampleId, SUM(s.percent) from hmfpatients_pilot.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.purity p on p.sampleId = s.sampleId
	where s.signature in ('BI_COMPOSITE_SNV_SBS6_S', 'BI_COMPOSITE_SNV_SBS15_S', 'BI_COMPOSITE_SNV_SBS21_S', 'BI_COMPOSITE_SNV_SBS26_S') and msStatus = 'MSI'
	group by s.sampleId
	having SUM(s.percent) >0.15
) v3 on v2.sampleId = v3.sampleId;

## Sum of signature 6/15/21/26 >15% in MS-stable samples v2 -> 57
select s.sampleId, SUM(s.percent) from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where s.signature in ('Sig6', 'Sig15', 'Sig21', 'Sig26') and msStatus = 'MSS'
group by s.sampleId
having SUM(s.percent) >0.15;

## Sum of signature 6/15/21/26 >15% in MS-stable samples v3 -> 13
select s.sampleId, SUM(s.percent) from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS6_S', 'BI_COMPOSITE_SNV_SBS15_S', 'BI_COMPOSITE_SNV_SBS21_S', 'BI_COMPOSITE_SNV_SBS26_S') and msStatus = 'MSS'
group by s.sampleId
having SUM(s.percent) >0.15;

## Sum of signature 6/15/21/26/44/73/76/79 >15% in MS-stable samples v3 -> 38
select s.sampleId, SUM(s.percent) from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.purity p on p.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS6_S', 'BI_COMPOSITE_SNV_SBS15_S', 'BI_COMPOSITE_SNV_SBS21_S', 'BI_COMPOSITE_SNV_SBS26_S', 'BI_COMPOSITE_SNV_SBS44_S', 'BI_COMPOSITE_SNV_SBS73_S', 'BI_COMPOSITE_SNV_SBS76_S', 'BI_COMPOSITE_SNV_SBS79_S') and msStatus = 'MSS'
group by s.sampleId
having SUM(s.percent) >0.15;

## Total HRD samples -> 392
select distinct s.sampleId from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where hrStatus = 'HR_DEFICIENT';

## Total HRP samples -> 5905
select distinct s.sampleId from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where hrStatus = 'HR_PROFICIENT';

## Signature 3 (HRD) >15% in HR-proficient samples - v2 -> 1313
select s.sampleId, signature, allocation, percent, hrStatus from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where s.signature = 'Sig3' and s.percent >=0.15 and hrStatus = 'HR_PROFICIENT';

## Signature 3 (HRD) in HR-proficient samples - v3 -> 16
select s.sampleId, signature, allocation, percent, hrStatus from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS3_P' and s.percent >=0.15 and hrStatus = 'HR_PROFICIENT';

## Signature 3 (HRD) in HR-deficient samples - v2 -> 369
select s.sampleId, signature, allocation, percent, hrStatus from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where s.signature = 'Sig3' and s.percent >=0.15 and hrStatus = 'HR_DEFICIENT';

## Signature 3 (HRD) in HR-deficient samples - v3 -> 3
select s.sampleId, signature, allocation, percent, hrStatus from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.chord c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS3_P' and s.percent >=0.15 and hrStatus = 'HR_DEFICIENT';

## Kidney cancer total in db -> 149
select distinct s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
and c.primaryTumorLocation = 'Kidney';

## Signature 40 (new in v3) - kidney cancer -> 50
select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS40_P'
and percent >= 0.15
and c.primaryTumorLocation = 'Kidney';

## Signature 40 (new in v3) - all except kidney cancer -> 2
select s.sampleId, signature, allocation, percent, primaryTumorLocation, primaryTumorType, primaryTumorSubType from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS40_P'
and percent >= 0.15
and c.primaryTumorLocation != 'Kidney' and c.primaryTumorLocation not like 'Unknown%';

## Prostate cancer total in db -> 475
select distinct s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
and c.primaryTumorLocation = 'Prostate';

## Signature 80 (new in v3) >15% - prostate cancer -> 
select distinct s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS80_P'
and percent >= 0.15
and c.primaryTumorLocation = 'Prostate';

## Signature 80 (new in v3) >15% - all except prostate cancer -> 
select distinct s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'BI_COMPOSITE_SNV_SBS80_P'
and percent >= 0.15
and c.primaryTumorLocation != 'Prostate' and c.primaryTumorLocation not like 'Unknown%';

## Total in db with HER2amp and/or PTEN loss
select distinct s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
where gene in ('ERBB2', 'PTEN') 
and driver in ('AMP', 'DEL');

## Sum of APOBEC signatures v2 (2, 13) >15% in HER2 and PTEN loss tumors -> 177
select s.sampleId, SUM(s.percent) from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
where s.signature in ('Sig2', 'Sig13')
and gene in ('ERBB2', 'PTEN') 
and driver in ('AMP', 'DEL')
group by sampleId
having SUM(s.percent) >0.15;

## Sum of APOBEC signatures v3 (2, 13) >15% in HER2 and PTEN loss tumors -> 150
select s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS2_P', 'BI_COMPOSITE_SNV_SBS13_P')
and gene in ('ERBB2', 'PTEN') 
and driver in ('AMP', 'DEL')
group by sampleId
having SUM(s.percent) >0.15;

## Sum of APOBEC signatures v3 (2, 13, 69) >15% in HER2 and PTEN loss tumors -> 192
select s.sampleId from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
where s.signature in ('BI_COMPOSITE_SNV_SBS2_P', 'BI_COMPOSITE_SNV_SBS13_P', 'BI_COMPOSITE_SNV_SBS69_P')
and gene in ('ERBB2', 'PTEN') 
and driver in ('AMP', 'DEL')
group by sampleId
having SUM(s.percent) >0.15;

## Overlap of v2 (2,13) and v3 (2,13,69) -> 171 of 177 samples in v2 are also among samples in v3
select count(distinct v2.sampleId) as common_samplecount
from (
    select s.sampleId, SUM(s.percent) from hmfpatients.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
	where s.signature in ('Sig2', 'Sig13')
	and gene in ('ERBB2', 'PTEN') 
	and driver in ('AMP', 'DEL')
	group by sampleId
	having SUM(s.percent) >0.15
) v2
inner join (
    select s.sampleId from hmfpatients_pilot.signature s
	inner join hmfpatients.hpc h on h.sampleId = s.sampleId
	inner join hmfpatients.driverCatalog d on d.sampleId = s.sampleId
	where s.signature in ('BI_COMPOSITE_SNV_SBS2_P', 'BI_COMPOSITE_SNV_SBS13_P', 'BI_COMPOSITE_SNV_SBS69_P')
	and gene in ('ERBB2', 'PTEN') 
	and driver in ('AMP', 'DEL')
	group by sampleId
	having SUM(s.percent) >0.15
) v3 on v2.sampleId = v3.sampleId;

## Percentage of SBS1 in v2 and v3
select percent from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
where s.signature = 'Sig1';

select signature, percent from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
where s.signature like '%BI_COMPOSITE_SNV_SBS1_P%';

select distinct sampleId from hmfpatients_pilot.signature s; ##8075
select distinct sampleId from hmfpatients.signature s; ##8079