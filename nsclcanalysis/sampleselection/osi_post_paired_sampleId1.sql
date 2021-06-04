SELECT t0.sampleId1, t0.driver, t0.qcStatus, t0.dateBiopsy1, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine1PdDate) as daysPdAfterBiopsy, 
pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine2PdDate) as daysPdAfterBiopsy, 
pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine3PdDate) as daysPdAfterBiopsy, 
pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine4PdDate) as daysPdAfterBiopsy, 
pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine5PdDate) as daysPdAfterBiopsy, 
pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine6PdDate) as daysPdAfterBiopsy, 
pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine7PdDate) as daysPdAfterBiopsy, 
pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine8PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy1Line1, posttreatmentbiopsy1Line1startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line1startDate) as daysStartAfterBiopsy, posttreatmentbiopsy1Line1PdDate, posttreatmentbiopsy1Line1StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line1PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy1Line2, posttreatmentbiopsy1Line2startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line2startDate) as daysStartAfterBiopsy, posttreatmentbiopsy1Line2PdDate, posttreatmentbiopsy1Line2StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line2PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy1Line3, posttreatmentbiopsy1Line3startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line3startDate) as daysStartAfterBiopsy, posttreatmentbiopsy1Line3PdDate, posttreatmentbiopsy1Line3StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line3PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy1Line4, posttreatmentbiopsy1Line4startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line4startDate) as daysStartAfterBiopsy, posttreatmentbiopsy1Line4PdDate, posttreatmentbiopsy1Line4StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line4PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy1Line5, posttreatmentbiopsy1Line5startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line5startDate) as daysStartAfterBiopsy, posttreatmentbiopsy1Line5PdDate, posttreatmentbiopsy1Line5StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy1Line5PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line1, posttreatmentbiopsy2Line1startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line1startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line1PdDate, posttreatmentbiopsy2Line1StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line1PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line2, posttreatmentbiopsy2Line2startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line2startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line2PdDate, posttreatmentbiopsy2Line2StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line2PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line3, posttreatmentbiopsy2Line3startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line3startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line3PdDate, posttreatmentbiopsy2Line3StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line3PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line4, posttreatmentbiopsy2Line4startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line4startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line4PdDate, posttreatmentbiopsy2Line4StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line4PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line5, posttreatmentbiopsy2Line5startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line5startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line5PdDate, posttreatmentbiopsy2Line5StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line5PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line6, posttreatmentbiopsy2Line6startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line6startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line6PdDate, posttreatmentbiopsy2Line6StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line6PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line7, posttreatmentbiopsy2Line7startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line7startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line7PdDate, posttreatmentbiopsy2Line7StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line7PdDate) as daysPdAfterBiopsy,
posttreatmentbiopsy2Line8, posttreatmentbiopsy2Line8startDate, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line8startDate) as daysStartAfterBiopsy, posttreatmentbiopsy2Line8PdDate, posttreatmentbiopsy2Line8StopReason, dateDiff(t0.dateBiopsy1,posttreatmentbiopsy2Line8PdDate) as daysPdAfterBiopsy

FROM
(select n.sampleId1, driver, qcStatus, n.dateBiopsy1 from nsclcPairedSample n inner join hmfpatients.purity p on n.sampleId1=p.sampleId inner join nsclcPairedDriver pd on n.sampleId1=pd.sampleId where (qcStatus NOT LIKE '%WARN_LOW_PURITY%' AND qcStatus <> 'FAIL_NO_TUMOR')) t0
left join 
(select sampleId1, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason from nsclcPairedSample where preTreatmentLine1 like '%Osi%') t1
on t0.sampleId1=t1.sampleId1
left join (
select sampleId1, pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason from nsclcPairedSample where preTreatmentLine2 like '%Osi%') t2
on t0.sampleId1=t2.sampleId1
left join (
select sampleId1, pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason from nsclcPairedSample where preTreatmentLine3 like '%Osi%') t3
on t0.sampleId1=t3.sampleId1
left join (
select sampleId1, pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason from nsclcPairedSample where preTreatmentLine4 like '%Osi%') t4
on t0.sampleId1=t4.sampleId1
left join (
select sampleId1, pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason from nsclcPairedSample where preTreatmentLine5 like '%Osi%') t5
on t0.sampleId1=t5.sampleId1
left join (
select sampleId1, pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason from nsclcPairedSample where preTreatmentLine6 like '%Osi%') t6
on t0.sampleId1=t6.sampleId1
left join (
select sampleId1, pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason from nsclcPairedSample where preTreatmentLine7 like '%Osi%') t7
on t0.sampleId1=t7.sampleId1
left join (
select sampleId1, pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason from nsclcPairedSample where preTreatmentLine8 like '%Osi%') t8
on t0.sampleId1=t8.sampleId1
left join 
(select sampleId1, posttreatmentbiopsy1Line1, posttreatmentbiopsy1Line1startDate, posttreatmentbiopsy1Line1PdDate, posttreatmentbiopsy1Line1StopReason from nsclcPairedSample where posttreatmentbiopsy1Line1 like '%Osi%') p1
on t0.sampleId1=p1.sampleId1
left join (
select sampleId1, posttreatmentbiopsy1Line2, posttreatmentbiopsy1Line2startDate, posttreatmentbiopsy1Line2PdDate, posttreatmentbiopsy1Line2StopReason from nsclcPairedSample where posttreatmentbiopsy1Line2 like '%Osi%') p2
on t0.sampleId1=p2.sampleId1
left join (
select sampleId1, posttreatmentbiopsy1Line3, posttreatmentbiopsy1Line3startDate, posttreatmentbiopsy1Line3PdDate, posttreatmentbiopsy1Line3StopReason from nsclcPairedSample where posttreatmentbiopsy1Line3 like '%Osi%') p3
on t0.sampleId1=p3.sampleId1
left join (
select sampleId1, posttreatmentbiopsy1Line4, posttreatmentbiopsy1Line4startDate, posttreatmentbiopsy1Line4PdDate, posttreatmentbiopsy1Line4StopReason from nsclcPairedSample where posttreatmentbiopsy1Line4 like '%Osi%') p4
on t0.sampleId1=p4.sampleId1
left join (
select sampleId1, posttreatmentbiopsy1Line5, posttreatmentbiopsy1Line5startDate, posttreatmentbiopsy1Line5PdDate, posttreatmentbiopsy1Line5StopReason from nsclcPairedSample where posttreatmentbiopsy1Line5 like '%Osi%') p5
on t0.sampleId1=p5.sampleId1
left join 
(select sampleId1, posttreatmentbiopsy2Line1, posttreatmentbiopsy2Line1startDate, posttreatmentbiopsy2Line1PdDate, posttreatmentbiopsy2Line1StopReason from nsclcPairedSample where posttreatmentbiopsy2Line1 like '%Osi%') pp1
on t0.sampleId1=pp1.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line2, posttreatmentbiopsy2Line2startDate, posttreatmentbiopsy2Line2PdDate, posttreatmentbiopsy2Line2StopReason from nsclcPairedSample where posttreatmentbiopsy2Line2 like '%Osi%') pp2
on t0.sampleId1=pp2.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line3, posttreatmentbiopsy2Line3startDate, posttreatmentbiopsy2Line3PdDate, posttreatmentbiopsy2Line3StopReason from nsclcPairedSample where posttreatmentbiopsy2Line3 like '%Osi%') pp3
on t0.sampleId1=pp3.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line4, posttreatmentbiopsy2Line4startDate, posttreatmentbiopsy2Line4PdDate, posttreatmentbiopsy2Line4StopReason from nsclcPairedSample where posttreatmentbiopsy2Line4 like '%Osi%') pp4
on t0.sampleId1=pp4.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line5, posttreatmentbiopsy2Line5startDate, posttreatmentbiopsy2Line5PdDate, posttreatmentbiopsy2Line5StopReason from nsclcPairedSample where posttreatmentbiopsy2Line5 like '%Osi%') pp5
on t0.sampleId1=pp5.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line6, posttreatmentbiopsy2Line6startDate, posttreatmentbiopsy2Line6PdDate, posttreatmentbiopsy2Line6StopReason from nsclcPairedSample where posttreatmentbiopsy2Line6 like '%Osi%') pp6
on t0.sampleId1=pp6.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line7, posttreatmentbiopsy2Line7startDate, posttreatmentbiopsy2Line7PdDate, posttreatmentbiopsy2Line7StopReason from nsclcPairedSample where posttreatmentbiopsy2Line7 like '%Osi%') pp7
on t0.sampleId1=pp7.sampleId1
left join (
select sampleId1, posttreatmentbiopsy2Line8, posttreatmentbiopsy2Line8startDate, posttreatmentbiopsy2Line8PdDate, posttreatmentbiopsy2Line8StopReason from nsclcPairedSample where posttreatmentbiopsy2Line8 like '%Osi%') pp8
on t0.sampleId1=pp8.sampleId1
WHERE 
(posttreatmentbiopsy1Line1 like '%Osi%' or posttreatmentbiopsy1Line2 like '%Osi%' or posttreatmentbiopsy1Line3 like '%Osi%' or posttreatmentbiopsy1Line4 like '%Osi%' or posttreatmentbiopsy1Line5 like '%Osi%'
or posttreatmentbiopsy2Line1 like '%Osi%' or posttreatmentbiopsy2Line2 like '%Osi%' or posttreatmentbiopsy2Line3 like '%Osi%' or posttreatmentbiopsy2Line4 like '%Osi%' or posttreatmentbiopsy2Line5 like '%Osi%' or posttreatmentbiopsy2Line6 like '%Osi%' or posttreatmentbiopsy2Line7 like '%Osi%' or posttreatmentbiopsy2Line8 like '%Osi%')
AND t0.sampleId1 not like 'CORE%';




