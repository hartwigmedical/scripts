SELECT t0.sampleId, t0.driver, t0.qcStatus, t0.dateBiopsy, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine1PdDate) as daysPdAfterBiopsy, 
pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine2PdDate) as daysPdAfterBiopsy, 
pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine3PdDate) as daysPdAfterBiopsy, 
pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine4PdDate) as daysPdAfterBiopsy, 
pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine5PdDate) as daysPdAfterBiopsy, 
pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine6PdDate) as daysPdAfterBiopsy, 
pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine7PdDate) as daysPdAfterBiopsy, 
pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine8PdDate) as daysPdAfterBiopsy,
posttreatmentLine1, posttreatmentLine1startDate, dateDiff(t0.dateBiopsy,posttreatmentLine1startDate) as daysStartAfterBiopsy, posttreatmentLine1PdDate, posttreatmentLine1StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine1PdDate) as daysPdAfterBiopsy,
posttreatmentLine2, posttreatmentLine2startDate, dateDiff(t0.dateBiopsy,posttreatmentLine2startDate) as daysStartAfterBiopsy, posttreatmentLine2PdDate, posttreatmentLine2StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine2PdDate) as daysPdAfterBiopsy,
posttreatmentLine3, posttreatmentLine3startDate, dateDiff(t0.dateBiopsy,posttreatmentLine3startDate) as daysStartAfterBiopsy, posttreatmentLine3PdDate, posttreatmentLine3StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine3PdDate) as daysPdAfterBiopsy,
posttreatmentLine4, posttreatmentLine4startDate, dateDiff(t0.dateBiopsy,posttreatmentLine4startDate) as daysStartAfterBiopsy, posttreatmentLine4PdDate, posttreatmentLine4StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine4PdDate) as daysPdAfterBiopsy,
posttreatmentLine5, posttreatmentLine5startDate, dateDiff(t0.dateBiopsy,posttreatmentLine5startDate) as daysStartAfterBiopsy, posttreatmentLine5PdDate, posttreatmentLine5StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine5PdDate) as daysPdAfterBiopsy,
posttreatmentLine6, posttreatmentLine6startDate, dateDiff(t0.dateBiopsy,posttreatmentLine6startDate) as daysStartAfterBiopsy, posttreatmentLine6PdDate, posttreatmentLine6StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine6PdDate) as daysPdAfterBiopsy,
posttreatmentLine7, posttreatmentLine7startDate, dateDiff(t0.dateBiopsy,posttreatmentLine7startDate) as daysStartAfterBiopsy, posttreatmentLine7PdDate, posttreatmentLine7StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine7PdDate) as daysPdAfterBiopsy,
posttreatmentLine8, posttreatmentLine8startDate, dateDiff(t0.dateBiopsy,posttreatmentLine8startDate) as daysStartAfterBiopsy, posttreatmentLine8PdDate, posttreatmentLine8StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine8PdDate) as daysPdAfterBiopsy,
posttreatmentLine9, posttreatmentLine9startDate, dateDiff(t0.dateBiopsy,posttreatmentLine9startDate) as daysStartAfterBiopsy, posttreatmentLine9PdDate, posttreatmentLine9StopReason, dateDiff(t0.dateBiopsy,posttreatmentLine9PdDate) as daysPdAfterBiopsy
FROM
(select n.sampleId, driver, qcStatus, n.dateBiopsy from nsclcSingleSample n inner join hmfpatients.purity p on n.sampleId=p.sampleId where (qcStatus NOT LIKE '%WARN_LOW_PURITY%' AND qcStatus <> 'FAIL_NO_TUMOR')) t0
left join 
(select sampleId, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason from nsclcSingleSample where preTreatmentLine1 like '%Osi%') t1
on t0.sampleId=t1.sampleId
left join (
select sampleId, pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason from nsclcSingleSample where preTreatmentLine2 like '%Osi%') t2
on t0.sampleId=t2.sampleId
left join (
select sampleId, pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason from nsclcSingleSample where preTreatmentLine3 like '%Osi%') t3
on t0.sampleId=t3.sampleId
left join (
select sampleId, pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason from nsclcSingleSample where preTreatmentLine4 like '%Osi%') t4
on t0.sampleId=t4.sampleId
left join (
select sampleId, pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason from nsclcSingleSample where preTreatmentLine5 like '%Osi%') t5
on t0.sampleId=t5.sampleId
left join (
select sampleId, pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason from nsclcSingleSample where preTreatmentLine6 like '%Osi%') t6
on t0.sampleId=t6.sampleId
left join (
select sampleId, pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason from nsclcSingleSample where preTreatmentLine7 like '%Osi%') t7
on t0.sampleId=t7.sampleId
left join (
select sampleId, pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason from nsclcSingleSample where preTreatmentLine8 like '%Osi%') t8
on t0.sampleId=t8.sampleId
left join 
(select sampleId, posttreatmentLine1, posttreatmentLine1startDate, posttreatmentLine1PdDate, posttreatmentLine1StopReason from nsclcSingleSample where posttreatmentLine1 like '%Osi%') p1
on t0.sampleId=p1.sampleId
left join (
select sampleId, posttreatmentLine2, posttreatmentLine2startDate, posttreatmentLine2PdDate, posttreatmentLine2StopReason from nsclcSingleSample where posttreatmentLine2 like '%Osi%') p2
on t0.sampleId=p2.sampleId
left join (
select sampleId, posttreatmentLine3, posttreatmentLine3startDate, posttreatmentLine3PdDate, posttreatmentLine3StopReason from nsclcSingleSample where posttreatmentLine3 like '%Osi%') p3
on t0.sampleId=p3.sampleId
left join (
select sampleId, posttreatmentLine4, posttreatmentLine4startDate, posttreatmentLine4PdDate, posttreatmentLine4StopReason from nsclcSingleSample where posttreatmentLine4 like '%Osi%') p4
on t0.sampleId=p4.sampleId
left join (
select sampleId, posttreatmentLine5, posttreatmentLine5startDate, posttreatmentLine5PdDate, posttreatmentLine5StopReason from nsclcSingleSample where posttreatmentLine5 like '%Osi%') p5
on t0.sampleId=p5.sampleId
left join (
select sampleId, posttreatmentLine6, posttreatmentLine6startDate, posttreatmentLine6PdDate, posttreatmentLine6StopReason from nsclcSingleSample where posttreatmentLine6 like '%Osi%') p6
on t0.sampleId=p6.sampleId
left join (
select sampleId, posttreatmentLine7, posttreatmentLine7startDate, posttreatmentLine7PdDate, posttreatmentLine7StopReason from nsclcSingleSample where posttreatmentLine7 like '%Osi%') p7
on t0.sampleId=p7.sampleId
left join (
select sampleId, posttreatmentLine8, posttreatmentLine8startDate, posttreatmentLine8PdDate, posttreatmentLine8StopReason from nsclcSingleSample where posttreatmentLine8 like '%Osi%') p8
on t0.sampleId=p8.sampleId
left join (
select sampleId, posttreatmentLine9, posttreatmentLine9startDate, posttreatmentLine9PdDate, posttreatmentLine9StopReason from nsclcSingleSample where posttreatmentLine9 like '%Osi%') p9
on t0.sampleId=p9.sampleId
WHERE 
(posttreatmentLine1 like '%Osi%' or posttreatmentLine2 like '%Osi%' or posttreatmentLine3 like '%Osi%' or posttreatmentLine4 like '%Osi%' or posttreatmentLine5 like '%Osi%' or posttreatmentLine6 like '%Osi%' or posttreatmentLine7 like '%Osi%' or posttreatmentLine8 like '%Osi%' or posttreatmentLine9 like '%Osi%')
AND t0.sampleId not like 'CORE%';


