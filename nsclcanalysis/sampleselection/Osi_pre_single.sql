SELECT t0.sampleId, t0.driver, t0.qcStatus, t0.dateBiopsy, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine1PdDate) as daysPdAfterBiopsy, 
pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine2PdDate) as daysPdAfterBiopsy, 
pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine3PdDate) as daysPdAfterBiopsy, 
pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine4PdDate) as daysPdAfterBiopsy, 
pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine5PdDate) as daysPdAfterBiopsy, 
pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine6PdDate) as daysPdAfterBiopsy, 
pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine7PdDate) as daysPdAfterBiopsy, 
pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(t0.dateBiopsy,pretreatmentLine8PdDate) as daysPdAfterBiopsy
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
WHERE (pretreatmentLine1 is not null or pretreatmentLine2 is not null or pretreatmentLine3 is not null or pretreatmentLine4 is not null or pretreatmentLine5 is not null or pretreatmentLine6 is not null or pretreatmentLine7 is not null or pretreatmentLine8 is not null)  
AND t0.sampleId not like 'CORE%';






