select t0.sampleId1, t0.driver, t0.qcStatus, t0.dateBiopsy1, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine1PdDate) as daysPdAfterBiopsy,
pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine2PdDate) as daysPdAfterBiopsy, 
pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine3PdDate) as daysPdAfterBiopsy, 
pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine4PdDate) as daysPdAfterBiopsy, 
pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine5PdDate) as daysPdAfterBiopsy, 
pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine6PdDate) as daysPdAfterBiopsy, 
pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine7PdDate) as daysPdAfterBiopsy, 
pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason, dateDiff(t0.dateBiopsy1,pretreatmentLine8PdDate) as daysPdAfterBiopsy
from
(select n.sampleId1, driver, qcStatus, n.dateBiopsy1 from nsclcPairedSample n inner join hmfpatients.purity p on n.sampleId1=p.sampleId inner join nsclcPairedDriver pd on n.sampleId1=pd.sampleId where (qcStatus NOT LIKE '%WARN_LOW_PURITY%' AND qcStatus <> 'FAIL_NO_TUMOR')) t0
left join 
(select sampleId1, pretreatmentLine1, pretreatmentLine1PdDate, pretreatmentLine1StopReason from nsclcPairedSample where preTreatmentLine1 like '%Osi%') t1
on t0.sampleId1=t1.sampleId1
left join 
(select sampleId1, pretreatmentLine2, pretreatmentLine2PdDate, pretreatmentLine2StopReason from nsclcPairedSample where preTreatmentLine2 like '%Osi%') t2
on t0.sampleId1=t2.sampleId1
left join 
(select sampleId1, pretreatmentLine3, pretreatmentLine3PdDate, pretreatmentLine3StopReason from nsclcPairedSample where preTreatmentLine3 like '%Osi%') t3
on t0.sampleId1=t3.sampleId1
left join 
(select sampleId1, pretreatmentLine4, pretreatmentLine4PdDate, pretreatmentLine4StopReason from nsclcPairedSample where preTreatmentLine4 like '%Osi%') t4
on t0.sampleId1=t4.sampleId1
left join 
(select sampleId1, pretreatmentLine5, pretreatmentLine5PdDate, pretreatmentLine5StopReason from nsclcPairedSample where preTreatmentLine5 like '%Osi%') t5
on t0.sampleId1=t5.sampleId1
left join 
(select sampleId1, pretreatmentLine6, pretreatmentLine6PdDate, pretreatmentLine6StopReason from nsclcPairedSample where preTreatmentLine6 like '%Osi%') t6
on t0.sampleId1=t6.sampleId1
left join 
(select sampleId1, pretreatmentLine7, pretreatmentLine7PdDate, pretreatmentLine7StopReason from nsclcPairedSample where preTreatmentLine7 like '%Osi%') t7
on t0.sampleId1=t7.sampleId1
left join 
(select sampleId1, pretreatmentLine8, pretreatmentLine8PdDate, pretreatmentLine8StopReason from nsclcPairedSample where preTreatmentLine8 like '%Osi%') t8
on t0.sampleId1=t8.sampleId1
WHERE (pretreatmentLine1 is not null or pretreatmentLine2 is not null or pretreatmentLine3 is not null or pretreatmentLine4 is not null or pretreatmentLine5 is not null or pretreatmentLine6 is not null or pretreatmentLine7 is not null or pretreatmentLine8 is not null) 
AND t0.sampleId1 not like 'CORE%';
