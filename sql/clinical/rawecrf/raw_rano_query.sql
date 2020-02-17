select patients.patientId,
RNTMRNYN, RNTMDT, RNTMTNR, RNTMTST, RNTMTLC, RNTMTMTHD, RNTMTTRGT1, RNTMTTRGT2, RNTMNWYN, RNTSER,
RNTSLC, RNTMSMTTL, RNDUMTMSUMTOTAL, RNTMTRGTT, RNTMNTNR, RNTMNTST, RNTMNTLC, RNTMNTMTHD, RNTMNTVL, RNTT2,
TERNTB, TERNCFB, TERNCFN, TERNOR, RNDNAD, RNCORTU, RNVRLL, RNRSPNTL, RNRSPTL
from
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as RNTMRNYN
    from cpctEcrf where item = 'FLD.TMRANO.RNTMRNYN' group by patientId) rntmrnyn
on patients.patientId = rntmrnyn.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMDT
     from cpctEcrf where item ='FLD.TMRANO.RNTMDT' group by patientId) rntmdt
on patients.patientId = rntmdt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTNR
     from cpctEcrf where item ='FLD.TMRANO.RNTMTNR' group by patientId) rntmtnr
on patients.patientId = rntmtnr.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTST
     from cpctEcrf where item ='FLD.TMRANO.RNTMTST' group by patientId) rntmtst
on patients.patientId = rntmtst.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTLC
     from cpctEcrf where item ='FLD.TMRANO.RNTMTLC' group by patientId) rntmtlc
on patients.patientId = rntmtlc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTMTHD
    from cpctEcrf where item ='FLD.TMRANO.RNTMTMTHD' group by patientId) rntmtmthd
on patients.patientId = rntmtmthd.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTTRGT1
     from cpctEcrf where item ='FLD.TMRANO.RNTMTTRGT1' group by patientId) rntmttrgt1
on patients.patientId = rntmttrgt1.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTTRGT2
     from cpctEcrf where item ='FLD.TMRANO.RNTMTTRGT2' group by patientId) rntmttrgt2
on patients.patientId = rntmttrgt2.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNWYN
     from cpctEcrf where item ='FLD.TMRANO.RNTMNWYN' group by patientId) rntmnwyn
on patients.patientId = rntmnwyn.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTSER
     from cpctEcrf where item ='FLD.TMRANO.RNTSER' group by patientId) rntser
on patients.patientId = rntser.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTSLC
     from cpctEcrf where item ='FLD.TMRANO.RNTSLC' group by patientId) rntslc
on patients.patientId = rntslc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMSMTTL
     from cpctEcrf where item ='FLD.TMRANO.RNTMSMTTL' group by patientId) rntmsmttl
on patients.patientId = rntmsmttl.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNDUMTMSUMTOTAL
     from cpctEcrf where item ='FLD.TMRANO.RNDUMTMSUMTOTAL' group by patientId) rndumtmsumtotal
on patients.patientId = rndumtmsumtotal.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMTRGTT
     from cpctEcrf where item ='FLD.TMRANO.RNTMTRGTT' group by patientId) rntmtrgtt
on patients.patientId = rntmtrgtt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNTNR
     from cpctEcrf where item ='FLD.TMRANO.RNTMNTNR' group by patientId) rntmntnr
on patients.patientId = rntmntnr.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNTST
     from cpctEcrf where item ='FLD.TMRANO.RNTMNTST' group by patientId) rntmntst
on patients.patientId = rntmntst.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNTLC
     from cpctEcrf where item ='FLD.TMRANO.RNTMNTLC' group by patientId) rntmntlc
on patients.patientId = rntmntlc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNTMTHD
     from cpctEcrf where item ='FLD.TMRANO.RNTMNTMTHD' group by patientId) rntmntmthd
on patients.patientId = rntmntmthd.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTMNTVL
     from cpctEcrf where item ='FLD.TMRANO.RNTMNTVL' group by patientId) rntmntvl
on patients.patientId = rntmntvl.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNTT2
     from cpctEcrf where item ='FLD.TMRANO.RNTT2' group by patientId) rntt2
on patients.patientId = rntt2.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TERNTB
     from cpctEcrf where item ='FLD.TMRANO.TERNTB' group by patientId) terntb
on patients.patientId = terntb.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TERNCFB
     from cpctEcrf where item ='FLD.TMRANO.TERNCFB' group by patientId) terncbf
on patients.patientId = terncbf.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TERNCFN
     from cpctEcrf where item ='FLD.TMRANO.TERNCFN' group by patientId) terncfn
on patients.patientId = terncfn.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TERNOR
     from cpctEcrf where item ='FLD.TMRANO.TERNOR' group by patientId) ternor
on patients.patientId = ternor.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNDNAD
     from cpctEcrf where item ='FLD.TMRANO.RNDNAD' group by patientId) rndnad
on patients.patientId = rndnad.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNCORTU
     from cpctEcrf where item ='FLD.TMRANO.RNCORTU' group by patientId) rncortu
on patients.patientId = rncortu.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNVRLL
     from cpctEcrf where item ='FLD.TMRANO.RNVRLL' group by patientId) rnvrll
on patients.patientId = rnvrll.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNRSPTL
     from cpctEcrf where item ='FLD.TMRANO.RNRSPTL' group by patientId) rnrsplt
on patients.patientId = rnrsplt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as RNRSPNTL
     from cpctEcrf where item ='FLD.TMRANO.RNRSPNTL' group by patientId) rnrspntl
on patients.patientId = rnrspntl.patientId
where patients.patientId in (select distinct patientId from clinical) order by patients.patientId;