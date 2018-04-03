select patientId, studyEvent, form, itemGroup, item, itemValue, fieldName from drupEcrf where patientId = 'DRUP01010001' and item in ('FLD. BAS.ICDTC','FLD.CSF.GEN',
'FLD.CSF.YOB','FLD.BAS.BASTTYP','FLD.BAS.BASTTOSP','FLD.BIOPT.TBDAT','FLD.BIOPT.TBSITE','FLD.BIOPT.TBLOC',
'FLD.BIOPT.TBTAKEN','FLD.REG.INST','FLD.EOT.DEATHDTC')