SELECT * FROM germlineVariant
WHERE (
(pathogenic = 'UNANNOTATED' AND codingEffect IN ('NONSENSE_OR_FRAMESHIFT','SPLICE') AND reported = 0)						                    #also check for filter here
OR (pathogenic IN ('CLINVAR_PATHOGENIC','CLINVAR_LIKELY_PATHOGENIC','WHITE_LIST') AND reported = 0)
OR (pathogenic = 'CLINVAR_CONFLICTING' AND filter = 'PASS')
OR (gene IN ('BMPR1A','EPCAM','FH','FLCN','MUTYH','NTHL1','POLD1','SDHA','SDHAF2','SDHB','SDHC','SDHD') AND biallelic = 0 AND reported = 1)			#these genes are not in gene panel AND require biallelic reporting. We should check somaticVariant in case the WT allele is not already lost.
OR (gene IN ('NTHL1','MUTYH') AND reported = 1 AND refStatus = 'HOM')                                                                                   #these genes should only be marked with patient notification in case of refStatus = homozygous
OR (gene IN ('CDK4','RET','MET','KIT') AND reported = 1 AND pathogenic = 'UNANNOTATED' AND codingEffect IN ('NONSENSE_OR_FRAMESHIFT','SPLICE')))                          #are unannotated frameshift/nonsense/splice in 'oncogenes' pathogenic?
AND sampleId IN ('X');
