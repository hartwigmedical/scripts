# curation primary tumor location
select * from clinicalFindings where ecrfItem ="FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS" and message = "Failed to curate primary tumor location.";

# curation treatment
select * from clinicalFindings where ecrfItem ='FLD.PRETHERAPY.SYSTEMICREG' and 
message ='Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous.';
