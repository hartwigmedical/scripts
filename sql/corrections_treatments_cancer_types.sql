# curation primary tumor location
select * from clinicalFindings where message ='failed to curate primary tumor location';

# curation pre-treatment
select * from clinicalFindings where
message ='Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous' and level ='preTreatmentCuration';

# curation treatment
select * from clinicalFindings where
message ='Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous' and level ='treatmentCuration';

select * from clinicalFindings where message ='tumor location search term not used' and level='tumorLocationCuration';
select * from clinicalFindings where message ='Treatment search term not used' and level='treatmentCuration' ;
