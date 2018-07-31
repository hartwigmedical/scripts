# LISC: curation primary tumor location
select * from clinicalFindings where message ='failed to curate primary tumor location';

# LISC: curation pre-treatment
select * from clinicalFindings where message ='Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous' and level ='preTreatmentCuration';

# LISC: curation treatment
select * from clinicalFindings where message ='Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous' and level ='treatmentCuration';

# KODU: unused search terms.
select * from clinicalFindings where message ='tumor location search term not used' and level='tumorLocationCuration';
select * from clinicalFindings where message ='Treatment search term not used' and level='treatmentCuration' ;

select * from clinicalFindings where message like '%Matched drugs are based on less than 90% of search term.%';