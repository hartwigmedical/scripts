CREATE OR REPLACE VIEW studyTumorLocations AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, tumorLocation
from study
inner join tumorLocations
on study.id= tumorLocations.tumorLocationId;
