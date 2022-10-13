CREATE OR REPLACE VIEW studyTumorLocation AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, tumorLocation
from study
left join tumorLocation
on study.id= tumorLocation.tumorLocationId;
