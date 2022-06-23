CREATE OR REPLACE VIEW studyBlacklistedTumorLocations AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, blacklistedTumorLocation
from study
left join blacklistedTumorLocations
on study.id= blacklistedTumorLocations.blacklistedTumorLocationId;