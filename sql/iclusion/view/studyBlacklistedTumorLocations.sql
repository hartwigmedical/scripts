CREATE OR REPLACE VIEW studyBlacklistedTumorLocation AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, blacklistedTumorLocation
from study
left join blacklistedTumorLocation
on study.id= blacklistedTumorLocation.blacklistedTumorLocationId;