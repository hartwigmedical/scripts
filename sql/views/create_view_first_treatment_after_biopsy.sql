create view first_treatment_after_biopsy as
select treatment.*
from treatment
         join
     -- If a patient has multiple treatments on the same day, we want to select the one with the lowest id
         (select min(treatment.id) as id, treatment.patientId, treatment.startDate
          from treatment
                   -- If a patient has multiple biopsies, we want to select the one with the lowest treatment date
                   join (select treatment.patientId, min(treatment.startDate) as startDate
                         from treatment
                                  left join biopsy on treatment.biopsyId = biopsy.id
                         -- If a patient has a biopsy, filter out treatments that are before the biopsy
                         where (biopsy.biopsyDate is null or biopsy.biopsyDate <= treatment.startDate)
                         group by patientId) as firstDate
                        on firstDate.patientId = treatment.patientId and firstDate.startDate = treatment.startDate
          group by treatment.patientId, treatment.startDate) as firstId on firstId.id = treatment.id;