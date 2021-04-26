# total count
select count(*) AS samples, 'all' AS category from clinical where sampleId like '%WIDE%'

union

# samples with primary tumor location
select count(*) AS samples, 'with curated primary tumor location' AS category from clinical where sampleId like '%WIDE%' AND NOT(isnull(primaryTumorLocation))

union

# samples with biopsy data
select count(*) AS samples, 'with matched biopsy' AS category from clinical where sampleId like '%WIDE%' AND NOT(isnull(biopsyDate))

union

# samples with treatments given in other hospital than NKI
select count(distinct sampleId, treatmentGiven) AS samples, 'treatment given in other hospital' AS category from treatment inner join sample on treatment.patientId=sample.patientId where sampleId like '%WIDE%' AND treatmentGiven = "no"

union

# samples with treatments given in NKI
select count(distinct sampleId, treatmentGiven) AS samples, 'treatment given at NKI' AS category from treatment inner join sample on treatment.patientId=sample.patientId where sampleId like '%WIDE%' AND treatmentGiven = "yes"

union
# samples with treatment responses
select count(distinct sampleId) AS samples, 'with treatment responses' AS category from treatmentResponse inner join sample on treatmentResponse.patientId=sample.patientId where sampleId like '%WIDE%' AND measurementDone = "yes";

