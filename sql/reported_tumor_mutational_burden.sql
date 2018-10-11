select count(*)/2859 as TMB, if (count(*)/2859 > 10, "High", "Low") as status
from somaticVariant where filter = 'PASS' and sampleId in ('XXX');
