select sampleId, count(*)
from somaticVariant where filter = "PASS" and worstCodingEffect = "MISSENSE" and sampleId in ('XXX')
group by 1;