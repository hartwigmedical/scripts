select sampleId,
sum(if(observedNormalRatio>0.9 and observedNormalRatio<1.1,bafCount,0))/sum(if(observedNormalRatio>0.65 and observedNormalRatio<1.35,bafCount,0)) as QCcheck
from copyNumberRegion
where sampleId in ('xxx')
group by 1 order by 2;