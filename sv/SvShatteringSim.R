

simResults = read.csv('~/logs/LNX_SIM_RESULTS.csv')
nrow(simResults)
View(simResults)
View(head(simResults,10000))

View(simResults %>% group_by(SegCount) %>% summarise(AvgLinkedPerc=round(sum(SegsLinked/SegCount*RepeatCount)/first(TestIterations),3)))

# graph of segments linked

View(simResults %>% group_by(SegCount) %>% summarise(AvgLinkedPerc=round(sum(SegsLinked/SegCount*RepeatCount)/first(TestIterations),3)))

print(ggplot(simResults %>% filter(SegCount<=10) %>% group_by(SegCount,SegsLinked) %>% summarise(Count=sum(RepeatCount)),aes(x=SegsLinked, y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~SegCount)
      + labs(title = "Shattering Simulation - number of segments relinked"))

print(ggplot(simResults %>% filter(SegCount==20) %>% group_by(SegCount,SegsLinked) %>% summarise(Count=sum(RepeatCount)),aes(x=SegsLinked))
      # + geom_line(stat='identity',colour='black')
      + stat_ecdf(geom='point')
      + facet_wrap(~SegCount)
      + labs(title = "Shattering Simulation - number of segments relinked"))

simResultsUngrp = read.csv('~/logs/LNX_SIM_RESULTS.csv')
nrow(simResultsUngrp)
View(simResultsUngrp %>% filter(SegCount==50))

View(simResultsUngrp %>% group_by(SegCount) %>% summarise(AvgLinkedPerc=round(sum(SegsLinked/SegCount)/first(TestIterations),3)))

print(ggplot(simResults %>% filter(SegCount>=40), aes(x=SegsLinked)) 
     + stat_ecdf(aes(group=SegCount,colour=SegCount),geom='point')
     + facet_wrap(~SegCount))

print(ggplot(simResults %>% filter(SegCount==50), aes(x=SegsLinked)) 
      + stat_ecdf(aes(group=SegCount,colour=SegCount),geom='point'))

print(ggplot(simResultsUngrp %>% filter(SegCount==50), aes(x=AdjacentSegs)) 
      + stat_ecdf(aes(group=SegCount,colour=SegCount),geom='point'))

print(ggplot(simResults %>% filter(SegCount<=10), aes(x=SegsLinked)) 
      + stat_ecdf(geom='point'))



# simResults = merge(simResults,simResults %>% group_by(SegCount) %>% summarise(TotalOutcomes=sum(RepeatCount)),by='SegCount',all.x=T)
simResults = simResults %>% mutate(RetainedPerc=round(SegsLinked/TestIterations,6))

# View(simResults %>% group_by(SegCount,SegsLinked) %>% summarise(Outcomes=n(),Percent=round(n()/first(TotalOutcomes),3)))

View(simResults %>% group_by(SegCount,SegsLinked) %>% summarise(Outcomes=n(),Percent=round(n()/first(TotalOutcomes),3)))


View(simResults %>% filter(SegCount==5) %>% group_by(SegsLinked,ExactRepairs,AdjacentSegs,InfLinks,InfLost) %>% count)

View(simResults %>% filter(SegCount==9&InfLinks==4&InfLost==5))

# statistics
View(simResults %>% group_by(SegCount) %>% summarise(LinkedPerc=round(sum(SegsLinked/SegCount*RepeatCount)/first(TestIterations),3),
                                                     AvgLinked=round(sum(SegsLinked*RepeatCount)/first(TestIterations),1),
                                                     NoLinks=round(sum(ifelse(SegsLinked==0,RepeatCount,0))/first(TestIterations),3),
                                                     AllLinks=round(sum(ifelse(SegsLinked==SegCount,RepeatCount,0))/first(TestIterations),3),
                                                     AllExact=round(sum(ifelse(ExactRepairs==SegCount+1,RepeatCount,0))/first(TestIterations),6),
                                                     AvgInfLinks=round(sum(InfLinks*RepeatCount)/first(TestIterations),1),
                                                     AvgInfLost=round(sum(InfLost*RepeatCount)/first(TestIterations),1),
                                                     AvgInfBreaks=round(sum((InfLinks+InfLost)*RepeatCount)/first(TestIterations),1)))

View(simResults %>% group_by(SegCount,RetainedPerc=round(SegsLinked/SegCount,1)) %>% summarise(Count=sum(RepeatCount)) %>% spread(RetainedPerc,Count,fill=0))
                                                    

print(ggplot(simResults %>% group_by(SegCountGrp=ifelse(SegCount<=10,SegCount,2**round(log(SegCount,2))),
                                     RetainedPerc=round(SegsLinked/SegCount,1)) %>% summarise(Perc=sum(RepeatCount)/first(TestIterations)),aes(x=RetainedPerc,y=Perc))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~SegCountGrp)
      + labs(title = "Shattering Sim - Retained Material as % of Total"))

print(ggplot(simResults %>% filter(SegCount %in% c(3,4,5,6,7,8,9,10)) %>% 
               group_by(SegCount,RetainedPerc=round(SegsLinked/SegCount,1)) %>% summarise(Count=sum(RepeatCount)),aes(x=RetainedPerc,y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~SegCount)
      + labs(title = "Shattering Sim - Retained Material as % of Total"))

print(ggplot(simResults %>% filter(SegCount %in% c(15,20,25,30,35,40,45,50)) %>% 
               group_by(SegCount,RetainedPerc=round(SegsLinked/SegCount,2)) %>% summarise(Count=sum(RepeatCount)),aes(x=RetainedPerc,y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~SegCount)
      + labs(title = "Shattering Sim - Retained Material as % of Total"))

# by inferred segs and lost
simResults = simResults %>% mutate(InferredSegs=InfLinks+InfLost,
                                   InferredRetainedPerc=round(1-(InfLost/InferredSegs),3))

View(simResults %>% filter(SegCount>=10))

print(ggplot(simResults %>% filter(SegCount>=10) %>% group_by(RetainedPerc=round(InferredRetainedPerc,2)) %>% summarise(Count=sum(RepeatCount)),aes(x=RetainedPerc,y=Count))
      + geom_bar(stat='identity',colour='black')
      + labs(title = "Shattering Sim - Inferred Segments Retained %"))

print(ggplot(simResults %>% filter(InferredSegs %in% c(3,4,5,6,7,8,9,10)) %>% 
               group_by(InferredSegs,RetainedPerc=round(InferredRetainedPerc,1)) %>% summarise(Count=sum(RepeatCount)),aes(x=RetainedPerc,y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~InferredSegs)
      + labs(title = "Shattering Sim - Retained Material as % of Total"))

print(ggplot(simResults %>% filter(InferredSegs %in% c(15,20,25,30,35,40,45,50)) %>% 
               group_by(InferredSegs,RetainedPerc=round(InferredRetainedPerc,2)) %>% summarise(Count=sum(RepeatCount)),aes(x=RetainedPerc,y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~InferredSegs)
      + labs(title = "Shattering Sim - Retained Material as % of Total"))


# grouped data equivalent with median values
View(simResults %>% group_by(SegCount) %>% summarise(Tests=first(TestIterations),
                                                     AvgLinked=round(mean(SegsLinked),3),
                                                     MedLinked=round(median(SegsLinked),3),
                                                     NoLinks=round(sum(ifelse(SegsLinked==0,RepeatCount,0))/first(TestIterations),3),
                                                     AllLinks=round(sum(SegsLinked==SegCount)/first(TestIterations),3),
                                                     AllExact=round(sum(ExactRepairs==SegCount+1)/first(TestIterations),6),
                                                     MedInfLinks=round(median(InfLinks),1),
                                                     MedInfLost=round(median(InfLost),1),
                                                     MedInfBreaks=round(median(InfLinks+InfLost),1)))

# DEBUG


# 2-segment theoretical outcomes
print(48/85) # 52%
print(13/85) # 15%
print(24/85) # 28%


tmp = read.csv('~/logs/LNX_SIM_RESULTS.csv')
View(tmp)

# tmp = merge(tmp,tmp %>% group_by(SegCount) %>% summarise(Outcomes=n()),by='SegCount',all.x=T)
# View(tmp %>% group_by(SegCount,SegsLinked) %>% summarise(RetainedPerc=round(n()/first(Outcomes),4)))

tmp20Random = tmp
tmp20Select100K = tmp

View(tmp %>% group_by(SegCount,SegsLinked) %>% summarise(RetainedPerc=round(n()/first(TestIterations),4)))
View(tmp20Random %>% group_by(SegCount,SegsLinked) %>% summarise(RetainedPerc=round(n()/first(TestIterations),4)))
View(tmp20Select100K %>% group_by(SegCount,SegsLinked) %>% summarise(RetainedPerc=round(n()/first(TestIterations),4)))

View(tmp %>% group_by(SegCount) %>% summarise(RetainedPerc=round(sum(SegsLinked/SegCount)/first(TestIterations),4)))

View(tmp %>% group_by(LinkStr,SegsLinked) %>% count)
View(tmp %>% group_by(LinkStr,SegsLinked) %>% count %>% group_by(SegsLinked) %>% summarise(LinkStr=n(),Count=sum(n)))

tmp = tmp %>% mutate(StrLen=stri_length(LinkStr),LinkStrExtra=LinkStr)
tmp = tmp %>% separate(LinkStrExtra,c('Rnd1','Rnd2','Rnd3'),sep=';')
# tmp = tmp %>% separate(LinkStrExtra,c('Rnd1','Rnd2'),sep=';')
View(tmp)

View(tmp %>% group_by(StrLen,SegsLinked) %>% count)
View(tmp %>% group_by(StrLen,TestIterations) %>% count %>% mutate(Percent=round(n/TestIterations,3)))
View(tmp %>% group_by(Rnd1) %>% count)
View(tmp %>% group_by(Rnd1,LinkStr) %>% count)
View(tmp %>% group_by(Rnd2) %>% count)
View(tmp %>% group_by(Rnd3) %>% count)

View(tmp %>% filter(SegsLinked==2) %>% group_by(LinkStr) %>% count)

View(tmp %>% filter(Rnd1!='0:e-3:s') %>% group_by(Rnd1,Rnd2) %>% count %>% group_by(Rnd1) %>% summarise(NumRnd2=n(),Count=sum(n)))


nrow(tmp %>% filter(LinkStr=='0:e-3:s'))

tmpGrp = read.csv('~/logs/LNX_SIM_RESULTS.csv')
View(tmpGrp)

View(tmpGrp %>% group_by(SegCount,SegsLinked) %>% summarise(AvgLinked=round(sum(RepeatCount)/first(TestIterations),4)))
                                                 

View(tmpGrp %>% group_by(SegCount) %>% summarise(Tests=first(TestIterations),
                                                     LinkedPerc=round(sum(SegsLinked/SegCount*RepeatCount)/first(TestIterations),3),
                                                     AvgLinked=round(sum(SegsLinked*RepeatCount)/first(TestIterations),1),
                                                     NoLinks=round(sum(ifelse(SegsLinked==0,RepeatCount,0))/first(TestIterations),3),
                                                     AllLinks=round(sum(ifelse(SegsLinked==SegCount,RepeatCount,0))/first(TestIterations),3),
                                                     AllExact=round(sum(ifelse(ExactRepairs==SegCount+1,RepeatCount,0))/first(TestIterations),6),
                                                     AvgInfLinks=round(sum(InfLinks*RepeatCount)/first(TestIterations),1),
                                                     AvgInfLost=round(sum(InfLost*RepeatCount)/first(TestIterations),1),
                                                     AvgInfBreaks=round(sum((InfLinks+InfLost)*RepeatCount)/first(TestIterations),1)))

View(tmp %>% group_by(SegCount,SegsLinked) %>% summarise(RetainedPerc=round(n()/first(TestIterations),4)))
View(tmp %>% group_by(SegCount) %>% summarise(RetainedPerc=round(sum(SegsLinked/SegCount)/first(TestIterations),4)))


View(simResults %>% group_by(SegCount,SegsLinked,ExactMatches,AdjacentSegs) %>% summarise(Outcomes=n(),
                                                                                          Percent=round(n()/first(TotalOutcomes),3)))
