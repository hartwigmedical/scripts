

# Simulation for skin samples

simSkinMatrixData = read.csv('~/data/sigs/input/snv_skin_sim_sc.csv')
View(simSkinMatrixData[,1:20])

simSkinSampleCounts = matrix_to_sample_counts(simSkinMatrixData,snvBuckets)
View(simSkinSampleCounts)
View(simSkinSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count)) %>% arrange(-SampleTotal))

cosmicSigs = read.csv('~/dev/nmf/snv_cosmic_sigs.csv')
View(cosmicSigs)

simSkinContributions = read.csv('~/data/sigs/input/snv_skin_sim_contributions.csv')
View(simSkinContributions[,1:10])
sigNames30 = get_signame_list(30,T)

tmp = cbind(sigNames30,simSkinContributions)
View(tmp[,1:10])
simSigCounts = gather(tmp, "SampleId", "Count", 2:ncol(tmp)) %>% select(SampleId,SigName=sigNames30,Count)

simSigCounts = merge(simSigCounts,simSigCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count)),by='SampleId',all.x=T)
simSigCounts$SigPercent = round(simSigCounts$Count/simSigCounts$SampleTotal,3)

View(simSigCounts)
View(simSigCounts %>% filter(Count>0))
View(simSigCounts %>% group_by(SigName) %>% summarise(SampleCount=sum(Count>0),
                                                      SampleTotal=sum(Count),
                                                      MedianCount=median(Count),
                                                      MaxCount=max(Count)) %>% filter(SampleCount>0))

View(simSigCounts %>% filter(SigName=='07'&Count>0) %>% arrange(Count))
View(simSigCounts %>% filter(SigName=='07'&Count==0) %>% arrange(SampleTotal))
View(simSigCounts %>% filter(SigName=='07') %>% arrange(Count))

singleSigSamples = simSigCounts %>% group_by(SampleId) %>% summarise(SigCount=sum(Count>0))
View(singleSigSamples)
View(singleSigSamples %>% group_by(SigCount) %>% count())

simSkinActualContribs = read.csv('~/data/sigs/logs/snv_skin_sim_ba_contribs.csv')
baSigNames = get_signame_list(nrow(simSkinActualContribs),T)
baSigNames = c('1_0', '7_123', '2_197', '7ext_203', '17_202', '7ext_188')
View(simSkinActualContribs[,1:10])
tmp = cbind(baSigNames,simSkinActualContribs)
simActualSigCounts = gather(tmp, "SampleId", "Count", 2:ncol(tmp)) %>% select(SampleId,SigName=baSigNames,Count)
View(simActualSigCounts)
simActualSigCounts = merge(simActualSigCounts,simActualSigCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count)),by='SampleId',all.x=T)
simActualSigCounts$SigPercent = round(simActualSigCounts$Count/simActualSigCounts$SampleTotal,3)
View(simActualSigCounts %>% filter(Count>0))

View(simActualSigCounts %>% group_by(SigName) %>% summarise(SampleCount=sum(Count>0),
                                                      SampleTotal=round(sum(Count),0),
                                                      MedianCount=round(median(Count),0),
                                                      MaxCount=round(max(Count),0)))


