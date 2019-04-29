

# Simulation for skin samples

simSkinMatrixData = read.csv('~/data/sigs/input/snv_skin_sim_sc.csv')
View(simSkinMatrixData[,1:20])

simSkinSampleCounts = matrix_to_sample_counts(simSkinMatrixData,snvBuckets)
View(simSkinSampleCounts)
View(simSkinSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count)) %>% arrange(-SampleTotal))
