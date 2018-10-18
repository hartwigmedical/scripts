
gridssSvData = read.csv('~/logs/CLUSTER_GRIDSS.csv')
nrow(gridssSvData)
View(gridssSvData)

gridssSamples = svData %>% group_by(SampleId) %>% summarise(Count=n())
View(gridssSamples)
# write.csv(gridssSvData, "~/logs/r_output/sv_gridss.csv", row.names=F, quote=F)


# filter to match the current set of GRIDSS samples
svGridssMatched = svData %>% filter(SampleId %in% gridssSamples$SampleId)
View(svGridssMatched %>% group_by(SampleId) %>% summarise(Count=n()))

# reverse filter to ensure GRIDSS doesn't have more samples
gridssSvData = gridssSvData %>% filter(SampleId %in% svGridssMatched$SampleId)
nrow(gridssSvData %>% group_by(SampleId) %>% count())
nrow(svGridssMatched %>% group_by(SampleId) %>% count())

svGridssMatched = within(svGridssMatched, rm(IsResolved))
gridssSvData = within(gridssSvData, rm(Imprecise))
svGridssMatched = within(svGridssMatched, rm(MantaPrecise))
ncol(svGridssMatched)
ncol(gridssSvData)


# GRIDSS Assembly vs SV-Analyser matching
gridssSvData = sv_set_common_fields(gridssSvData)
gridssSvData = set_sv_stressed_state(gridssSvData)
gridssSvData$ResolvedType = 'NONE'
gridssSvData = set_sv_non_clustered_types(gridssSvData)
allClusterData = get_sv_clustered_data(gridssSvData)
gridssSvData = set_sv_line_types(gridssSvData, allClusterData)
gridssSvData = set_sv_clustered_types(gridssSvData, allClusterData)

assemblyMatching = gridssSvData %>% filter(IsDB | IsTI | AsmbStart!="" | AsmbEnd!="")
assemblyMatching = assemblyMatching %>% filter(Type!='SGL'&ClusterCount<100&ClusterCount>1)
assemblyMatching = assemblyMatching %>% filter(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd)|grepl('dsb',AsmbStart)|grepl('dsb',AsmbEnd))
assemblyMatching = assemblyMatching %>% filter(IsLINE==F&ResolvedType!='LINE_CLUST')
assemblyMatching = assemblyMatching %>% filter(SingleDupBE==F&DoubleDupBE==F)

assemblyOnly = assemblyMatching %>% filter((AsmbMatchStart=="ASMB_ONLY"&grepl('asm',AsmbStart)) | (AsmbMatchEnd=="ASMB_ONLY"&grepl('asm',AsmbEnd)))
View(assemblyOnly)

nrow(gridssSvData %>% filter(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd)|grepl('dsb',AsmbStart)|grepl('dsb',AsmbEnd)))


gridssCohortVariants = gridssCohortVariants %>%  filter(!(type=='DEL'&(InexactHOEnd-InexactHOStart>=6|nchar(Homology)<6)&endPosition-startPosition<1000))

View(assemblyMatching)

# split into start and end to make analysis easier?
assemblyMatchStats = assemblyMatching %>% group_by(AsmbMatchStart, AsmbMatchEnd) %>% summarise(Count=n())
View(assemblyMatchStats)

assemblyStarts = gridssSvData %>% filter(grepl(';asm',AsmbStart))
assemblyStarts$IsEnd = F
assemblyStarts$Assembly = assemblyStarts$AsmbStart
assemblyEnds = gridssSvData %>% filter(grepl(';asm',AsmbEnd))
assemblyEnds$IsEnd = T
assemblyEnds$Assembly = assemblyEnds$AsmbEnd
multipleAssemblies = rbind(assemblyStarts, assemblyEnds)

multipleAssemblies = multipleAssemblies %>% filter(!grepl("dsb",Assembly)&!grepl("bp",Assembly))
multipleAssemblyGrouped = (multipleAssemblies %>% group_by(Assembly) %>% summarise(Type=first(Type), PosStart=first(PosStart), PosEnd=first(PosEnd), Count=n()))
multipleAssemblyGrouped$Length = ifelse(multipleAssemblyGrouped$Type!='BND',multipleAssemblyGrouped$PosEnd-multipleAssemblyGrouped$PosStart,0)
View(multipleAssemblyGrouped %>% filter(Count==1))




View(gridssSvData %>% filter(SampleId=='CPCT02010267TIII'&ClusterId==20))

View(gridssSvData)


svResolvedSummary = (svData %>% group_by(ResolvedType,IsStressed,ClusterSize)
                     %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize))

View(svResolvedSummary)

svGridssResolvedSummary = (gridssSvData %>% group_by(ResolvedType,IsStressed,ClusterSize)
                           %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize))

View(svGridssResolvedSummary)

# nrow(svData %>% filter(SampleId=="CPCT02050018T"))

svGridssMatched$Source = "MANTA"
gridssSvData$Source = "GRIDSS"
svGridsAndManta = rbind(svGridssMatched,gridssSvData)
nrow(svGridsAndManta)

svCompareResolvedSummary = (svGridsAndManta %>% group_by(Source,ResolvedType,ClusterSize)
                            %>% summarise(Count=n()) %>% arrange(Source,ResolvedType,ClusterSize))

View(svCompareResolvedSummary)

svCompareResolvedSum2 = svCompareResolvedSummary %>% spread(Source,Count)
View(svCompareResolvedSum2)
svCompareResolvedSum2[is.na(svCompareResolvedSum2)] = 0
svCompareResolvedSum2$Diff = svCompareResolvedSum2$MANTA - svCompareResolvedSum2$GRIDSS
View(svCompareResolvedSum2)

svCompareResolvedSum3 = (svGridsAndManta %>% group_by(Source,SampleId,ResolvedType,ClusterSize)
                         %>% summarise(Count=n()) %>% arrange(Source,SampleId,ResolvedType,ClusterSize) %>% spread(Source,Count))

svCompareResolvedSum3[is.na(svCompareResolvedSum3)] = 0
svCompareResolvedSum3$Diff = abs(svCompareResolvedSum3$MANTA - svCompareResolvedSum3$GRIDSS)
View(svCompareResolvedSum3)

write.csv(svCompareResolvedSum3, "~/logs/r_output/gridss_vs_manta_resolved_types.csv", quote=F, row.names=F)

# matching SVs where possible
svGridsAndMantaMatched = (svGridsAndManta %>% group_by(SampleId,Type,ChrStart,ChrEnd,PosStart,PosEnd)
                          %>% summarise(Count=n(),
                                        NGCount=sum(Source=="MANTA"),
                                        GCount=sum(Source=="GRIDSS"),
                                        Source1=first(Source),
                                        Source2=last(Source),
                                        Id1=first(Id),
                                        Id2=last(Id)))

svGridsAndMantaMatched = svGridsAndMantaMatched %>% filter(Count==2&NGCount==1&GCount==1)
View(svGridsAndMantaMatched)
exactMatchedIds = svGridsAndMantaMatched %>% ungroup() %>% select(Source1,Source2,Id1,Id2)
exactMatchedIds$MatchType = "Exact"
rowIndex = data.frame(as.numeric(as.character(rownames(exactMatchedIds))))
colnames(rowIndex) <- c("RowIndex")
exactMatchedIds = cbind(exactMatchedIds, rowIndex)

View(exactMatchedIds)

svGridsAndManta$MatchType = "None"
svGridsAndManta = within(svGridsAndManta, rm(MatchType))
View(svGridsAndManta)
svGridsAndManta2 = merge(svGridsAndManta, exactMatchedIds %>% select(Id1,Id2,RowIndex,MatchType), by.x="Id", by.y="Id1", all.x=T)
View(svGridsAndManta2)
svGridsAndManta2 = merge(svGridsAndManta2, exactMatchedIds %>% select(Id1,Id2,RowIndex,MatchType), by.x="Id", by.y="Id2", all.x=T)

svGridsAndManta2$MatchSvId = ifelse(!is.na(svGridsAndManta2$Id1),svGridsAndManta2$Id1,svGridsAndManta2$Id2)
svGridsAndManta2$MatchId = ifelse(!is.na(svGridsAndManta2$RowIndex.x),svGridsAndManta2$RowIndex.x,svGridsAndManta2$RowIndex.y)
svGridsAndManta2$Match = ifelse(!is.na(svGridsAndManta2$MatchId),"Exact","None")
View(svGridsAndManta2)
svGridsAndManta = svGridsAndManta2 %>% select(-RowIndex.x,-RowIndex.y,-Id1,-Id2,-MatchType.x,-MatchType.y)

View(exactMatchedIds %>% group_by(Id1) %>% count())
View(svGridsAndManta %>% group_by(Id) %>% count())


View(svGridsAndManta %>% group_by(MatchType) %>% count())

# match based on proximity
svGridsAndManta$PosApproxStart = round(svGridsAndManta$PosStart/10)*10
svGridsAndManta$PosApproxEnd = round(svGridsAndManta$PosEnd/10)*10

svGnMMatchApprox = (svGridsAndManta %>% filter(Match=="None")
                    %>% group_by(SampleId,Type,ChrStart,ChrEnd,PosApproxStart,PosApproxEnd)
                    %>% summarise(Count=n(),
                                  NGCount=sum(Source=="MANTA"),
                                  GCount=sum(Source=="GRIDSS"),
                                  Source1=first(Source),
                                  Source2=last(Source),
                                  Id1=first(Id),
                                  Id2=last(Id)))

svGnMMatchApprox = svGnMMatchApprox %>% filter(Count==2&NGCount==1&GCount==1)
View(svGnMMatchApprox)

approxMatchedIds = svGnMMatchApprox %>% ungroup() %>% select(Source1,Source2,Id1,Id2)
approxMatchedIds$MatchType = "Approx"
rowIndex = data.frame(as.numeric(as.character(rownames(approxMatchedIds))))
colnames(rowIndex) <- c("RowIndex")
approxMatchedIds = cbind(approxMatchedIds, rowIndex)
View(approxMatchedIds)

svGridsAndManta2 = merge(svGridsAndManta, approxMatchedIds %>% select(Id1,Id2,RowIndex,MatchType), by.x="Id", by.y="Id1", all.x=T)
svGridsAndManta2 = merge(svGridsAndManta2, approxMatchedIds %>% select(Id1,Id2,RowIndex,MatchType), by.x="Id", by.y="Id2", all.x=T)
View(svGridsAndManta2)

svGridsAndManta2$MatchSvId = ifelse(!is.na(svGridsAndManta2$Id1),svGridsAndManta2$Id1,ifelse(!is.na(svGridsAndManta2$Id2), svGridsAndManta2$Id2, svGridsAndManta2$MatchSvId))
svGridsAndManta2$MatchId = ifelse(!is.na(svGridsAndManta2$RowIndex.x),svGridsAndManta2$RowIndex.x,ifelse(!is.na(svGridsAndManta2$RowIndex.y),svGridsAndManta2$RowIndex.y, svGridsAndManta2$MatchId))
svGridsAndManta2$Match = ifelse(svGridsAndManta2$Match=="None"&!is.na(svGridsAndManta2$MatchId),"Approx",svGridsAndManta2$Match)

svGridsAndManta = svGridsAndManta2 %>% select(-RowIndex.x,-RowIndex.y,-Id1,-Id2,-MatchType.x,-MatchType.y)
nrow(svGridsAndManta %>% filter(Match=="Exact"))
nrow(svGridsAndManta %>% filter(Match=="Approx"))
nrow(svGridsAndManta %>% filter(Match=="None"))
View(svGridsAndManta)

save(svGridsAndManta, file="~/logs/r_output/svGridsAndManta.RData")



svGridsAndMantaById = svGridsAndManta %>% group_by(Id) %>% count()
View(svGridsAndMantaById)





