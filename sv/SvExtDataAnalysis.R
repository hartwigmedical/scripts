library(dplyr)
library(tidyr)


##############################
# EXTERNAL GROUP DATA ANALYSIS


purplePath='~/data/ext_runs/garvan/garvan001T/purple/' 
linxPath='~/data/ext_runs/garvan/garvan001T/linx/' 
sampleId='garvan001T'

## PURPLE DATA
samplePurity = read.csv(paste0(purplePath,sampleId,'.purple.purity.tsv'),sep='\t')
View(samplePurity)

copyNumber = read.csv(paste0(purplePath,sampleId,'.purple.segment.tsv'),sep='\t')
View(copyNumber)

geneCopyNumber = read.csv(paste0(purplePath,sampleId,'.purple.cnv.gene.tsv'),sep='\t')
View(geneCopyNumber)


## LINX DATA
svs = read.csv(paste0(linxPath,sampleId,'.svs.tsv'),sep='\t')
svData = read.csv(paste0(linxPath,sampleId,'.linx.svs.tsv'),sep='\t')
clusters = read.csv(paste0(linxPath,sampleId,'.linx.clusters.tsv'),sep='\t')
tiLinks = read.csv(paste0(linxPath,sampleId,'.linx.links.tsv'),sep='\t')
fusions = read.csv(paste0(linxPath,sampleId,'.linx.fusion.tsv'),sep='\t')
transBreakends = read.csv(paste0(linxPath,sampleId,'.linx.breakend.tsv'),sep='\t')
viralInserts = read.csv(paste0(linxPath,sampleId,'.linx.viral_inserts.tsv'),sep='\t')
driverCatalog = read.csv(paste0(linxPath,sampleId,'.driver.catalog.tsv'),sep='\t')
linxDrivers = read.csv(paste0(linxPath,sampleId,'.linx.drivers.tsv'),sep='\t')

# consolidated SV view

linxFusions = merge(fusions,
                    transBreakends %>% filter(IsUpstream=='true') %>% 
                      select(FivePrimeBreakendId=Id,GeneUp=Gene,SvIdUp=SvId,IsStartUp=IsStart,TransIdUp=TranscriptId,CanonicalUp=Canonical,
                             DisruptiveUp=Disruptive,RegionTypeUp=RegionType,CodingContextUp=CodingContext,BiotypeUp=Biotype,
                             ExonBasePhaseUp=ExonBasePhase,NextSpliceRankUp=NextSpliceRank,NextSplicePhaseUp=NextSplicePHase,NextSpliceDistanceUp=NextSpliceDistance),
                    by='FivePrimeBreakendId',all.x=T)

linxFusions = merge(linxFusions,
                    transBreakends %>% filter(IsUpstream=='false') %>% 
                      select(ThreePrimeBreakendId=Id,GeneDown=Gene,SvIdDown=SvId,IsStartDown=IsStart,TransIdDown=TranscriptId,CanonicalDown=Canonical,
                             DisruptiveDown=Disruptive,RegionTypeDown=RegionType,CodingContextDown=CodingContext,BiotypeDown=Biotype,
                             ExonBasePhaseDown=ExonBasePhase,NextSpliceRankDown=NextSpliceRank,NextSplicePhaseDown=NextSplicePHase,NextSpliceDistanceDown=NextSpliceDistance),
                    by='ThreePrimeBreakendId',all.x=T)

View(linxFusions)

linxDrivers = merge(linxDrivers,driverCatalog,by='gene',all=T)
View(linxDrivers)

View(clusters)
View(svData)


## CLEAN-UP
rm(samplePurity)
rm(copyNumber)
rm(geneCopyNumber)
rm(svData)
rm(clusters)
rm(tiLinks)
rm(fusions)
rm(breakends)
rm(linkxFusions)
rm(linxDrivers)



######
## SPECIFIC GROUPS


## CMRI
egSvData = read.csv('~/data/ext_runs/cmri/LNX_SVS.csv')
View(egSvData)
sgDrivers = read.csv('~/data/ext_runs/cmri/LNX_DRIVERS.csv')
View(sgDrivers)

View(egSvData %>% filter(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n'))))
View(egSvData %>% filter(Type=='SGL'&(grepl('GGGTTA',InsertSeq)|grepl('CCCTAA',InsertSeq))))
                               













