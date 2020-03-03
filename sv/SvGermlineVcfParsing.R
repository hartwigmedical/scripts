library(dplyr)

# create the gene panel
dpGenes = read.csv('~/dev/bachelor/drivers_152_genes.csv')
tsgs = read.csv('~/data/gene_panel.csv') # genePanel DB table where reportablePointMutation = 'tsg'
germlineGenes = rbind(dpGenes %>% select(Gene),tsgs)

# remove duplicates and add geneId
View(germlineGenes %>% group_by(Gene) %>% count)
germlineGenes = germlineGenes %>% group_by(GeneName=Gene) %>% count %>% ungroup() %>% select(GeneName)

ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData %>% filter(GeneName=='APC'))
germlineGenes = merge(germlineGenes,ensemblGeneData %>% select(GeneId,GeneName),by='GeneName',all.x = T)

write.csv(germlineGenes %>% select(GeneId,GeneName),'~/data/sv/germline_gene_panel.csv',row.names = F, quote = F)

germlineGenes = read.csv('~/data/sv/germline_gene_panel.csv')
nrow(germlineGenes)
View(germlineGenes)


# clean-up
rm(germlineSVs)
rm(gmDisruptions)
rm(gmSvDisruptions)


## Reprocessed SVs
germlineSVs = read.csv('~/logs/LNX_GERMLINE_SVS.csv')
nrow(germlineSVs)
View(germlineSVs)

# check for duplicates - currently only 1, not sure why
View(germlineSVs %>% group_by(SampleId,Id) %>% count %>% filter(n>1))
View(germlineSVs %>% filter(Id=='gridss15_19451o'))

# high level stats

View(germlineSVs %>% group_by(Summary="Summary",SampleId) %>% summarise(SvCount=n()) %>%
       group_by(Summary) %>% summarise(Samples=n(),TotalSVs=sum(SvCount),
                                       AvgSvCount=round(mean(SvCount),1),MedianSvCount=median(SvCount),
                                       MinSvCount=min(SvCount),MaxSvCount=max(SvCount)))

# lengths of variants
germlineSVs = germlineSVs %>% mutate(Length=ifelse(Type=='DEL'|Type=='DUP'|Type=='INV',PosEnd-PosStart,0))

print(ggplot(germlineSVs %>% filter(Type=='DEL'|Type=='DUP'|Type=='INV') %>% group_by(Type,Length=2**round(log(Length,2))) %>% count %>%
               spread(Type,n,fill=9), aes(x=Length))
      + geom_line(aes(y=DEL, colour='DEL'))
      + geom_line(aes(y=DUP, colour='DUP'))
      + geom_line(aes(y=INV, colour='INV'))
      + scale_x_log10()
      + labs(title='Germline SV lengths',y='Count'))




## Disruptions
gmDisruptions = read.csv('~/logs/LNX_GERMLINE_DISRUPTIONS.csv')
View(gmDisruptions)
View(gmDisruptions %>% group_by(DisruptionType) %>% count)
View(gmDisruptions %>% group_by(DisruptionType,Type) %>% count %>% spread(Type,n,fill=0))

gmDisruptions = gmDisruptions %>% mutate(IsDisruptive=DisruptionType!='SHARD'&DisruptionType!='RECIP_INV')

gmSvDisruptions = gmDisruptions %>% filter(DisruptionType!='SHARD'&DisruptionType!='RECIP_INV') %>% group_by(SampleId,SvId) %>% summarise(DisruptionCount=n())
View(gmSvDisruptions)

# most common genes
View(gmDisruptions %>% filter(IsDisruptive) %>% group_by(GeneName) %>% count %>% arrange(-n))
View(gmDisruptions %>% filter(IsDisruptive) %>% group_by(GeneName,DisruptionType) %>% count %>% spread(DisruptionType,n,fill=0))


# key germline SNP genes
View(gmDisruptions %>% filter(IsDisruptive&GeneName %in% c('CHEK2','BRCA1','BRCA2','ATM','MUTYH','TP53','STK11','PMS2','RAD50')) %>% 
       group_by(GeneName,DisruptionType) %>% count %>% arrange(-n))



## Repeated SVs
germlineSVs = merge(germlineSVs,gmSvDisruptions %>% select(SampleId,Id=SvId,DisruptionCount),by=c('SampleId','Id'),all.x=T)
repeatedSVs = germlineSVs %>% group_by(Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd) %>% 
  summarise(Count=n(),DisruptedSVs=sum(!is.na(DisruptionCount)&DisruptionCount>0)) %>% filter(Count>=5)
View(repeatedSVs)
View(germlineSVs %>% filter(PosStart==57434663))


View(germlineSVs)



View(gmGeneSVs %>% filter(!(GenesStart %in% germlineGenes$GeneName)&!(GenesEnd %in% germlineGenes$GeneName)))


# all gene-related passing SVs
gmGeneSVs= read.csv('~/data/sv/LNX_GERMLINE_SVS.csv')
View(gmGeneSVs)
View(gmGeneSVs %>% filter(GenePanelOverlaps!=''))
nrow(gmGeneSVs)
rm(gmGeneSVs)

View(gmGeneSVs %>% filter(!(Type %in% c('BND','DEL','INS','DUP','SGL','INV'))))
gmGeneSVs = gmGeneSVs %>% filter(Type %in% c('BND','DEL','INS','DUP','SGL','INV'))


gmGeneSVs = gmGeneSVs %>% mutate(SampleId=stri_replace_all_fixed(SampleId,'.bam',''),
                                 SampleId=stri_replace_all_fixed(SampleId,'_dedup.realigned',''),
                                 SampleId=stri_replace_all_fixed(SampleId,'.sorted',''))



# most common genes
gmGeneSVs = gmGeneSVs %>% mutate(GeneStartCount=stri_count_fixed(GenesStart,';')+1,GeneEndCount=stri_count_fixed(GenesStart,';')+1)
View(gmGeneSVs)
View(gmGeneSVs %>% group_by(GeneStartCount) %>% count)
View(gmGeneSVs %>% filter(GeneStartCount>1) %>% group_by(GenesStart) %>% count)

# likely PON
View(gmGeneSVs %>% group_by(GenesStart,GenesEnd,Type,PosStart,PosEnd) %>% count %>% filter(n>5) %>% arrange(-n))


# gene frequency
genesStart = gmGeneSVs %>% filter(GenesStart!='') %>% mutate(Gene=as.character(GenesStart),
                                                             Chr=as.character(ChrStart),
                                                             SampleId=as.character(SampleId),
                                                             Type=as.character(Type),
                                                             Pos=as.numeric(PosStart)) %>% select(SampleId,Type,Chr,Pos,Gene)

genesEnd = gmGeneSVs %>% filter(GenesEnd!='') %>% mutate(Gene=as.character(GenesEnd),
                                                         Chr=as.character(ChrEnd),
                                                         SampleId=as.character(SampleId),
                                                         Type=as.character(Type),
                                                         Pos=as.numeric(PosEnd)) %>% select(SampleId,Type,Chr,Pos,Gene)

sampleGenes = rbind(genesStart,genesEnd) %>% group_by(SampleId,Gene,Chr,) %>% summarise(Count=1)
View(sampleGenes)
rm(genesStart)
rm(genesEnd)

geneFrequency = sampleGenes %>% group_by(Gene) %>% summarise(SampleCount=n())
View(geneFrequency %>% arrange(-SampleCount))




## DBEUG

# nrow(gmGeneSVs %>% filter(GenesStart=='NRG1'|GenesEnd=='NRG1') %>% group_by(SampleId) %>% count)





## LINX parsing
sampleCsv = read.csv('~/data/runs/CPCT02010944/sv/CPCT02010944R_filtered_germline_svs.csv')
sampleCsv = read.csv('~/data/runs/WIDE01010356T/WIDE01010356R_filtered_germline_svs.csv')
sampleCsv = read.csv('~/data/runs/CPCT02030278T/CPCT02030278R_filtered_germline_svs.csv')
View(sampleCsv)

print(4900*30/60/60)
print(4900*60/60/60)




sampleCsv = read.csv('~/logs/CPCT02010944R_filtered_germline_svs.csv')
View(sampleCsv)
View(sampleCsv %>% filter(is.na(ChrEnd)))

View(sampleCsv %>% filter(GenesStart!=''|GenesEnd!='') %>% group_by(Gene=ifelse(GenesStart!='',as.character(GenesStart),as.character(GenesEnd))) %>% count %>% arrange(-n))

View(sampleCsv %>% filter(ChrStart==17&ChrEnd==17&PosStart<42e6&PosEnd>41e6))
View(sampleCsv %>% filter(ChrStart==17&ChrEnd==17&PosStart>40e6&PosEnd<42e6))




## Replication of Gridss Somatic VCF parsing and filtering

#BiocManager::install("VariantAnnotation")
#BiocManager::install("StructuralVariantAnnotation")
# install.packages('BSgenome.Hsapiens.UCSC.hg19')

# library(VariantAnnotation)
# library(StructuralVariantAnnotation)

raw_vcf = readVcf("~/data/runs/CPCT02010944/CPCT02010944R_CPCT02010944T.gridss.vcf")
raw_vcf = readVcf("~/data/runs/CPCT02010944/CPCT02010944_RT_10K.gridss.vcf")
View(raw_vcf)
rm(raw_vcf)

# Tumor is orderinal #2
tumourordinal = seq(ncol(geno(raw_vcf)$VF))[-1]
normalordinal = seq(ncol(geno(raw_vcf)$VF))[-2]
print(tumourordinal)
print(normalordinal)
print(ncol(geno(raw_vcf)$VF))

# understanding the VCF

# list of available attributes
geno(raw_vcf)

# print the headers
geno(header(raw_vcf))

header(raw_vcf)

# SV ID aka 'names()' eg gridss0_11782o
View(names(full_vcf))

# total Qual scores
totalQuals = qual(raw_vcf)
View(totalQuals)

# access all values for a sample eg NORMAL Qual scores
tumorQuals = vcf.geno$QUAL[,2]
View(tumorQuals)

# ??
vcf.rowRanges = rowRanges(raw_vcf)

View(vcf.alt)


vcf.info = info(raw_vcf)
vcf.alt = CharacterList(alt(raw_vcf))
vcf.qual = qual(raw_vcf)
vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
vcf.geno = geno(raw_vcf)
vcf.ad = vcf.geno$AD
normalAD = unlist(lapply(vcf.ad[ ,1], paste0, collapse=","))
tumorAD = unlist(lapply(vcf.ad[ ,2], paste0, collapse=","))
View(tumorAD[1:10,])
View(vcf.geno$QUAL)

print(gridss.allowable_normal_contamination)


# the fields
# GT ASQ ASRP ASSR BANRP BANRPQ BANSR BANSRQ BAQ BASRP BASSR BQ BSC BSCQ BUM BUMQ BVF CASQ IC IQ QUAL RASQ REF REFPAIR RF RP RPQ SR SRQ VF


#####
## FILTERS

# 1. hard filter variants that are obviously not somatic
# only take records where the normalQual/totalQual < 4 * 0.03
gridss.allowable_normal_contamination=0.03
View(VariantAnnotation::fixed(raw_vcf)$QUAL)
full_vcf = raw_vcf[geno(raw_vcf)$QUAL[,1] / VariantAnnotation::fixed(raw_vcf)$QUAL < 4 * gridss.allowable_normal_contamination]
other_vcf = raw_vcf[geno(raw_vcf)$QUAL[,1] / VariantAnnotation::fixed(raw_vcf)$QUAL > 4 * gridss.allowable_normal_contamination]
View(geno(full_vcf)$QUAL)
View(geno(other_vcf)$QUAL)

# 2. hard filter unpaired breakpoints (caused by inconsistent scoring across the two breakends)
# check that each ID has a value in the attributes field PAIR_ID matching the ID field (top level)
tmp_vcf = full_vcf[is.na(info(full_vcf)$PARID) | info(full_vcf)$PARID %in% names(full_vcf)]
View(names(tmp_vcf))

other_vcf = full_vcf[is.na(info(full_vcf)$PARID)]
View(names(other_vcf))

paired_ids = full_vcf[info(full_vcf)$PARID %in% names(full_vcf)]
View(names(paired_ids))

# unsure how alignment works
full_vcf = align_breakpoints(full_vcf)

info(full_vcf)$BPI_AF = rep("", length(full_vcf))

# create a filters string list to store filter info
filters = rep("", length(full_vcf))
names(filters) = names(full_vcf)
View(filters)

View(names(full_vcf))

# get breakpoint GRanges data

bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
View(bpgr)
View(seqnames(bpgr))

source("~/hmf/repos/scripts/gridss/libgridss.R")


# 3. Single breakend filters
befiltered = gridss_breakend_filter(begr, full_vcf, pon_dir='~/logs/', 1, 2)
filters[names(begr)] = befiltered

filtered = rep("", length(bpgr))

# is short DEL or DUP

isShort = is_short_deldup(bpgr)
View(isShort)
# ihomlen = gridss_inexact_homology_length(bpgr, vcf)
homlen = elementExtract(info(full_vcf)$HOMLEN, 1)
homlen[is.na(homlen)] = 0
View(homlen)

bp_short_deldup = strand(bpgr) != strand(partner(bpgr)) &
  seqnames(bpgr) == seqnames(partner(bpgr)) &
  abs(start(bpgr)-start(partner(bpgr))) < gridss.short_event_size_threshold

View(seqnames(bpgr))
View(start(bpgr))
View(start(bpgr[1,]))
View(strand(bpgr[1,]))


# 3a. PON filter...
filtered = .addFilter(filtered, "PON", gridss_overlaps_breakend_pon(gr, pon_dir))

# 3b. Strand bias - for DELs and DUPs < 1K 
filtered = .addFilter(filtered, "strand_bias", !is.na(strandbias) & isShort & strandbias > gridss.max_allowable_short_event_strand_bias)

# 3c. AF
tmp = gridss_bp_af(bpgr, full_vcf, 2) < gridss.min_af
View(tmp)

# .gridss_af = function(gr, vcf, ordinal, includeRefPair, no_coverage_af=0, includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE) {
  
# 2. Imprecise

is_imprecise = !(is.na(info(full_vcf)$IMPRECISE) | !info(full_vcf)$IMPRECISE) |
  !((!is.na(info(full_vcf)$PARID) & info(full_vcf)$ASSR + info(full_vcf)$SR + info(full_vcf)$IC > 0) |
      (is.na(info(full_vcf)$PARID) & info(full_vcf)$BASSR + info(full_vcf)$BSC > 0))
filters[names(full_vcf)[is_imprecise]] = paste0(filters[names(full_vcf)[is_imprecise]], ";imprecise")
View(filters)

# 3. “NO_ASRP”: Filter single breakend variants without an assembly containing at least one discordant read pair

basrp = info(full_vcf)$BASRP
View(basrp %>% filter(V1>0))
filtered = .addFilter(filtered, "NO_ASRP", i$BASRP == 0)

# 4. Filter single breakends with a poly-C or poly-G run at least 16bp in the breakend sequence

View(bpgr$insSeq)

filtered = .addFilter(filtered, "LongPolyC", str_detect(gr$insSeq, "CCCCCCCCCCCCCCCC") | str_detect(gr$insSeq, "GGGGGGGGGGGGGGGG"))

  
genotype = geno(full_vcf[names(bpgr)])
View(genotype$QUAL[,2])


View(genotype)

afSum = .genosum(genotype[['af']], 2)
g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], 2) } else { genotype[[field]] } })
View(g)

print(length(bpgr))
names(g) <- names(genotype)
support = rep(0, length(bpgr))
support = support + g$VF
support = support + g$BVF
View(support)
View(g$REF)
View(g$REFPAIR)
View(g$BVF)
View(g$VF)
vf_af = (VF + BVF) / (VF + BVF + REF + REFPAIR)
vf_af[is.nan(vf_af)] = no_coverage_af
View(vf_af)
print((g[1,])$REF)

print(107 / (107 + 83 + 91))
print(55 / (55 + 83 + 91))


bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
View(bpgr)
bpgr$af = round(gridss_bp_af(bpgr, full_vcf, 2), 5)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
if (length(bpgr) > 0) {
  info(full_vcf[names(bpgr)])$BPI_AF = bpgr$af_str
}

bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, bsgenome=refgenome, pon_dir='`/logs/', 1, 2)
filters[names(bpgr)] = bpfiltered




vcf_data_frame<- function(vcf) 
{
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.qual = qual(vcf)
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  vcf.geno = geno(vcf)
  vcf.ad = vcf.geno$AD
  normalAD = unlist(lapply(vcf.ad[ ,1], paste0, collapse=","))
  tumorAD = unlist(lapply(vcf.ad[ ,2], paste0, collapse=","))
  
  vcf.df = data.frame(
    chromosome = seqnames(vcf), 
    pos = start(vcf), 
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),
    qual = as.numeric(vcf.qual),
    strelka = vcf.info$STRELKA,
    af = unlist(vcf.info$AF),
    mappability = vcf.info$MAPPABILITY,
    germlinePonCount = vcf.info$GERMLINE_PON_COUNT,
    somaticPonCount = vcf.info$SOMATIC_PON_COUNT,
    normalAD = normalAD,
    tumorAD = tumorAD )
  
  vcf.df[is.na(vcf.df)] <- 0
  
  return (vcf.df)
}

vcf = readVcf("~/data/runs/CPCT02010944/CPCT02010944R_CPCT02010944T.gridss.pass_only.vcf")
View(vcf)
vcfData = vcf_data_frame(vcf)

vcf = readVcf("~/hmf/analyses/sage/colo829.sage.final.vcf.gz")
vcf.df = vcf_data_frame(vcf)

sage = vcf.df %>% separate(normalAD,c('germlineRefSupport','germlineAltSupport'),sep=',') %>% separate(tumorAD,c('tumorRefSupport','tumorAltSupport'),sep=',')  %>% 
  mutate(tumorRefSupport=as.numeric(tumorRefSupport),tumorAltSupport=as.numeric(tumorAltSupport),germlineRefSupport=as.numeric(germlineRefSupport),germlineAltSupport=as.numeric(germlineAltSupport))
View(sage %>% filter(strelka==T))
View(sage %>% filter(germlineAltSupport<0.05*tumorAltSupport,tumorAltSupport>=3,germlineRefSupport<100,germlineRefSupport>10,germlinePonCount<6,somaticPonCount<5,strelka==F,tumorAltSupport>0.03*tumorRefSupport,qual>120))
