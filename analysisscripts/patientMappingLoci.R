######################## Generate BED file of AMBER heterozygous locations to use
oldHetLocations = read.table(file = "/Users/jon/hmf/resources/CytoScanHD_hg19_SNPs_sorted.bed", sep = "\t", header = F, stringsAsFactors = F) %>% select(chr = V1, start = V2, end = V3) %>% mutate(old = T)
newHetLocations = read.table(file = "/Users/jon/hmf/resources/GermlineHetPon.hg19.bed", sep = "\t", header = F, stringsAsFactors = F)  %>% select(chr = V1, start = V2, end = V3) %>% mutate(new = T)
combinedHetLocations = inner_join(oldHetLocations, newHetLocations, by = c("chr", "start", "end"))
amberBed = combinedHetLocations %>% filter(chr != 'X') %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start) %>%
    select(-old, -new) %>%
    mutate(distance = lead(start) - start, source = "Amber") %>%
    filter(distance < 0 | distance > 200000) %>%
    select(-distance)

## Add in snp check locations
snpDesign = read.table(file = "~/hmf/repos/scripts/snpcheck/26SNPtaq_design.vcf") %>% select(1, 2) %>% mutate(source = "26SNPtaq_design")
names(snpDesign) <- c("chr", "pos","source")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3978886/
snpProfilingPanel = read.csv(file = '~/Downloads/gm492-S1.csv') %>% select(1, 2, 3)
names(snpProfilingPanel) <- c("Optimized", "chr", "pos")
snpProfilingPanel = snpProfilingPanel %>% filter(Optimized== 1) %>% mutate(source = "PengellyEtAl") %>% select(chr, pos, source)

externalDesign = bind_rows(snpDesign, snpProfilingPanel) %>% mutate(start = pos - 1) %>% select(chr, start, end = pos, source)
verifyExternalDesign = left_join(externalDesign, newHetLocations %>% mutate(Amber = T, chr = as.numeric(chr)), by = c("chr","start", "end"))

patientMappingBed = bind_rows(amberBed, externalDesign) %>%
  arrange(chr, start) %>%
  select(-source) %>%
  distinct()

write.table(patientMappingBed, file = "~/hmf/resources/patientMapping.bed", sep = "\t", row.names = F, col.names = F)
