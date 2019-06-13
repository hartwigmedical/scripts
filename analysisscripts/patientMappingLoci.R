######################## Generate BED file of AMBER heterozygous locations to use
oldHetLocations = read.table(file = "/Users/jon/hmf/resources/CytoScanHD_hg19_SNPs_sorted.bed", sep = "\t", header = F, stringsAsFactors = F) %>% select(chr = V1, start = V2, end = V3) %>% mutate(old = T)
newHetLocations = read.table(file = "/Users/jon/hmf/resources/GermlineHetPon.hg19.bed", sep = "\t", header = F, stringsAsFactors = F)  %>% select(chr = V1, start = V2, end = V3) %>% mutate(new = T)
combinedHetLocations = inner_join(oldHetLocations, newHetLocations, by = c("chr", "start", "end"))
patientMappingBed = combinedHetLocations %>% filter(chr != 'X') %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start) %>%
    select(-old, -new) %>%
    mutate(distance = lead(start) - start, source = "Amber") %>%
    filter(distance < 0 | distance > 200000) %>%
    select(chr, start, end)

write.table(patientMappingBed, file = "~/hmf/resources/GermlineSnp.hg19.bed", sep = "\t", row.names = F, col.names = F)
