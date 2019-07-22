load("/Users/jon/hmf/analysis/mantaVgridss/rawData.RData")
rm(manta, reportableDels, strelka)
gridss = gridss %>% filter(sampleId == 'CPCT02450014T')

save(gridss, file = "/Users/jon/hmf/analysis/mantaVgridss/CPCT02450014T.RData")






library(BSgenome.Hsapiens.UCSC.hg19)
load(file = "~/hmf/analysis/mantaVgridss/CPCT02450014T.RData")

primeLength = 1000

gridss = gridss %>% filter(sampleId == 'CPCT02450014T', type == 'DUP') %>%
  mutate(endPosition = as.numeric(endPosition), length = endPosition - startPosition + 1 + nchar(insertSequence)) %>%
  filter(length >= 30, length <= 50) %>%
  select(sampleId, startChromosome, startPosition, endPosition, type, length, startHomologySequence, insertSequence, endHomologySequence) %>%
  mutate(
    fivePrimeSequence = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), startPosition - primeLength, startPosition-1)),
    dupSequence = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), startPosition, endPosition)),
    threePrimeSequence = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), endPosition+1, endPosition + primeLength)),
    ref = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), startPosition - primeLength, endPosition + primeLength)),
    alt = paste0(fivePrimeSequence, dupSequence, insertSequence, dupSequence, threePrimeSequence),
    match = ref == paste0(fivePrimeSequence, dupSequence, threePrimeSequence))

viewBam = gridss %>%
  mutate(
    cmd = paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02450014R_CPCT02450014T.assembly.bam.sv.bam ", startChromosome, ":" ,startPosition, "-", endPosition, "> verification.sam"),
    cmd2 = paste0("grep ", alt, " verification.sam")
    )

write.csv(gridss, file = "/Users/jon/hmf/analysis/mantaVgridss/CPCT02450014T.tandem.dups.csv")



library(RMySQL)
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
currentProd = dbGetQuery(dbProd, "SELECT * from structuralVariant where sampleId = 'CPCT02450014T' and type = 'DUP' and filter = 'PASS'")
dbDisconnect(dbProd)
rm(dbProd)

new = currentProd %>% mutate(length = endPosition - startPosition + 1) %>%
  filter(length >= 30, length <= 50) %>% select(startChromosome, startPosition, startHomologySequence, length)

old = gridss %>% select(startChromosome, startPosition, startHomologySequence)

combined = full_join(new, old, by = c("startChromosome", "startPosition"), suffix = c(".prod", ".compare")) %>%
  mutate(match = startHomologySequence.prod == startHomologySequence.compare)


