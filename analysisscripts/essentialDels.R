library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(scales)

approximate_distance <- function(start, end) {
  difference = end - start;
  if (difference > 1000000) {
    return (paste0(round(difference/1000000), "Mb"))
  }
  
  return (paste0(round(difference/1000), "kb"))
}


create_histogram_data <- function (geneCopyNumbers) {
    minStart = min(geneCopyNumbers$start)
    maxEnd = max(geneCopyNumbers$end)

    vec = vector()
    vec[maxEnd - minStart + 1] <- 0
    vec[] <- 0

    for (i in 1:nrow(geneCopyNumbers)) {
        a = geneCopyNumbers[i,"start"] - minStart + 1
        b = geneCopyNumbers[i,"end"] - minStart + 1

        vec[a:b] <- vec[a:b] + 1
    }

    df = data.frame(x = c(minStart:maxEnd), y = vec)
    result = df %>% mutate(lag_y = lag(y), lead_y = lead(y)) %>% filter(is.na(lag_y) | is.na(lead_y) | lag_y != lead_y) %>%
        mutate(x_end = lead(x)) %>%
        select(x, x_end, y) %>%
        filter(!is.na(x_end))


    return (result)

}



create_plot <- function(geneName, canonicalTranscripts, hpcCopyNumbers) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.5)
    strictGeneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.2)
    
    minStart = min(geneCopyNumbers$start)
    maxEnd = max(geneCopyNumbers$end)

    geneHistogramData = create_histogram_data(geneCopyNumbers)
    strictHistogramData = create_histogram_data(strictGeneCopyNumbers)
    
    geneHistogramMax = max(geneHistogramData$y)

    overlappingGenes = canonicalTranscripts %>% filter(chromosome == primaryGene$chromosome, geneStart <= maxEnd + 20000, geneEnd >= minStart - 20000) %>%
        mutate(value = (row_number() %% 15 + 1) * geneHistogramMax / 15 ) %>%
        mutate(gene = ifelse(is.na(pValue), gene, paste0(gene, "(", pValue, ")")))

    p = ggplot() +
        geom_rect(data=geneHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.6) +
        geom_rect(data=strictHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 1) +
        geom_rect(data=overlappingGenes, aes(xmin = geneStart, xmax = geneEnd, ymin = 0, ymax = value, fill = gene), alpha = 0.5) +
        geom_text(data=overlappingGenes, aes(x = (geneStart + geneEnd)/2, y = value, label = gene), hjust = 0.5, size = 2, nudge_y = 0.5) +
        theme(legend.position="none") +  xlab("Position") + ylab("Deletes") + ggtitle(paste0(geneName, " (range=", approximate_distance(minStart, maxEnd), ")"))

    return (p)
}



load("~/hmf/RData/Reference/hpcCopyNumbers.RData")
load("~/hmf/RData/Reference/canonicalTranscripts.RData")
essential = read.csv(file = "/Users/jon/hmf/analysis/essential/NIHMS732683-supplement-supp_table_3.csv", stringsAsFactors = F) %>%
    select(gene = Gene, pValue = KBM7.adjusted.p.value) %>% mutate(pValue = paste0(round(pValue, 2)))
canonicalTranscripts = canonicalTranscripts %>% left_join(essential, by = "gene")

load("~/hmf/RData/Processed/hpcDriversByGene.RData")
genesToExamine = hpcDriversByGene %>% filter(driver %in% c("Deletion","FragileDel")) %>% select(gene) %>% distinct() %>% filter(!grepl("telomere", gene), !grepl("centromere", gene))  %>% arrange(gene)


plots = list()
for (geneName in genesToExamine$gene) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.5)
    if (nrow(geneCopyNumbers) > 10) {
        cat (geneName, "\n")
        plots[[geneName]] <- create_plot(geneName , canonicalTranscripts, hpcCopyNumbers)
    }
}

pdf(file="/Users/jon/hmf/analysis/essential/EssentialGenes.pdf",width=15, height = 6)
for (i in 1:length(plots)) {
    print(plots[[i]])
}
dev.off()


geneName = "TP53"
create_plot("TP53" , canonicalTranscripts, hpcCopyNumbers)
create_plot("PTEN" , canonicalTranscripts, hpcCopyNumbers)
create_plot("CDKN1B" , canonicalTranscripts, hpcCopyNumbers)
create_plot("RAD51B" , canonicalTranscripts, hpcCopyNumbers)



#######  VIOLIN PLOT OF LENGTH OF DELS
load("~/hmf/RData/Processed/hpcDriversByGene.RData")
load("~/hmf/RData/Reference/hpcCopyNumbers.RData")
load("~/hmf/RData/Reference/canonicalTranscripts.RData")
genesToExamine = hpcDriversByGene %>% 
  filter(driver %in% c("Deletion","FragileDel"), !grepl("telomere", gene), !grepl("centromere", gene)) %>% 
  group_by(gene, driver) %>% 
  count() %>% 
  ungroup() %>%
  top_n(200, n) %>%
  left_join(canonicalTranscripts, by = "gene") %>% arrange(-n)

hpcDels = hpcCopyNumbers %>% filter(copyNumber <= 0.5)
hpcDelRanges = GRanges(hpcDels$chromosome, IRanges(hpcDels$start, hpcDels$end))
geneRanges = GRanges(genesToExamine$chromosome, IRanges(genesToExamine$geneStart, genesToExamine$geneEnd))
ol = as.matrix(findOverlaps(hpcDelRanges, geneRanges, type="any", select="all"))

hpcDelsInGenes = hpcDels[ol[, 1], ]
hpcDelsInGenes$gene <- genesToExamine[ol[,2], "gene"]$gene
hpcDelsInGenes$driver <- genesToExamine[ol[,2], "driver"]$driver

hpcDelsInGenes = hpcDelsInGenes %>% mutate(
  length = end - start, 
  lengthFactor = cut(length, breaks = c(0,1,10,100,1000,10000,100000,1000000,10000000000)),
  geneFactor = factor(gene, levels = genesToExamine$gene, ordered = T)) %>%
  filter(length > 0)

## VALIDATION
#jon = hpcDelsInGenes %>% group_by(sampleId, chromosome, start, end) %>% summarise(n = n(), genes = paste0(gene, collapse = ",") ) %>% filter(n > 1)
#jon2 = hpcDelsInGenes %>% group_by(sampleId, gene) %>% count()  %>% filter(n > 1)
  
filteredGenes = hpcDelsInGenes %>% filter(gene %in% genesToExamine[1:30, ]$gene) 
plotViolin <- function(df) {
  ggplot(df, aes(geneFactor, length)) + 
  geom_violin( aes(fill=driver), draw_quantiles = c(0.25, 0.5, 0.75), scale = "count") +
  xlab("Gene") + ylab("DelLength") +
  scale_y_continuous(trans="log10",labels = comma)
}

p1 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[1:30, ]$gene))
p2 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[31:60, ]$gene))
p3 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[61:90, ]$gene))

cowplot::plot_grid(p1, p2, p3, ncol = 1)

p4 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[91:120, ]$gene))
p5 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[121:150, ]$gene))
p6 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[151:183, ]$gene))

cowplot::plot_grid(p4, p5, p6, ncol = 1)




  