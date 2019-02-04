library(dplyr)
library(ggplot2)

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

    minStart = min(geneCopyNumbers$start)
    maxEnd = max(geneCopyNumbers$end)

    geneHistogramData = create_histogram_data(geneCopyNumbers)
    geneHistogramMax = max(geneHistogramData$y)

    overlappingGenes = canonicalTranscripts %>% filter(chromosome == primaryGene$chromosome, geneStart <= maxEnd, geneEnd >= minStart) %>%
        mutate(value = (row_number() %% 10 + 1) * geneHistogramMax / 10 ) %>%
        mutate(geneEnd = pmin(geneEnd, maxEnd), gene = ifelse(is.na(pValue), gene, paste0(gene, "-", pValue)))

    p = ggplot() +
        geom_rect(data=geneHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 1) +
        geom_rect(data=overlappingGenes, aes(xmin = geneStart, xmax = geneEnd, ymin = 0, ymax = value, fill = gene), alpha = 0.5) +
        geom_text(data=overlappingGenes, aes(x = (geneStart + geneEnd)/2, y = value, label = gene), hjust = 0.5, size = 2, nudge_y = 0.5) +
        theme(legend.position="none") +  xlab("Position") + ylab("Deletes") + ggtitle(geneName)

    return (p)
}


load("~/hmf/RData/Reference/hpcCopyNumbers.RData")
load("~/hmf/RData/Reference/canonicalTranscripts.RData")
load("~/hmf/RData/Processed/driverGenes.RData")
essential = read.csv(file = "/Users/jon/hmf/analysis/essential/NIHMS732683-supplement-supp_table_3.csv", stringsAsFactors = F) %>%
    select(gene = Gene, pValue = KBM7.adjusted.p.value) %>% mutate(pValue = paste0(round(pValue, 2)))
canonicalTranscripts = canonicalTranscripts %>% left_join(essential, by = "gene")


genesToExamine = hpcCopyNumbers %>% filter(copyNumber <= 0.5) %>% group_by(gene) %>% count()

plots = list()
for (geneName in tsGenes$gene_name) {
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
create_plot("SMAD4" , canonicalTranscripts, hpcCopyNumbers)




