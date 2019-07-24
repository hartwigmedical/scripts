library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(scales)
library(RMySQL)

approximate_distance <- function(start, end) {
    difference = end - start;
    if (difference > 1000000) {
        return (paste0(round(difference/1000000), "Mb"))
    }

    return (paste0(round(difference/1000), "kb"))
}


create_histogram_data <- function (geneCopyNumbers) {
    if (nrow(geneCopyNumbers) == 0) {
        return (data.frame())
    }

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

########## AMPLIFICATIONS ########## 
load("~/hmf/analysis/cohort/cohort.RData")
load("~/hmf/analysis/cohort/hpcCopyNumbers.RData")
load("~/hmf/analysis/cohort/reference/canonicalTranscripts.RData")
essential = read.csv(file = "/Users/jon/hmf/analysis/essential/NIHMS732683-supplement-supp_table_3.csv", stringsAsFactors = F) %>%
    select(gene = Gene, pValue = KBM7.adjusted.p.value) %>% mutate(pValue = paste0(round(pValue, 2)))
canonicalTranscripts = canonicalTranscripts %>% left_join(essential, by = "gene")

load("~/hmf/RData/Processed/hpcDriversByGene.RData")
genesToExamine = hpcDriversByGene %>% filter(driver %in% c("Amplification"), !grepl("telomere", gene), !grepl("centromere", gene)) %>%
    group_by(gene) %>% count()  %>% arrange(gene)

hpcCopyNumbers = hpcCopyNumbers %>% left_join(highestPurityCohort[, c("sampleId", "ploidy")], by = "sampleId")

create_plot <- function(geneName, canonicalTranscripts, hpcCopyNumbers) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber >= 1 * ploidy)

    geneHistogramData = create_histogram_data(geneCopyNumbers %>% filter(copyNumber >= 3 * ploidy))
    geneHistogramData5 = create_histogram_data(geneCopyNumbers %>% filter(copyNumber >= 5 * ploidy))
    geneHistogramData10 = create_histogram_data(geneCopyNumbers %>% filter(copyNumber >= 10 * ploidy))
    geneHistogramData20 = create_histogram_data(geneCopyNumbers %>% filter(copyNumber >= 20 * ploidy))

    minStart = min((geneCopyNumbers %>% filter(copyNumber >= 5 * ploidy))$start)
    maxEnd = max((geneCopyNumbers %>% filter(copyNumber >= 5 * ploidy))$end)

    geneHistogramMax = max(geneHistogramData$y)
    geneHistogram5Max = max(geneHistogramData5$y)

    lowerSide = geneHistogramData5 %>% filter(y < 0.05 * geneHistogram5Max, x < primaryGene$geneStart) %>% select(x)
    if (nrow(lowerSide) == 0) {
        xMin = minStart
    } else {
        xMin = pmin(primaryGene$geneStart - 10000, max(lowerSide))
    }

    upperSide = geneHistogramData5 %>% filter(y < 0.05 * geneHistogram5Max, x_end > primaryGene$geneEnd) %>% select(x)
    if (nrow(upperSide) == 0) {
        xMax = maxEnd
    } else {
        xMax = pmax(primaryGene$geneStart - 10000, min(upperSide))
    }

    overlappingGenes = canonicalTranscripts %>% filter(chromosome == primaryGene$chromosome, geneStart <= maxEnd + 20000, geneEnd >= minStart - 20000) %>%
        mutate(value = (row_number() %% 15 + 1) * geneHistogramMax / 15 ) %>%
        mutate(gene = ifelse(is.na(pValue), gene, paste0(gene, "(", pValue, ")")))

    p = ggplot() +
        geom_rect(data=geneHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.3) +
        geom_rect(data=geneHistogramData5, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.5)

    if (nrow(geneHistogramData10) != 0) {
        p = p + geom_rect(data=geneHistogramData10, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.7)
    }

    if (nrow(geneHistogramData20) != 0) {
        p = p + geom_rect(data=geneHistogramData20, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 1)
    }

    p = p +
        geom_rect(data=overlappingGenes, aes(xmin = geneStart, xmax = geneEnd, ymin = 0, ymax = value, fill = gene), alpha = 0.5) +
        geom_text(data=overlappingGenes, aes(x = (geneStart + geneEnd)/2, y = value, label = gene), hjust = 0.5, size = 2, nudge_y = 0.5) +
        theme(legend.position="none") +  xlab("Position") + ylab("Amps") + ggtitle(paste0(geneName, " (range=", approximate_distance(minStart, maxEnd), ")")) +
        coord_cartesian(xlim = c(xMin, xMax))

    return (p)
}

geneName = "BRAF"
create_plot("BRAF", canonicalTranscripts, hpcCopyNumbers)

plots = list()
for (geneName in genesToExamine$gene) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber >= 5 * ploidy)
    if (nrow(geneCopyNumbers) > 10) {
        cat (geneName, "\n")
        plots[[geneName]] <- create_plot(geneName , canonicalTranscripts, hpcCopyNumbers)
    }
}


pdf(file="/Users/jon/hmf/analysis/essential/EssentialAmps.3.5.10.20.pdf",width=15, height = 6)
for (i in 1:length(plots)) {
    print(plots[[i]])
}
dev.off()


########## DELETES ########## 

maxEnd = 61467686

create_plot <- function(geneName, canonicalTranscripts, hpcCopyNumbers, maxEnd = NA) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.5)
    strictGeneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.2)
    geneCopyNumberDeleteCount = geneCopyNumbers %>% group_by(sampleId) %>% count()
    nrow(geneCopyNumberDeleteCount)

    minStart = min(geneCopyNumbers$start)
    maxEnd = ifelse(is.na(maxEnd), max(geneCopyNumbers$end), maxEnd)

    strictHistogramData = create_histogram_data(strictGeneCopyNumbers)
    geneHistogramData = create_histogram_data(geneCopyNumbers)
    geneHistogramMax = max(geneHistogramData$y)

    overlappingGenes = canonicalTranscripts %>% filter(chromosome == primaryGene$chromosome, geneStart <= maxEnd + 20000, geneEnd >= minStart - 20000) %>%
        mutate(value = (row_number() %% 15 + 1) * geneHistogramMax / 15 ) %>%
        mutate(gene = ifelse(is.na(pValue), gene, paste0(gene, "(", pValue, ")")))

    p = ggplot() +
        geom_rect(data=geneHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.6) +
        geom_rect(data=strictHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 1) +
        geom_rect(data=overlappingGenes, aes(xmin = geneStart, xmax = geneEnd, ymin = 0, ymax = value, fill = gene), alpha = 0.5) +
        geom_text(data=overlappingGenes, aes(x = (geneStart + geneEnd)/2, y = value, label = gene), hjust = 0.5, size = 2, nudge_y = 0.5) +
        theme(legend.position="none") +  xlab("Position") + ylab("Deletes") + ggtitle(paste0(geneName, " (range=", approximate_distance(minStart, maxEnd), " n=", nrow(geneCopyNumberDeleteCount), ")"))

    return (p)
}

geneName = "TP53"
maxEnd = NA
create_plot(geneName, canonicalTranscripts, hpcCopyNumbers)

pFIT = create_plot("FHIT", canonicalTranscripts %>% filter(gene == "FHIT"), hpcCopyNumbers, 61467686)  +  coord_cartesian(xlim = c(59425278, 61403445))
pMacrod = create_plot("MACROD2", canonicalTranscripts%>% filter(grepl("MACROD", gene)), hpcCopyNumbers)
pwwox = create_plot("WWOX", canonicalTranscripts %>% filter(gene == "WWOX"), hpcCopyNumbers) +  coord_cartesian(xlim = c(78134148, 79285923))

pDelSpots = plot_grid(pFIT, pwwox, pMacrod, nrow = 1)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DelSpots.pdf", pDelSpots, width = 189, height = 80, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DelSpots.png", pDelSpots, width = 189, height = 80, units = "mm", dpi = 300)


#load("~/hmf/RData/Reference/allPurity.RData")
#load("~/hmf/RData/Reference/hpcCopyNumbers.RData")
#load("~/hmf/RData/Reference/canonicalTranscripts.RData")
#essential = read.csv(file = "/Users/jon/hmf/analysis/essential/NIHMS732683-supplement-supp_table_3.csv", stringsAsFactors = F) %>%
#    select(gene = Gene, pValue = KBM7.adjusted.p.value) %>% mutate(pValue = paste0(round(pValue, 2)))
#canonicalTranscripts = canonicalTranscripts %>% left_join(essential, by = "gene")

load("~/hmf/RData/Processed/hpcDriversByGene.RData")
#genesToExamine = hpcDriversByGene %>% filter(driver %in% c("Deletion","FragileDel")) %>% select(gene) %>% distinct() %>% filter(!grepl("telomere", gene), !grepl("centromere", gene))  %>% arrange(gene)
genesToExamine = hpcDelsInGenes %>% select(gene, sampleId) %>% distinct() %>% group_by(gene) %>% count() %>% filter(n >= 20) %>% arrange(-n)



plots = list()
for (geneName in genesToExamine$gene) {
    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.5)
    if (nrow(geneCopyNumbers) > 10) {
        cat (geneName, "\n")
        plots[[geneName]] <- create_plot(geneName , canonicalTranscripts, hpcCopyNumbers)
    }
}

pdf(file="~/hmf/analysis/cohort/EssentialDels.pdf",width=15, height = 6)
for (i in 1:length(plots)) {
    print(plots[[i]])
}
dev.off()


geneName = "FHIT"
create_plot("FHIT" , canonicalTranscripts, hpcCopyNumbers)
create_plot("FHIT" , canonicalTranscripts, hpcCopyNumbers, 61467686)
create_plot("TP53" , canonicalTranscripts, hpcCopyNumbers)
create_plot("PTEN" , canonicalTranscripts, hpcCopyNumbers)
create_plot("CDKN1B" , canonicalTranscripts, hpcCopyNumbers)
create_plot("RAD51B" , canonicalTranscripts, hpcCopyNumbers)


####### HISTOGRAM WITH TRANSCRIPTS

query_transcripts <- function(dbConnect, genes) {
    geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
    query = paste0(
    "select gx.display_label as gene, tx.display_label, t.* from gene g, xref gx, transcript t, xref tx ",
    "where g.display_xref_id = gx.xref_id and gx.display_label in (", geneString,") and t.gene_id = g.gene_id ",
    "and t.biotype = 'protein_coding' and tx.xref_id = t.display_xref_id")
    return (dbGetQuery(dbConnect, query))
}

#load("~/hmf/RData/Processed/hpcDriversByGene.RData")
#genes = hpcDriversByGene %>% filter(driver == "FragileDel") %>% select(gene) %>% distinct()
#dbProd = dbConnect(MySQL(), dbname='homo_sapiens_core_89_37', groups="RAnalysis")
#transcripts = query_transcripts(dbProd, genes)
#dbDisconnect(dbProd)
#rm(dbProd)
#save(transcripts, file = "/Users/jon/hmf/analysis/essential/transcripts.RData")

load("~/hmf/analysis/cohort/cohort.RData")
load("~/hmf/analysis/cohort/hpcCopyNumbers.RData")
load("~/hmf/analysis/cohort/canonicalTranscripts.RData")
load(file = "/Users/jon/hmf/analysis/essential/transcripts.RData")

hpcCopyNumbers = hpcCopyNumbers %>% left_join(cohort[, c("sampleId", "cancerType")], by = "sampleId")

create_plot <- function(geneName, canonicalTranscripts, hpcCopyNumbers, transcripts) {

    primaryGene = canonicalTranscripts %>% filter(gene == geneName)
    geneCopyNumbers = hpcCopyNumbers %>% filter(chromosome == primaryGene$chromosome, start <= primaryGene$geneEnd, end >= primaryGene$geneStart, copyNumber <= 0.5)
    transcripts = transcripts %>% filter(gene == geneName)

    transcriptStart = min(transcripts$seq_region_start)
    transcriptEnd = max(transcripts$seq_region_end)
    minStart = transcriptStart - 0.1 * (transcriptEnd - transcriptStart)
    maxEnd = transcriptEnd + 0.1 * (transcriptEnd - transcriptStart)

    cancerTypes = geneCopyNumbers %>% select(cancerType) %>% distinct()
    geneHistogramData = data.frame()
    for (singleCancerType in cancerTypes$cancerType) {
        cancerTypeResult = create_histogram_data(geneCopyNumbers %>% filter(cancerType == singleCancerType))
        cancerTypeResult$cancerType = singleCancerType
        geneHistogramData = bind_rows(cancerTypeResult, geneHistogramData)
    }

    cancersToInclude = geneHistogramData %>% group_by(cancerType) %>% summarise(n = max(y)) %>% filter(n > 4)
    geneHistogramData = geneHistogramData %>% filter(cancerType %in% cancersToInclude$cancerType)

    geneHistogramMax = max(geneHistogramData$y)

    overlappingTranscripts = transcripts %>%
        mutate(length =  seq_region_end - seq_region_start, n = n()) %>%
        arrange(seq_region_start, -length) %>%
        mutate(value = (row_number() %% (n + 1)) * geneHistogramMax / n, ymin = ((row_number() - 1) %% (n + 1)) * geneHistogramMax / n)

    p = ggplot() +
        geom_rect(data=geneHistogramData, aes(xmin = x, xmax = x_end, ymin = 0, ymax = y), alpha = 0.6) +
        geom_rect(data=overlappingTranscripts, aes(xmin = seq_region_start, xmax = seq_region_end, ymin = ymin + 0.1, ymax = value - 0.1, fill = gene), alpha = 0.5) +
        geom_text(data=overlappingTranscripts, aes(x = (seq_region_start + seq_region_end)/2, y = value, label = display_label), vjust = 1, hjust = 0.5, size = 2) +
        theme(legend.position="none") +  xlab("Position") + ylab("Deletes") + ggtitle(paste0(geneName, " (range=", approximate_distance(transcriptStart, transcriptEnd), ")")) +
        facet_wrap(~cancerType) + coord_cartesian(xlim = c(minStart, maxEnd))
    p


    return (p)
}

create_plot("DMD", canonicalTranscripts, hpcCopyNumbers, transcripts)

plots = list()
for (geneName in unique(transcripts$gene)) {
    cat (geneName, "\n")
    plots[[geneName]] <- create_plot(geneName , canonicalTranscripts, hpcCopyNumbers, transcripts)
}

pdf(file="/Users/jon/hmf/analysis/essential/FragileGenesWithTranscripts.pdf",width=10, height = 10)
for (i in 1:length(plots)) {
    print(plots[[i]])
}
dev.off()


load("~/hmf/analysis/cohort/reference/cohort.RData")
sampleIdString = paste("'", highestPurityCohort$sampleId, "'", collapse = ",", sep = "")

library(RMySQL)
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
geneCopyNumberDels = dbGetQuery(dbProd, paste0("select gc.* from geneCopyNumber gc where gc.gene in (select gene from genePanel where reportableDel) and minCopyNumber <= 0.5 and sampleId in (", sampleIdString, ")"))

driverCatalog = dbGetQuery(dbProd, paste0("SELECT * from driverCatalog where sampleId in (", sampleIdString, ")"))
dbDisconnect(dbProd)
rm(dbProd)
save(driverCatalog, file = "~/hmf/analysis/cohort/reference/driverCatalog.RData")
save(geneCopyNumberDels, file = "~/hmf/analysis/cohort/reference/geneCopyNumberDels.RData")


#######  VIOLIN PLOT OF LENGTH OF DELS
load(file = "~/hmf/analysis/cohort/reference/geneCopyNumberDels.RData")
load(file = "~/hmf/analysis/cohort/reference/driverCatalog.RData")
load(file = "~/hmf/RData/Processed/fragileGenes.RData")
load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")
load("~/hmf/analysis/cohort/hpcCopyNumbers.RData")
load("~/hmf/analysis/cohort/reference/canonicalTranscripts.RData")

reportableDels = genePanel %>% filter(reportableDel) %>% select(gene, reportablePointMutation) %>%
    left_join(canonicalTranscripts %>% select(gene, chromosome, geneStart, geneEnd, codingBases), by = "gene") %>%
    left_join(fragileGenes, by = c("gene" = "gene_name")) %>%
    mutate(fragile = ifelse(is.na(fragile), F, fragile),
    status = ifelse(is.na(reportablePointMutation), "HmfDel", "PointMutation"),
    status = ifelse(fragile, "Fragile", status))

pointMutations = driverCatalog %>%
    filter(gene %in% reportableDels$gene, driver == 'MUTATION') %>%
    mutate(driver = "mutationDriverLikelihood") %>%
    group_by(gene, driver) %>%
    summarise(driverLikelihood = sum(driverLikelihood)) %>%
    spread(driver, driverLikelihood, fill = 0)


hpcDels = hpcCopyNumbers %>% filter(copyNumber <= 0.5)
hpcDelRanges = GRanges(hpcDels$chromosome, IRanges(hpcDels$start, hpcDels$end))
geneRanges = GRanges(reportableDels$chromosome, IRanges(reportableDels$geneStart, reportableDels$geneEnd))
ol = as.matrix(findOverlaps(hpcDelRanges, geneRanges, type="any", select="all"))

hpcDelsInGenes = hpcDels[ol[, 1], ]
hpcDelsInGenes$gene <- reportableDels[ol[,2], ]$gene
hpcDelsInGenes$fragile <- reportableDels[ol[,2],]$fragile
hpcDelsInGenes$hmfDeletion <- reportableDels[ol[,2],]$hmfDeletion
hpcDelsInGenes$status <- reportableDels[ol[,2],]$status
hpcDelsInGenes = hpcDelsInGenes %>% mutate(length = end - start + 1) %>% filter(length > 0)

geneCopyNumberDels = geneCopyNumberDels %>% mutate(length = minRegionEnd - minRegionStart + 1)

del_lengths <- function(hpcDelsInGenes) {
    result = data.frame()
    for (selectedGene in unique(hpcDelsInGenes$gene)) {
        geneDels = hpcDelsInGenes %>% filter(gene == selectedGene)
        qValues = data.frame(gene = selectedGene, quantile =  round(t(quantile(geneDels$length, type = 8, name = F))))
        colnames(qValues) <- c("gene", "len_0", "len_25", "len_50", "len_75", "len_100")
        result = bind_rows(qValues, result)
    }

    return (result)
}

del_counts <- function(hpcDelsInGenes) {
    hpcDelsInGenes %>% select(sampleId, gene) %>% distinct() %>% group_by(gene) %>% summarise(delDriverLikelihood = n())
}


lengths = del_lengths(hpcDelsInGenes)
counts = del_counts(hpcDelsInGenes)
delSummary =reportableDels %>%
    left_join(lengths, by = "gene") %>%
    left_join(counts, by = "gene") %>%
    left_join(pointMutations, by = "gene")

delSummary[is.na(delSummary)] <- 0
save(delSummary, file = "~/hmf/analysis/cohort/processed/delSummary.RData")

View(delSummary %>% filter(reportablePointMutation!=0,delDriverLikelihood+mutationDriverLikelihood>20) %>% mutate(geneLength=geneEnd-geneStart,delProp=round(delDriverLikelihood/(delDriverLikelihood+mutationDriverLikelihood),2)))


hpcDelsInGenesCount = hpcDelsInGenes %>% group_by(gene, status) %>% count() %>% arrange(-n)
hpcDelsInGenes = hpcDelsInGenes %>% mutate(geneFactor = factor(gene, levels = hpcDelsInGenesCount$gene, ordered = T))

str(hpcDelsInGenes)
## VALIDATION
#jon = hpcDelsInGenes %>% group_by(sampleId, chromosome, start, end) %>% summarise(n = n(), genes = paste0(gene, collapse = ",") ) %>% filter(n > 1)
#jon2 = hpcDelsInGenes %>% group_by(sampleId, gene) %>% count()  %>% filter(n > 1)

filteredGenes = hpcDelsInGenes %>% filter(gene %in% hpcDelsInGenesCount[1:30, ]$gene)
plotViolin <- function(df) {
    singleBlue = "#6baed6"

    geneDelSummary = delSummary %>% filter(gene %in% df$gene) %>% mutate(geneFactor = factor(gene, levels = unique(df$gene))) %>%
        mutate(delProportion = 7.0 * delDriverLikelihood / (delDriverLikelihood + mutationDriverLikelihood))
    geneLengths = canonicalTranscripts %>% filter(gene %in% df$gene) %>% mutate(length = geneEnd - geneStart + 1, geneFactor = factor(gene, levels = unique(df$gene)))

    ggplot() +
        geom_bar(data = geneDelSummary, mapping = aes(geneFactor, delProportion), stat = "identity") +
        geom_violin(data = df, mapping = aes(geneFactor, length, fill=status), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
        geom_point(data = geneLengths, mapping = aes(geneFactor, length), shape = "cross", size = 5) +
        geom_point(data = geneLengths, mapping = aes(geneFactor, codingBases), size = 5, shape = 1) +
        xlab("Gene") + ylab("DelLength") +
        scale_y_continuous(trans="log10",labels = comma, sec.axis = sec_axis(~.^10))
    #scale_y_continuous(trans="log10",labels = comma, sec.axis = sec_axis(~. / 1, breaks = c(10, 1000, 100000, 10000000), labels = c("25%", "50%", "75%", "100%")))

    p1 = ggplot() +
        geom_bar(data = geneDelSummary, mapping = aes(geneFactor, delProportion), stat = "identity", fill = singleBlue) +
        geom_violin(data = df, mapping = aes(geneFactor, log10(length), fill=status), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
        geom_point(data = geneLengths, mapping = aes(geneFactor, log10(length)), shape = "cross", size = 5) +
        geom_point(data = geneLengths, mapping = aes(geneFactor, log10(codingBases)), size = 5, shape = 1) +
        xlab("Gene") + ylab("DelLength") +
        scale_y_continuous(breaks = c(1,3, 5, 7), labels = c(1, "1,000", "10,000", "10,000,000"), sec.axis = sec_axis(~. / 7, breaks = c(0.25, 0.50, 0.75, 1), labels = c("25%", "50%", "75%", "100%")))

    return (p1)

}

fragileGenes = hpcDelsInGenesCount %>% filter(status == "Fragile")
hmfDelGenes = hpcDelsInGenesCount %>% filter(status == "HmfDel")
pointMutationGenes = hpcDelsInGenesCount %>% filter(status == "PointMutation")

p1 = plotViolin(hpcDelsInGenes %>% filter(gene %in% fragileGenes[1:15, ]$gene))
p2 = plotViolin(hpcDelsInGenes %>% filter(gene %in% hmfDelGenes[1:15, ]$gene))
p3 = plotViolin(hpcDelsInGenes %>% filter(gene %in% pointMutationGenes[1:15, ]$gene))

cowplot::plot_grid(p1, p2, p3, ncol = 1)

p4 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[91:120, ]$gene))
p5 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[121:150, ]$gene))
p6 = plotViolin(hpcDelsInGenes %>% filter(gene %in% genesToExamine[151:183, ]$gene))

cowplot::plot_grid(p4, p5, p6, ncol = 1)




  