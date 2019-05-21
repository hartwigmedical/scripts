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
load("~/hmf/analysis/cohort/canonicalTranscripts.RData")
essential = read.csv(file = "/Users/jon/hmf/analysis/essential/NIHMS732683-supplement-supp_table_3.csv", stringsAsFactors = F) %>%
  select(gene = Gene, pValue = KBM7.adjusted.p.value) %>% mutate(pValue = paste0(round(pValue, 2)))
canonicalTranscripts = canonicalTranscripts %>% left_join(essential, by = "gene")

load("~/hmf/RData/Processed/hpcDriversByGene.RData")
genesToExamine = hpcDriversByGene %>% filter(driver %in% c("Amplification"), !grepl("telomere", gene), !grepl("centromere", gene)) %>%
  group_by(gene) %>% count()  %>% arrange(gene)

hpcCopyNumbers = hpcCopyNumbers %>% left_join(allPurity[, c("sampleId", "ploidy")], by = "sampleId")

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



load("~/hmf/RData/Reference/allPurity.RData")
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

load("~/hmf/RData/Reference/allPurity.RData")
load("~/hmf/RData/Reference/hpcCopyNumbers.RData")
load("~/hmf/RData/Reference/canonicalTranscripts.RData")
load(file = "/Users/jon/hmf/analysis/essential/transcripts.RData")

hpcCopyNumbers = hpcCopyNumbers %>% left_join(allPurity[, c("sampleId", "cancerType")], by = "sampleId")

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

create_plot("FHIT", canonicalTranscripts, hpcCopyNumbers, transcripts)

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




  