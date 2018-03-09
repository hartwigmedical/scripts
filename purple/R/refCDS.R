rm(list=ls())

library(GenomicRanges)
library(foreach)
library(parallel)
library(doParallel)
library(seqinr)
library(BSgenome.Hsapiens.UCSC.hg19)


createGenomicRangesForRefCDS <- function(RefCDS) {
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores, type="FORK")
  geneRanges = parLapply(cl, RefCDS, function(x) {createGenomicRangesForGene(x)})
  stopCluster(cl)
  result = do.call(c, geneRanges)
  return (result)
}

createGenomicRangesForGene <- function(gene) {
  cds_ranges = apply(gene$intervals_cds, 1 , function(x) {c(GRanges(gene$chr, strand = "*", IRanges(x[1], x[2]), names = gene$gene_name))})
  splice_ranges = lapply(gene$intervals_splice, function(x) {c(GRanges(gene$chr, strand = "*", IRanges(x, x), names = gene$gene_name))})
  combinedRanges = append(cds_ranges, splice_ranges)
  return (do.call(c, combinedRanges))
}

saveEnsemblData <- function() {
  dbEnsembl = dbConnect(MySQL(), dbname='homo_sapiens_core_89_37', groups="RAnalysisWrite")
  ensemblExons = query_exons_from_ensembl(dbEnsembl)
  save(ensemblExons, file = "~/hmf/RData/ensemblExons.RData")
  dbDisconnect(dbEnsembl)
  rm(dbEnsembl)
}

createGeneRef<-function(geneExons, refGenome) {
  firstRecord = geneExons[1,]
  strand = firstRecord$strand
  chromosome = firstRecord$chromosome
  gene = list(gene_name = firstRecord$gene_name, gene_id = firstRecord$gene_id, transcript_id = firstRecord$transcript_id, chr = chromosome, strand = strand)

  # Intervals CDS
  exons = geneExons[, c("exon_start", "exon_end")]
  codingStarts = pmax(exons$exon_start, firstRecord$coding_start)
  codingEnds = pmin(exons$exon_end, firstRecord$coding_end)
  codingExons = data.frame(start = codingStarts, end = codingEnds)
  intervals_cds = as.matrix(codingExons[codingExons$start < codingExons$end, ])
  dimnames(intervals_cds) <- NULL
  gene$intervals_cds = intervals_cds

  # Intervals Splice
  refGenomeChromosome = paste("chr", chromosome, sep ="")
  intervals_splice = createIntervalsSplice(strand, intervals_cds)
  gene$intervals_splice = intervals_splice


  # CDS
  seq_cds = createSeqCDS(refGenomeChromosome, strand, intervals_cds, refGenome)
  gene$CDS_length <- length(seq_cds[["seq_cds"]])
  gene$seq_cds1up <- seq_cds[["seq_cds1up"]]
  gene$seq_cds <- seq_cds[["seq_cds"]]
  gene$seq_cds1down <- seq_cds[["seq_cds1down"]]

  # Splice CDS
  if (length(intervals_splice) > 0) {
    splice_cds = createSpliceCDS(refGenomeChromosome, strand, intervals_splice, refGenome)
    gene$seq_splice <- splice_cds[["seq_splice"]]
    gene$seq_splice1up <- splice_cds[["seq_splice1up"]]
    gene$seq_splice1down <- splice_cds[["seq_splice1down"]]
  }

  return (gene)
}


compareRefCDS <- function(old, new) {
  allIndentical =
    identical(old$gene_name, new$gene_name) &&
    identical(old$CDS_length, new$CDS_length) &&
    identical(old$chr, new$chr)  &&
    identical(old$strand, new$strand) &&
    identical(old$intervals_splice, new$intervals_splice) &&
    identical(as.character(old$seq_cds), as.character(new$seq_cds)) &&
    identical(as.character(old$seq_cds1up), as.character(new$seq_cds1up)) &&
    identical(as.character(old$seq_cds1down), as.character(new$seq_cds1down)) &&
    identical(as.character(old$seq_splice), as.character(new$seq_splice)) &&
    identical(as.character(old$seq_splice1up), as.character(new$seq_splice1up)) &&
    identical(as.character(old$seq_splice1down), as.character(new$seq_splice1down))

  return (allIndentical)
}


createSpliceCDS<-function(chromosome, strand, intervals_splice, refGenome) {

  allSequences = sapply(intervals_splice, function (x) {as.character(getSeq(refGenome, GRanges(chromosome, strand = strand, IRanges(x-1, x+1))))})
  allSequencesMatrix = sapply(strsplit(allSequences, ""), function(x) {x})

  seq_splice1up = DNAString(paste(allSequencesMatrix[1, ], collapse = ""))
  seq_splice = DNAString(paste(allSequencesMatrix[2, ], collapse = ""))
  seq_splice1down = DNAString(paste(allSequencesMatrix[3, ], collapse = ""))

  return (list(seq_splice1up = seq_splice1up, seq_splice = seq_splice, seq_splice1down = seq_splice1down))
}

createIntervalsSplice <-function(strand, intervals_cds, spliceOffsets = c(-2, -1, 1, 2, 5)) {
  correctedSpliceOffsets = strand * spliceOffsets

  positiveOffsets = correctedSpliceOffsets[correctedSpliceOffsets > 0]
  positivePositions = unlist(lapply(intervals_cds[,2], function (x) {x + positiveOffsets}))

  negativeOffsets = correctedSpliceOffsets[correctedSpliceOffsets < 0]
  negativePositions = unlist(lapply(intervals_cds[,1], function (x) {x + negativeOffsets}))

  allPositions = c(positivePositions, negativePositions)
  intervalsSplice = sort(allPositions[allPositions > min(intervals_cds) & allPositions < max(intervals_cds)])
  return (intervalsSplice)
}

createSeqCDS <-function(chromosome, strand, interval_cds, refGenome) {
  priorPosition = min(interval_cds) - 1
  postPosition = max(interval_cds) + 1

  seq_cds = DNAString(paste(apply(interval_cds, 1, function (x) {as.character(getSeq(refGenome, GRanges(chromosome, IRanges(x[1],x[2]))))}), collapse = ""))
  seq_cds1up = DNAString(paste(apply(interval_cds, 1, function (x) {as.character(getSeq(refGenome, GRanges(chromosome, IRanges(x[1]-strand,x[2]-strand))))}), collapse = ""))
  seq_cds1down = DNAString(paste(apply(interval_cds, 1, function (x) {as.character(getSeq(refGenome, GRanges(chromosome, IRanges(x[1]+strand,x[2]+strand))))}), collapse = ""))

  if (strand == -1) {
    seq_cds = reverseComplement(seq_cds)
    seq_cds1up = reverseComplement(seq_cds1up)
    seq_cds1down = reverseComplement(seq_cds1down)
  }

  return (list(seq_cds1up = seq_cds1up, seq_cds = seq_cds, seq_cds1down = seq_cds1down))
}


createRefCDS <- function() {
  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)

  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)

  ### DNDS
  data("refcds_hg19", package="dndscv")
  KRAS = RefCDS[[9224]]
  AL606500=RefCDS[[940]]
  RP11 = RefCDS[[14756]]

  AL606500$strand

  krasL = createLmatrix(KRAS, nt, trinucsubsind)
  AL606500L = createLmatrix(AL606500, nt, trinucsubsind)
  RP11L = createLmatrix(RP11, nt, trinucsubsind)

  identical(krasL,KRAS$L)
  identical(AL606500L,AL606500$L)
  identical(RP11L,RP11$L)

  all(krasL == AL606500$L)

}

createLmatrix <-function(geneRefCDS, nt, trinucsubsind) {
  L = array(0, dim=c(192,4))
  for (cdsPos in (1:length(geneRefCDS$seq_cds))) {
    ref = as.character(geneRefCDS$seq_cds[cdsPos])
    alts = setdiff(nt, ref)
    for (alt in alts) {
      altImpact = exonImpactFromCDSPos(cdsPos, ref, alt, geneRefCDS$seq_cds1up, geneRefCDS$seq_cds, geneRefCDS$seq_cds1down)
      impactInd = impactIndex(altImpact["impact"])
      if (!is.na(impactInd)) {
        triInd = trinucsubsind[ paste(altImpact["refTri"],altImpact["altTri"], sep=">") ]
        L[triInd, impactInd] = L[triInd, impactInd] + 1
      }
    }
  }

  # Splice
  for (cdsPos in (1:length(geneRefCDS$seq_splice))) {
    ref = as.character(geneRefCDS$seq_splice[cdsPos])
    alts = setdiff(nt, ref)
    for (alt in alts) {
      refTri = paste(geneRefCDS$seq_splice1up[cdsPos], geneRefCDS$seq_splice[cdsPos], geneRefCDS$seq_splice1down[cdsPos], sep = "")
      altTri = paste(geneRefCDS$seq_splice1up[cdsPos], alt, geneRefCDS$seq_splice1down[cdsPos], sep = "")
      triInd = trinucsubsind[ paste(refTri,altTri, sep=">") ]
      L[triInd, 4] = L[triInd, 4] + 1
    }
  }

  return (L)
}

exonImpactFromCDSPos <-function(cdsPos, ref, alt, seq_cds1up, seq_cds, seq_cds1down) {

  refTri = paste(seq_cds1up[cdsPos], seq_cds[cdsPos], seq_cds1down[cdsPos], sep = "")
  altTri = paste(seq_cds1up[cdsPos], alt, seq_cds1down[cdsPos], sep = "")

  codonPos = c(ceiling(cdsPos/3)*3-2, ceiling(cdsPos/3)*3-1, ceiling(cdsPos/3)*3)
  if (max(codonPos) <= length(seq_cds)) {
    refCodon = as.character(as.vector(seq_cds[codonPos]))

    altCodon = refCodon;
    altCodon[cdsPos-(ceiling(cdsPos/3)-1)*3] = alt

    refProtein = seqinr::translate(refCodon)
    altProtein = seqinr::translate(altCodon)
    impact = impact(refProtein, altProtein)

    return (setNames(c(cdsPos, ref, alt, refTri, altTri, impact), c("cdsPos", "ref", "alt", "refTri", "altTri", "impact")))
  } else {
    return (setNames(c(cdsPos, ref, alt, refTri, altTri, "Unknown"), c("cdsPos", "ref", "alt", "refTri", "altTri", "impact")))
  }

}

impact <- function(old_aa, new_aa) {
  # Annotating the impact of the mutation
  if (new_aa == old_aa){
    return("Synonymous")
  } else if (new_aa == "*"){
    return("Nonsense")
  } else if (old_aa != "*"){
    return("Missense")
  } else if (old_aa=="*") {
    return ("Stop_loss")
  }
}

impactIndex <- function(impact) {
  if (impact == "Synonymous") {
    return (1)
  } else if (impact == "Nonsense") {
    return (3)
  } else if (impact == "Missense") {
    return (2)
  } else if (impact == "Stop_loss") {
    return (3)
  }

  return (NA)
}

chr2cds = function(pos,cds_int,strand) {
  if (strand==1) {
    return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
  } else if (strand==-1) {
    return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
  }
}

generateHmfRefCDS <- function(ensemblExons, refGenome) {
  codingGenes = sort(unique(ensemblExons[!is.na(ensemblExons$coding_start), c("gene_name") ]))
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores, type="FORK")
  registerDoParallel(cl)
  date()
  HmfRefCDS = foreach(gene_name = codingGenes,
                      .combine = list,
                      .maxcombine = 50000,
                      .multicombine = TRUE)  %dopar%
                      {geneExons = ensemblExons[ensemblExons$gene_name == gene_name, ]; return (createGeneRef(geneExons, refGenome))}
  date()
  stopCluster(cl)
  return(HmfRefCDS)
}

generateHmfRefCDSL <-function(HmfRefCDS) {

  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)
  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)

  HmfRefCDSNames = sapply(HmfRefCDS, function (x) {x$gene_name})
  HmfRefCDSInd = setNames(1:length(HmfRefCDSNames), HmfRefCDSNames)

  date()
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores, type="FORK")
  HmfRefCDSL = parLapply(cl, HmfRefCDS, function (x) {(list(gene_name=x$gene_name, L = createLmatrix(HmfRefCDS[[HmfRefCDSInd[x$gene_name]]], nt, trinucsubsind)))})
  stopCluster(cl)
  date()

  return (HmfRefCDSL)
}

attachLToRefCDS <- function(HmfRefCDS, HmfRefCDSL) {

  HmfRefCDSNames = sapply(HmfRefCDS, function (x) {x$gene_name})
  HmfRefCDSInd = setNames(1:length(HmfRefCDSNames), HmfRefCDSNames)

  HmfRefCDSLNames = sapply(HmfRefCDSL, function (x) {x$gene_name})
  HmfRefCDSLInd = setNames(1:length(HmfRefCDSLNames), HmfRefCDSLNames)

  for (gene_name in HmfRefCDSLNames) {
    HmfRefCDS[[HmfRefCDSInd[gene_name]]]$L <- HmfRefCDSL[[HmfRefCDSLInd[gene_name]]]$L
  }

  return (HmfRefCDS)
}


createCleanRefCDS <-function(HmfRefCDS, HmfRefCDSL, HmfRefCDSL2) {

  data("covariates_hg19", package="dndscv")
  covNames = rownames(covs)

  result = attachLToRefCDS(HmfRefCDS, HmfRefCDSL)
  result = attachLToRefCDS(result, HmfRefCDSL2)

  validCodons = sapply(result, function (x) {x$CDS_length %% 3 == 0})
  result = result[validCodons]

  validCovNames = sapply(result, function (x) {x$gene_name %in% covNames})
  result = result[validCovNames]

  return (result)
}