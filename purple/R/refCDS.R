library(GenomicRanges)


library(seqinr)
library(BSgenome.Hsapiens.UCSC.hg19)

refGenome <- BSgenome.Hsapiens.UCSC.hg19
interval_cds = KRAS$intervals_cds
strand = KRAS$strand
str(KRAS$seq_cds)
str(interval_cds)
KRAS
chromosome = "chr12"
KRAS$seq_cds

myKRASCDS = createExonCDS("chr12", KRAS$strand, KRAS$intervals_cds, refGenome)
myKRASCDS$seq_cds1down
KRAS$seq_cds1down

AL606500$chr
myAL606500 = createExonCDS("chr1", AL606500$strand, AL606500$intervals_cds, refGenome)
AL606500$intervals_cds
myAL606500$intervals_splice

KRAS$intervals_cds
KRAS$intervals_splice
intervals_splice

createIntervalsSplice <-function(chromosome, strand, interval_cds, spliceOffsets = c(-2, -1, 1, 2, 5)) {
  correctedSpliceOffsets = strand * spliceOffsets

  positiveOffsets = correctedSpliceOffsets[correctedSpliceOffsets > 0]
  positivePositions = unlist(lapply(interval_cds[,2], function (x) {x + positiveOffsets}))

  negativeOffsets = correctedSpliceOffsets[correctedSpliceOffsets < 0]
  negativePositions = unlist(lapply(interval_cds[,1], function (x) {x + negativeOffsets}))

  allPositions = c(positivePositions, negativePositions)
  intervalsSplice = sort(allPositions[allPositions > min(interval_cds) & allPositions < max(interval_cds)])
}

createSeqCDS <-function(chromosome, strand, interval_cds, refGenome) {
  priorPosition = min(interval_cds) - 1
  postPosition = max(interval_cds) + 1

  sequence = DNAString(paste(apply(interval_cds, 1, function (x) {as.character(getSeq(refGenome, GRanges(chromosome, IRanges(x[1],x[2]))))}), collapse = ""))
  priorSequence = as.character(getSeq(refGenome, GRanges(chromosome, IRanges(priorPosition,priorPosition))))
  postSequence = as.character(getSeq(refGenome, GRanges(chromosome, IRanges(postPosition,postPosition))))

  cds = DNAString(paste(priorSequence, sequence, postSequence, sep = ""))
  if (strand == -1) {
    cds = reverseComplement(cds)
  }

  seq_cds1up = cds[1:(length(cds)-2)]
  seq_cds = cds[2:(length(cds)-1)]
  seq_cds1down = cds[3:length(cds)]

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
  refCodon = as.character(as.vector(seq_cds[codonPos]))

  altCodon = refCodon;
  altCodon[cdsPos-(ceiling(cdsPos/3)-1)*3] = alt

  refProtein = seqinr::translate(refCodon)
  altProtein = seqinr::translate(altCodon)
  impact = impact(refProtein, altProtein)

  return (setNames(c(cdsPos, ref, alt, refTri, altTri, impact), c("cdsPos", "ref", "alt", "refTri", "altTri", "impact")))
}

## OLD
exonImpact <-function(position, ref, alt, geneRefCDS) {
  cdsPos = chr2cds(position, geneRefCDS$intervals_cds, geneRefCDS$strand)

  refTri = paste(geneRefCDS$seq_cds1up[cdsPos], geneRefCDS$seq_cds[cdsPos], geneRefCDS$seq_cds1down[cdsPos], sep = "")
  altTri = paste(geneRefCDS$seq_cds1up[cdsPos], alt, geneRefCDS$seq_cds1down[cdsPos], sep = "")

  codonPos = c(ceiling(cdsPos/3)*3-2, ceiling(cdsPos/3)*3-1, ceiling(cdsPos/3)*3)
  refCodon = as.character(as.vector(geneRefCDS$seq_cds[codonPos]))

  altCodon = refCodon;
  altCodon[cdsPos-(ceiling(cdsPos/3)-1)*3] = alt

  refProtein = seqinr::translate(refCodon)
  altProtein = seqinr::translate(altCodon)
  impact = impact(refProtein, altProtein)

  return (c(position, ref, alt, refTri, altTri, impact))
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









