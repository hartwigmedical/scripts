library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(dplyr)
library(stringr)
gridss.short_deldup_size_threshold = 1000
gridss.allowable_normal_contamination=0.03
gridss.use_read_threshold=3
gridss.min_breakpoint_depth = 10
gridss.max_homology_length = 50

#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_filter = function(gr, vcf, somatic_filters=TRUE) {
	vcf = vcf[names(gr)]
	i = info(vcf)
	g = geno(vcf)
	isShort = is_short_deldup(gr)
	ihomlen <- rowSums(abs(as.matrix(info(vcf)$IHOMPOS)))
	ihomlen[is.na(ihomlen)] = 0
	filtered <- str_detect(gr$FILTER, "NO_ASSEMBLY") |
		str_detect(gr$FILTER, "LOW_QUAL") |
		ihomlen > gridss.max_homology_length |
		# BPI.Filter.MinDepth
		(i$SR + i$RP < gridss.min_breakpoint_depth) |
		# BPI.Filter.PRSupportZero
		(!isShort & i$RP == 0) |
		# BPI.Filter.SRSupportZero
		(isShort & i$SR == 0)
	if (somatic_filters) {
		# see https://github.com/PapenfussLab/gridss/issues/114 for why we can't default to using QUAL
		support = g$SR
		rpsupport = g$RP
		rpsupport[as.logical(is_short_deldup(gr))] = 0
		support = support + rpsupport
		has_sufficient_read_support = rowSums(support) >= gridss.use_read_threshold
		filtered <- filtered |
			(has_sufficient_read_support & support[,1] > 0.03 * support[,2]) |
			(!has_sufficient_read_support & g$QUAL[,1] > 0.03 * g$QUAL[,2])
			# Filter.SRNormalSupport is already covered by above logic
			#(isShort & g$SR[,1] != 0)
	}
	return(as.logical(filtered))
}
is_short_deldup = function(gr) {
	strand(gr) != strand(partner(gr)) &
		seqnames(gr) == seqnames(partner(gr)) &
		abs(start(gr)-start(partner(gr))) < gridss.short_deldup_size_threshold
}
gridss_af = function(gr, vcf, i=1) {
	g = geno(vcf[names(gr)])
	ref <- g$REF[,i]
	refpair <- g$REFPAIR[,i]
	# assembly counts must be halved for reciprocal assemblies since
	# we independently eassemble the same reads on both sides
	assembled_read_weight = ifelse(str_detect(gr$FILTER, "SINGLE_ASSEMBLY"), 1, 0.5)
	assrrp_af <- ((g$ASSR[,i] + g$ASRP[,i]) * assembled_read_weight)/((g$ASSR[,i] + g$ASRP[,i]) * assembled_read_weight + ref + refpair)
	srrp_af <- (g$SR[,i] + g$RP[,i]) / (g$SR[,i] + g$RP[,i] + ref + refpair)
	assr_af <- (g$ASSR[,i] * assembled_read_weight) / (g$ASSR[,i] * assembled_read_weight + ref)
	# it is questionably whether we should incoproate BSC.
	# It's definitely incorrect in the case of promiscuous breakpoints
	sr_af <- (g$SR[,i]+g$BSC[,i])/(g$SR[,i]+g$BSC[,i]+ref)
	
	# a complete assembly will give the most accurate AF, but
	# a fragmented assembly will underestimate
	ifelse(is_short_deldup(gr), pmax(assr_af, sr_af), pmax(assrrp_af, srrp_af))
}
gridss_somatic_af = function(gr, vcf) {
	return(gridss_af(gr, vcf, 2))
}
annotate_overlaps = function(query, subject, ...) {
	hits = findBreakpointOverlaps(query, subject, ...)
	hits = hits %>%
		mutate(queryQUAL = query$QUAL[queryHits]) %>%
		mutate(subjectQUAL = subject$QUAL[subjectHits]) %>%
		mutate(QUAL=ifelse(is.na(queryQUAL), subjectQUAL, ifelse(is.na(subjectQUAL), queryQUAL, pmax(queryQUAL, subjectQUAL)))) %>%
		replace_na(list(QUAL=0)) %>%
		arrange(QUAL)
	result = rep(NA_character_, length(query))
	result[hits$queryHits] = names(subject)[hits$subjectHits]
	return(result)
}
gr_join_to_df = function(gr1, vcf1, gr2, vcf2, ..., suffix=c(".1", ".2")) {
	gr1$key = paste0(names(gr1), "_", annotate_overlaps(gr1, gr2, ...))
	gr2$key = paste0(annotate_overlaps(gr2, gr1, ...), "_", names(gr2))
	df1 = as.data.frame(gr1)
	df2 = as.data.frame(gr2)
	df1$orientation = paste0(strand(gr1), strand(partner(gr1)))
	df2$orientation = paste0(strand(gr2), strand(partner(gr2)))
	if (!is.null(vcf1)) {
		df1 %>% bind_cols(as.data.frame(info(vcf1[gr1$vcfId])))
	}
	if (!is.null(vcf2)) {
		df2 %>% bind_cols(as.data.frame(info(vcf2[gr2$vcfId])))
	}
	merged_df = full_join(df1, df2, by="key", suffix=suffix)
	return(merged_df)
}
query_structural_variants_for_sample_as_granges = function(dbConnect, sampleId) {
	sampleIdString = paste("'", sampleId, "'", collapse = ",", sep = "")
	query = paste(
		"SELECT *",
		"FROM structuralVariant",
		"WHERE sampleId in (", sampleIdString, ")",
		sep = " ")
	df = dbGetQuery(dbConnect, query)
	hmf_structuralVariants_to_granges(df)
}

hmf_structuralVariants_to_granges = function(df) {
	gro = GRanges(
		seqnames=df$startChromosome,
		ranges=IRanges(start=df$startPosition, width=1),
		strand=ifelse(df$startOrientation == 1, "+", "-"),
		id=df$id,
		sampleId=df$sampleId,
		ploidy=df$ploidy,
		insertSequence=df$insertSequence,
		type=df$type,
		af=df$startAF,
		homseq=df$startHomologySequence,
		adjustedaf=df$adjustedStartAF,
		adjustedcn=df$adjustedStartCopyNumber,
		adjustedcn_delta=df$adjustedStartCopyNumberChange,
		partner=paste0(df$id, "h"))
	names(gro)=paste0(df$id, "o")
	grh = GRanges(
		seqnames=df$endChromosome,
		ranges=IRanges(start=df$endPosition, width=1),
		strand=ifelse(df$endOrientation == 1, "+", "-"),
		id=df$id,
		sampleId=df$sampleId,
		ploidy=df$ploidy,
		insertSequence=df$insertSequence,
		type=df$type,
		af=df$startAF,
		homseq=df$endHomologySequence,
		adjustedaf=df$adjustedEndAF,
		adjustedcn=df$adjustedEndCopyNumber,
		adjustedcn_delta=df$adjustedEndCopyNumberChange,
		partner=paste0(df$id, "o"))
	names(grh)=paste0(df$id, "h")
	c(gro, grh)
}

load_full_gridss_gr = function(filename, directory="", filter=c("somatic", "qual")) {
	vcf = readVcf(paste0(directory, filename), "hg19")
	if ("somatic" %in% filter) {
		vcf = gridss_somatic_filter(vcf)
	}
	if ("qual" %in% filter) {
		vcf = gridss_qual_filter(vcf)
	}
	gr = breakpointRanges(vcf)
	gr$filename = filename
	gr$af = gridss_somatic_af(gr, vcf)
	info_df = info(vcf[gr$vcfId,])
	info_df$filtered = should_filter(vcf[gr$vcfId,])
	mcols(gr) = bind_cols(as.data.frame(mcols(gr)), as.data.frame(info_df))
	return(gr)
}
