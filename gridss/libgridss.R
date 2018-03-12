library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(dplyr)
library(stringr)
library(testthat)
gridss.short_deldup_size_threshold = 1000
gridss.allowable_normal_contamination=0.03
gridss.use_read_threshold=3
gridss.min_breakpoint_depth = 10 # BPI default. This could be significantly lower
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
	ihomlen = rowSums(abs(as.matrix(info(vcf)$IHOMPOS)))
	ihomlen[is.na(ihomlen)] = 0
	filtered = str_detect(gr$FILTER, "NO_ASSEMBLY") |
		str_detect(gr$FILTER, "LOW_QUAL") | # exactly the same as QUAL >= 500
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
		filtered = filtered |
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
	ref = g$REF[,i]
	refpair = g$REFPAIR[,i]
	# assembly counts must be halved for reciprocal assemblies since
	# we independently eassemble the same reads on both sides
	assembled_read_weight = ifelse(str_detect(gr$FILTER, "SINGLE_ASSEMBLY"), 1, 0.5)
	assrrp_af = ((g$ASSR[,i] + g$ASRP[,i]) * assembled_read_weight)/((g$ASSR[,i] + g$ASRP[,i]) * assembled_read_weight + ref + refpair)
	srrp_af = (g$SR[,i] + g$RP[,i]) / (g$SR[,i] + g$RP[,i] + ref + refpair)
	assr_af = (g$ASSR[,i] * assembled_read_weight) / (g$ASSR[,i] * assembled_read_weight + ref)
	# it is questionably whether we should incoproate BSC.
	# It's definitely incorrect in the case of promiscuous breakpoints
	sr_af = (g$SR[,i]+g$BSC[,i])/(g$SR[,i]+g$BSC[,i]+ref)

	# a complete assembly will give the most accurate AF, but
	# a fragmented assembly will underestimate
	ifelse(is_short_deldup(gr), pmax(assr_af, sr_af), pmax(assrrp_af, srrp_af))
}
gridss_somatic_af = function(gr, vcf) {
	return(gridss_af(gr, vcf, 2))
}
annotate_overlaps = function(query, subject, ..., group_by_col_name="sampleId") {
	hits = findBreakpointOverlaps(query, subject, ...)
	if (group_by_col_name %in% names(mcols(query)) & group_by_col_name %in% names(mcols(subject))) {
		hits = hits %>% filter(mcols(query)[[group_by_col_name]][queryHits] == mcols(subject)[[group_by_col_name]][subjectHits])
	}
	hits = hits %>%
		mutate(queryQUAL = query$QUAL[queryHits]) %>%
		mutate(subjectQUAL = subject$QUAL[subjectHits]) %>%
		mutate(QUAL=ifelse(is.na(queryQUAL), subjectQUAL, ifelse(is.na(subjectQUAL), queryQUAL, pmax(queryQUAL, subjectQUAL)))) %>%
		replace_na(list(QUAL=-1)) %>%
		arrange(desc(QUAL)) %>%
		# only allow 1-1 mappings
		filter(!duplicated(queryHits)) %>%
		filter(!duplicated(subjectHits))
	result = rep(NA_character_, length(query))
	result[hits$queryHits] = names(subject)[hits$subjectHits]
	return(result)
}
gr_join_to_df = function(gr1, vcf1, gr2, vcf2, ..., suffix=c(".1", ".2"), group_by_col_name="sampleId") {
	gr1$key = paste0(names(gr1), "_", annotate_overlaps(gr1, gr2, ..., group_by_col_name=group_by_col_name))
	gr2$key = paste0(annotate_overlaps(gr2, gr1, ..., group_by_col_name=group_by_col_name), "_", names(gr2))
	if (group_by_col_name %in% names(mcols(gr1)) & group_by_col_name %in% names(mcols(gr2))) {
		gr1$key = paste(mcols(gr1)[[group_by_col_name]], gr1$key)
		gr2$key = paste(mcols(gr2)[[group_by_col_name]], gr2$key)
	}
	df1 = as.data.frame(gr1)
	df2 = as.data.frame(gr2)
	df1$orientation = paste0(strand(gr1), strand(partner(gr1)))
	df2$orientation = paste0(strand(gr2), strand(partner(gr2)))
	if (!is.null(vcf1)) {
		df1 <- df1 %>% bind_cols(as.data.frame(info(vcf1[gr1$vcfId])))
	}
	if (!is.null(vcf2)) {
		df2 <- df2 %>% bind_cols(as.data.frame(info(vcf2[gr2$vcfId])))
	}
	merged_df = full_join(df1, df2, by="key", suffix=suffix)
	return(merged_df)
}
query_structural_variants_samples = function(dbConnect) {
	query = paste(
		"SELECT DISTINCT sampleId",
		"FROM structuralVariant",
		sep = " ")
	df = dbGetQuery(dbConnect, query)
	return(df$sampleId)
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
		partner=paste0(df$id, "h"),
		svlen=df$endPosition-df$startPosition)
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
		partner=paste0(df$id, "o"),
		svlen=df$endPosition-df$startPosition)
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
transitive_of = function(gr, max_traversal_length, max_positional_error, ...) {
	paths = transitive_paths(gr, max_traversal_length, max_positional_error, ...)
	result = rep(NA_character_, length(gr))
	result[paths$ordinal] = paths$path
	return(result)
}
#' Calculates all transitive paths for the given set of breakpoints
#' @param gr breakpoint GRanges
#' @param max_traversal_length max length of sequence to traverse
#' @param max_positional_error maximum position error in breakpoint calls
#' @param max_depth maximum number of breakpoints to traverse.
#' Highly connected breakpoint sets are O(e^max_depth)
#' @param group_by_col_name column name to separate GRanges by.
#' The default of sampleId causes each sample to calculated independently
transitive_paths = function(gr, max_traversal_length, max_positional_error, max_depth=4, group_by_col_name="sampleId") {
	gr$ordinal = seq(length(gr))
	hits = as.data.frame(findOverlaps(gr, gr, maxgap=max_traversal_length, ignore.strand=TRUE))
	hits = hits %>% filter(hits$subjectHits != hits$queryHits)
	if (group_by_col_name %in% names(mcols(gr))) {
		hits = hits %>% filter(mcols(gr)[[group_by_col_name]][queryHits] == mcols(gr)[[group_by_col_name]][subjectHits])
	}
	hits = hits %>% mutate(
			query_strand = as.character(strand(gr)[queryHits]),
			subject_strand = as.character(strand(gr)[subjectHits])) %>%
		mutate(
			upstream_distance=(start(gr)[queryHits] - start(gr)[subjectHits]) * ifelse(query_strand =="+", 1, -1))

	terminal_hits = hits %>%
		filter(query_strand == subject_strand) %>%
		filter(abs(upstream_distance) < max_positional_error) %>%
		mutate(current=queryHits, terminal_ordinal=subjectHits, terminal_distance=upstream_distance) %>%
		dplyr::select(current, terminal_ordinal, terminal_distance)
	transitive_hits = hits %>%
		filter(query_strand != subject_strand) %>%
		filter(upstream_distance > -max_positional_error) %>%
		mutate(
			transitive_start=queryHits,
			transitive_end=partner(gr)$ordinal[subjectHits],
			transitive_distance=upstream_distance,
			path_name=names(gr)[subjectHits]) %>%
		dplyr::select(transitive_start, transitive_end, transitive_distance, path_name)

	activedf = data.frame(
			ordinal=gr$ordinal,
			terminal_ordinal=partner(gr)$ordinal) %>%
		dplyr::inner_join(terminal_hits, by=c("ordinal"="current"), suffix=c("", ".first")) %>%
		mutate(
			distance=-terminal_distance,
			path=names(gr)[terminal_ordinal.first],
			current=partner(gr)$ordinal[terminal_ordinal.first],
			path_length=1) %>%
		dplyr::select(ordinal, current, path, path_length, terminal_ordinal, distance)

	resultdf = NULL
	while(nrow(activedf) > 0 & max_depth > 0) {
		activedf = activedf %>% dplyr::inner_join(transitive_hits, by=c("current"="transitive_start")) %>%
			mutate(
				distance = distance + transitive_distance,
				current=transitive_end,
				path=paste(path, path_name),
				path_length=path_length + 1) %>%
			dplyr::select(ordinal, current, path, path_length, terminal_ordinal, distance) %>%
			filter(distance <= max_traversal_length + max_positional_error & distance >= -max_positional_error) %>%
			filter(current != ordinal & current != terminal_ordinal) # don't follow loops
		current_terminal = activedf %>%
			dplyr::inner_join(terminal_hits, by=c("current"="current", "terminal_ordinal"="terminal_ordinal")) %>%
			mutate(distance = distance + terminal_distance,
						 name=names(gr)[ordinal]) %>%
			dplyr::select(ordinal, name, path, path_length, distance)
		resultdf = bind_rows(resultdf, current_terminal)
		max_depth = max_depth - 1
	}
	return(resultdf)
}
test_that("transitive_paths", {
	#debugonce(transitive_paths)
	gr = GRanges(seqnames=c("A", "C", "A", "B", "B", "C", "B", "C"),
							 ranges=IRanges(start=c(1, 5000, 1, 1000, 2000, 5100, 2000, 5100), width=1),
							 strand=c("+", "-", "+", "-", "+", "-", "+", "-"))
	names(gr) = c("AC1", "AC2", "AB1", "AB2", "BC1", "BC2", "diffSampleBC1", "diffSampleBC2")
	gr$partner = c("AC2", "AC1", "AB2", "AB1", "BC2", "BC1", "diffSampleBC2", "diffSampleBC1")
	gr$sampleId = c(1, 1, 1, 1, 1, 1, 2, 2)
	paths = transitive_paths(gr, 2000, 200)
	expect_equal(2, nrow(paths))
	expect_equal(rep(1000-100, 2), paths$distance)
	expect_equal(c("AC1", "AC2"), paths$name)
})





