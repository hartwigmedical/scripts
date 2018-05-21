library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(dplyr)
library(stringr)
library(testthat)
gridss.short_event_size_threshold = 1000
# 0.5% of the supporting fragments can come from the normal
gridss.allowable_normal_contamination=0.005
gridss.use_read_threshold=3
gridss.min_breakpoint_depth = 5 # Half BPI default
gridss.min_normal_depth = 8
gridss.max_homology_length = 50
gridss.max_allowable_strand_bias = 0.95

#' sum of genotype fields
.genosum <- function(genotypeField, columns) {
	rowSums(genotypeField[,columns, drop=FALSE])
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, support_quality_filters=TRUE, normalOrdinal=1, tumourOrdinal=2) {
	vcf = vcf[names(gr)]
	i = info(vcf)
	g = geno(vcf)
	isShort = is_short_deldup(gr)
	ihomlen = rowSums(abs(as.matrix(info(vcf)$IHOMPOS)))
	ihomlen[is.na(ihomlen)] = 0
	filtered = rep(FALSE, length(gr))
	if (support_quality_filters) {
	  # TODO: update this to a binomial test so we don't filter low confidence
	  # variants that are strand biased by chance
	  filtered = filtered | pmax(i$SB, 1 - i$SB) > gridss.max_allowable_strand_bias
	}
	if (min_support_filters) {
		filtered = filtered |
			# str_detect(gr$FILTER, "NO_ASSEMBLY") | # very high coverage hits assembly threshold
			str_detect(gr$FILTER, "LOW_QUAL") | # exactly the same as QUAL >= 500
			# ihomlen > gridss.max_homology_length | # homology FPs handled by normal and/or PON
			# BPI.Filter.MinDepth
			(.genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$SR, c(normalOrdinal, tumourOrdinal)) < gridss.min_breakpoint_depth) |
			# BPI.Filter.PRSupportZero
			# Added ASRP into filter otherwise breakpoint chains don't get called
			(!isShort & .genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$ASRP,c(normalOrdinal, tumourOrdinal)) == 0) |
			# BPI.Filter.SRSupportZero
			(isShort & .genosum(g$SR,c(normalOrdinal, tumourOrdinal)) == 0)
	}
	if (somatic_filters) {
		# TODO fix https://github.com/PapenfussLab/gridss/issues/114
		normalaf <- gridss_af(gr, vcf, normalOrdinal)
		filtered = filtered |
	    .genosum(g$VF,normalOrdinal) > gridss.allowable_normal_contamination *  .genosum(g$VF,tumourOrdinal)
			# Filter.SRNormalSupport
			(isShort & .genosum(g$SR, normalOrdinal) != 0) |
			.genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$VF) < gridss.min_normal_depth
	}
	return(as.logical(filtered))
}
is_short_deldup = function(gr) {
	strand(gr) != strand(partner(gr)) &
		seqnames(gr) == seqnames(partner(gr)) &
		abs(start(gr)-start(partner(gr))) < gridss.short_event_size_threshold
}
is_short_event = function(gr) {
  seqnames(gr) == seqnames(partner(gr)) &
    abs(start(gr)-start(partner(gr))) < gridss.short_event_size_threshold
}
gridss_af = function(gr, vcf, i=c(1)) {
	genotype = geno(vcf[names(gr)])
	g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], i) } else { genotype[[field]] } })
	names(g) <- names(genotype)
	ref = g$REF
	refpair = g$REFPAIR
	# assembly counts must be halved for reciprocal assemblies since
	# we independently assemble the same reads on both sides
	assembled_read_weight = ifelse(str_detect(gr$FILTER, "SINGLE_ASSEMBLY"), 1, 0.5)
	assrrp_af = ((g$ASSR + g$ASRP) * assembled_read_weight)/((g$ASSR + g$ASRP) * assembled_read_weight + ref + refpair)
	srrp_af = (g$SR + g$RP) / (g$SR + g$RP + ref + refpair)
	assr_af = (g$ASSR * assembled_read_weight) / (g$ASSR * assembled_read_weight + ref)
	# it is questionably whether we should incoproate BSC.
	# It's definitely incorrect in the case of promiscuous breakpoints
	sr_af = (g$SR+g$BSC)/(g$SR+g$BSC+ref)

	# a complete assembly will give the most accurate AF, but
	# a fragmented assembly will underestimate
	component_af = ifelse(is_short_deldup(gr), pmax(assr_af, sr_af), pmax(assrrp_af, srrp_af))

	vf_af = g$VF / (g$VF + ref + ifelse(is_short_deldup(gr), 0, refpair))
	return(vf_af)
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
		partner=paste0(df$id, "h")[seq_along(df$id)],
		name=paste0(df$id, "o")[seq_along(df$id)],
		svlen=df$endPosition-df$startPosition)
	names(gro)=gro$name
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
		partner=paste0(df$id, "o")[seq_along(df$id)],
		name=paste0(df$id, "h")[seq_along(df$id)],
		svlen=df$endPosition-df$startPosition)
	names(grh)=grh$name
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
full_gridss_annotate_gr = function(gr, vcf, geno_suffix=c(".normal", ".tumour")) {
	vcf = vcf[gr$vcfId,]
	g = geno(vcf)
	extra_df <- as.data.frame(as.data.frame(info(vcf)))
	for (i in seq_along(geno_suffix)) {
		if (!is.na(geno_suffix[i])) {
			gdf = data.frame(QUAL=g$QUAL[,i])
			for (col in names(g)) {
				gdf[[paste0(col, geno_suffix[i])]] = g[[col]][,i]
			}
			extra_df <- bind_cols(extra_df, gdf)
		}
	}
	mcols(gr) = bind_cols(as.data.frame(mcols(gr)), extra_df)
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

get_db_comparision_df = function(dbExisting, dbNew, suffix=c(".old", ".new"), sampleIds=NULL, line_annotation_bed=NULL) {
  if (is.null(sampleIds)) {
    common_sample_ids = query_structural_variants_samples(dbExisting)
    common_sample_ids <- common_sample_ids[common_sample_ids %in% query_structural_variants_samples(dbNew)]
  } else {
    common_sample_ids = sampleIds
  }
  grex <- query_structural_variants_for_sample_as_granges(dbExisting, common_sample_ids)
  grnew <- query_structural_variants_for_sample_as_granges(dbNew, common_sample_ids)
  grex$transitive=transitive_of(grex, 2000, 300)
  grnew$transitive=transitive_of(grnew, 2000, 300)
  # make a proxy QUAL since gr_join_to_df needs it to resolve matches in favour of the 'better' one
  grex$QUAL <- ifelse(is.na(grex$ploidy), grex$af, grex$ploidy)
  grnew$QUAL <- ifelse(is.na(grnew$ploidy), grnew$af, grnew$ploidy)
  simpleEventType <- function(gr) {
    return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
                  ifelse(str_length(gr$insertSequence) >= abs(gr$svlen) * 0.7, "INS", # TODO: improve classification of complex events
                         ifelse(strand(gr) == strand(partner(gr)), "INV",
                                ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                       "DUP")))))
  }
  grex$type = simpleEventType(grex)
  grnew$type = simpleEventType(grnew)

  if (!is.null(line_annotation_bed)) {
    grex$isline <- overlapsAny(grex, line_annotation_bed, maxgap=10000)
    grnew$isline <- overlapsAny(grnew, line_annotation_bed, maxgap=10000)
  }
  fullmatchdf <- gr_join_to_df(grex, NULL, grnew, NULL, maxgap=500, sizemargin=2, suffix=suffix)
  fullmatchdf = fullmatchdf %>%
    mutate(called=c("Neither", "Existing", "New", "Both")
           [1+ifelse(!is.na(fullmatchdf[[paste0("start", suffix[1])]]), 1, 0) + ifelse(!is.na(fullmatchdf[[paste0("start", suffix[2])]]), 2, 0)]) %>%
    mutate(transitive=c("Neither", "Existing", "New", "Both")
           [1+ifelse(!is.na(fullmatchdf[[paste0("transitive", suffix[1])]]), 1, 0) + ifelse(!is.na(fullmatchdf[[paste0("transitive", suffix[2])]]), 2, 0)]) %>%
    mutate(type=ifelse(!is.na(fullmatchdf[[paste0("type", suffix[1])]]), fullmatchdf[[paste0("type", suffix[1])]], fullmatchdf[[paste0("type", suffix[2])]]),
           sampleId=ifelse(!is.na(fullmatchdf[[paste0("sampleId", suffix[1])]]), fullmatchdf[[paste0("sampleId", suffix[1])]], fullmatchdf[[paste0("sampleId", suffix[2])]]),
           orientation=ifelse(!is.na(fullmatchdf[[paste0("orientation", suffix[1])]]), fullmatchdf[[paste0("orientation", suffix[1])]], fullmatchdf[[paste0("orientation", suffix[2])]])) %>%
    mutate(svlen=ifelse(!is.na(fullmatchdf[[paste0("svlen", suffix[2])]]), fullmatchdf[[paste0("svlen", suffix[2])]], fullmatchdf[[paste0("svlen", suffix[1])]]))
  if (!is.null(line_annotation_bed)) {
    line1 <- fullmatchdf[[paste0("isline", suffix[1])]]
    line2 <- fullmatchdf[[paste0("isline", suffix[2])]]
    fullmatchdf = fullmatchdf %>% mutate(isline=(!is.na(line1) & line1) | (!is.na(line2) & line2))
  }
  return(fullmatchdf)
}



