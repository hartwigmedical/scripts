library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(tidyverse)
library(stringr)
library(testthat)
source("gridss.config.R")

#' sum of genotype fields
.genosum <- function(genotypeField, columns) {
	rowSums(genotypeField[,columns, drop=FALSE])
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakpoint_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, support_quality_filters=TRUE, normalOrdinal=1, tumourOrdinal=2) {
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
	  # long variants we expect to be heavily strand biased as RP support (including via assembly)
	  # is fully strand biased when originating from one side of a breakend
	  filtered = filtered | (isShort & pmax(i$SB, 1 - i$SB) > gridss.max_allowable_shot_event_strand_bias)
	}
	if (min_support_filters) {
		filtered = filtered |
			# str_detect(gr$FILTER, "NO_ASSEMBLY") | # very high coverage hits assembly threshold; we also need to keep transitive calls so we reallocate them to get the correct VF
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
		#normalaf <- gridss_af(gr, vcf, normalOrdinal)
		filtered = filtered |
	    .genosum(g$VF,normalOrdinal) > gridss.allowable_normal_contamination *  .genosum(g$VF,tumourOrdinal)
			# Filter.SRNormalSupport
			(isShort & .genosum(g$SR, normalOrdinal) != 0) |
			.genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$VF, normalOrdinal) < gridss.min_normal_depth
	}
	return(as.logical(filtered))
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakend_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, normalOrdinal=1, tumourOrdinal=2) {
  vcf = vcf[names(gr)]
  i = info(vcf)
  g = geno(vcf)
  filtered = rep(FALSE, length(gr))
  if (min_support_filters) {
    filtered = filtered |
      str_detect(gr$FILTER, "NO_ASSEMBLY") |
      str_detect(gr$FILTER, "ASSEMBLY_ONLY") |
      # BPI.Filter.MinDepth
      .genosum(g$BVF,c(normalOrdinal, tumourOrdinal)) < gridss.min_breakpoint_depth |
      # require some sort of breakend anchoring
      i$IMPRECISE |
      # require at least one read pair included in the assembly
      # this is a relatively strict filter but does filter out most of the
      # noise from microsatellite sequences
      i$BASRP == 0
  }
  if (somatic_filters) {
    filtered = filtered |
      .genosum(g$BVF,normalOrdinal) > gridss.allowable_normal_contamination * .genosum(g$BVF,tumourOrdinal) |
    .genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$BVF, normalOrdinal) < gridss.min_normal_depth
  }
  return(as.logical(filtered))
}
is_short_deldup = function(gr) {
  is_deldup = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr <- gr[isbp]
  	bp_short_deldup = strand(bpgr) != strand(partner(bpgr)) &
  		seqnames(bpgr) == seqnames(partner(bpgr)) &
  		abs(start(bpgr)-start(partner(bpgr))) < gridss.short_event_size_threshold
  	is_deldup[isbp] = bp_short_deldup
  }
	return(is_deldup)
}
is_short_event = function(gr) {
  seqnames(gr) == seqnames(partner(gr)) &
    abs(start(gr)-start(partner(gr))) < gridss.short_event_size_threshold
}
gridss_bp_af = function(gr, vcf, i=c(1)) {
	genotype = geno(vcf[names(gr)])
	g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], i) } else { genotype[[field]] } })
	names(g) <- names(genotype)
	ref = g$REF
	refpair = g$REFPAIR
	vf_af = g$VF / (g$VF + ref + ifelse(is_short_deldup(gr), 0, refpair))
	return(vf_af)
}
gridss_somatic_bp_af = function(gr, vcf) {
	return(gridss_bp_af(gr, vcf, 2))
}
gridss_be_af = function(gr, vcf, i=c(1)) {
  genotype = geno(vcf[names(gr)])
  g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], i) } else { genotype[[field]] } })
  names(g) <- names(genotype)
  ref = g$REF
  refpair = g$REFPAIR
  vf_af = g$BVF / (g$BVF + ref + ifelse(is_short_deldup(gr), 0, refpair))
  return(vf_af)
}
gridss_somatic_be_af = function(gr, vcf) {
  return(gridss_be_af(gr, vcf, 2))
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
	result[paths$ordinal] = paths$bp_path
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
			bp_path=names(gr)[terminal_ordinal.first],
			current=partner(gr)$ordinal[terminal_ordinal.first],
			path_length=1) %>%
		dplyr::select(ordinal, current, bp_path, path_length, terminal_ordinal, distance)

	resultdf = NULL
	while(nrow(activedf) > 0 & max_depth > 0) {
		activedf = activedf %>% dplyr::inner_join(transitive_hits, by=c("current"="transitive_start")) %>%
			mutate(
				distance = distance + transitive_distance,
				current=transitive_end,
				bp_path=paste(bp_path, path_name),
				path_length=path_length + 1) %>%
			dplyr::select(ordinal, current, bp_path, path_length, terminal_ordinal, distance) %>%
			filter(distance <= max_traversal_length + max_positional_error & distance >= -max_positional_error) %>%
			filter(current != ordinal & current != terminal_ordinal) # don't follow loops
		current_terminal = activedf %>%
			dplyr::inner_join(terminal_hits, by=c("current"="current", "terminal_ordinal"="terminal_ordinal")) %>%
			mutate(distance = distance + terminal_distance,
						 name=names(gr)[ordinal]) %>%
			dplyr::select(ordinal, name, bp_path, path_length, distance)
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

linked_asm <- function(vcf) {
  asm_linked_df <- data.frame(
    event=info(vcf)$EVENT,
    beid=sapply(info(vcf)$BEID, function(x) paste0(c("", unlist(x)), collapse=";"))) %>%
    separate_rows(beid, sep=";") %>%
    filter(!is.na(beid) & nchar(beid) > 0) %>%
    group_by(event, beid) %>%
    distinct() %>%
    group_by(beid) %>%
    filter(n() > 1) %>%
    mutate(linked_by=beid) %>%
    dplyr::select(event, linked_by=beid)
  return (asm_linked_df)
}

#' @description Determines which A-C transitives calls can be explained by an A-B-C.
#' Transitive calls are caused by fragments that completely span the B fragment.
#' For paired end sequencing, transitive calls will be IMPRECISE as all split reads
#' support the A-B and B-C breakpoints.
#' @param gr input breakpoint GRanges object
#' @param max_traversed_length maximum length of sequence in the traversed breakpoint chain.
#' @param min_segment_length minimum length of traversed sequence. Aligners typically cannot align less than around 20 bases.
#' @param allow_loops allow breakpoints to be traversed multiple times
#' @param max_hops maximum number of breakends to traverse
#' @param report When report is "shortest", only the shortest transitive breakpoint chain is reported.
#' @return id of the
transitive_breakpoints <- function(gr, max_traversed_length=1000, min_segment_length=0, transitive_call_slop=100, allow_loops=FALSE, max_hops=8, report=c("shortest", "all")) {
  ordinal_lookup = seq_len(length(gr))
  names(ordinal_lookup) = names(gr)
  partner_lookup = ordinal_lookup[gr$partner]
  next_df = .traversal_next_breakpoint(gr, max_traversed_length, min_segment_length) %>%
    # prevent self-intersections
    filter(queryHits != subjectHits & queryHits != partner_lookup[subjectHits]) %>%
    # incorporate breakpoint inserted sequence into the traversal length
    mutate(
      min_traversed=min_traversed + gr$insLen[subjectHits],
      max_traversed=max_traversed + gr$insLen[subjectHits],
      source_from=partner_lookup[queryHits],
      source_to=queryHits,
      dest_from=subjectHits,
      dest_to=partner_lookup[subjectHits]) %>%
    dplyr::select(-queryHits, -subjectHits)
  terminal_df = .adjacent_breakends(gr, gr, maxgap=transitive_call_slop, allowed_orientation=c("--", "++")) %>%
    filter(queryHits != subjectHits & queryHits != partner_lookup[subjectHits]) %>%
    # [min_traversed, max_traversed] overlaps [-transitive_call_slop, transitive_call_slop]
    filter(min_traversed < transitive_call_slop & max_traversed > -transitive_call_slop)

  result_df = NULL
  # start with the first traversal away from the putative transitive call
  active_df = terminal_df %>%
    mutate(
      terminal_start=queryHits,
      terminal_end=partner_lookup[queryHits],
      current_to=partner_lookup[subjectHits],
      bp_path=names(gr)[subjectHits],
      min_length=min_traversed + gr$insLen[subjectHits],
      max_length=max_traversed + gr$insLen[subjectHits]) %>%
    dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length)
  i = 0
  while (nrow(active_df) > 0 & i < max_hops) {
    # continue traversing
    active_df = active_df %>%
      dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length) %>%
      inner_join(next_df, by=c("current_to"="source_to")) %>%
      filter(allow_loops | !str_detect(bp_path, stringr::fixed(names(gr)[dest_from]))) %>%
      mutate(
        bp_path=paste0(bp_path, ";", names(gr)[dest_from]),
        current_to=dest_to,
        min_length=min_length + min_traversed,
        max_length=max_length + max_traversed) %>%
      dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length) %>%
      filter(min_length < max_traversed_length)
    # check for terminal completion
    active_df = active_df %>% left_join(terminal_df, by=c("current_to"="queryHits", "terminal_end"="subjectHits"))
    result_df = active_df %>%
      filter(!is.na(min_traversed)) %>%
      mutate(min_length=min_length + min_traversed,
             max_length=max_length + max_traversed,
             transitive=names(gr)[terminal_start]) %>%
      dplyr::select(transitive, bp_path, min_length, max_length) %>%
      bind_rows(result_df)
    if (report == "shortest") {
      # any further solutions will be longer so we don't need to consider them
      active_df = active_df %>% filter(is.na(min_traversed))
      best_min = result_df
    }
    i = i + 1
  }
  if (report == "shortest") {
    # Just the shortest result for each transitive call
    result_df = result_df %>%
      group_by(transitive) %>%
      top_n(1, min_length) %>%
      ungroup()
  }
  return(result_df)
}
.traversal_next_breakpoint <- function(gr, max_traversed_length, min_segment_length) {
  .adjacent_breakends(gr, gr, max_traversed_length, allowed_orientation=c("-+", "+-")) %>%
    filter(max_traversed > min_segment_length) %>%
    mutate(min_traversed=pmax(min_segment_length, min_traversed))
}
#' @description determines which breakends are near the given breakend
#' @param maxgap maximum distance between adjacent breakends
#' @param allowed_orientation allowed breakend orientation of the query then subject breakends
#' @return adjacent breakends of the given orientations.
#' Traversal distances are defined as the number of number of bases traversed into the DNA segment
#' arriving from the query breakend until the subject breakend position is reached.
#' Note that the subject breakend orientation does not affect the traversal distance.
.adjacent_breakends <- function(query, subject, maxgap, allowed_orientation=c("--", "-+", "+-", "++")) {
  allowed_orientation = match.arg(allowed_orientation, several.ok=TRUE)
  findOverlaps(query, subject, maxgap=maxgap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    dplyr::filter(paste0(strand(query[queryHits]), strand(subject[subjectHits])) %in% allowed_orientation) %>%
    dplyr::mutate(
      start_end_traversed=(end(subject[subjectHits]) - start(query[queryHits])) * ifelse(as.logical(strand(query[queryHits])=="-"), 1, -1),
      end_start_traversed=(start(subject[subjectHits]) - end(query[queryHits])) * ifelse(as.logical(strand(query[queryHits])=="-"), 1, -1)) %>%
    dplyr::mutate(min_traversed=pmin(start_end_traversed, end_start_traversed)) %>%
    dplyr::mutate(max_traversed=pmax(start_end_traversed, end_start_traversed)) %>%
    dplyr::select(-start_end_traversed, -end_start_traversed)
}




