library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(rtracklayer)
library(tidyverse)
library(stringr)
library(testthat)
source("gridss.config.R")

#' sum of genotype fields
.genosum <- function(genotypeField, columns) {
	rowSums(genotypeField[,columns, drop=FALSE])
}
.addFilter = function(existing, filterName, appliesTo) {
  existing[appliesTo] = paste(existing[appliesTo], filterName, sep=";")
  return(existing)
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakpoint_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, support_quality_filters=TRUE, normalOrdinal=1, tumourOrdinal=2,
                                    pon_dir=NULL,
                                    pongr=read_gridss_breakpoint_pon(paste(pon_dir, "gridss_pon_breakpoint.bedpe", sep="/"))) {
	vcf = vcf[names(gr)]
	i = info(vcf)
	g = geno(vcf)
	isShort = is_short_deldup(gr)
	ihomlen = rowSums(abs(as.matrix(info(vcf)$IHOMPOS)))
	ihomlen[is.na(ihomlen)] = 0
	filtered = rep("", length(gr))
	if (!is.null(pongr)) {
	  filtered = .addFilter(filtered, "pon", findBreakpointOverlaps(gr, pongr)$queryHits)
	}
	if (support_quality_filters) {
	  # TODO: update this to a binomial test so we don't filter low confidence
	  # variants that are strand biased by chance
	  # long variants we expect to be heavily strand biased as RP support (including via assembly)
	  # is fully strand biased when originating from one side of a breakend
	  filtered = .addFilter(filtered, "strand_bias", isShort & pmax(i$SB, 1 - i$SB) > gridss.max_allowable_short_event_strand_bias)


	  #filtered = .addFilter(filtered, "FlankingHighQualIndel", str_detect(gridss_gr$FILTER, "SINGLE_ASSEMBLY") & (i$RP + i$SR) / i%VF < 0.1 & i$RP >= 2 & i$SR >= 2 # fixed in GRIDSSv1.8.0
	}
	if (min_support_filters) {
	  filtered = .addFilter(filtered, "af", gridss_somatic_bp_af(gr, vcf) < gridss.min_af)
	  # Multiple biopsy concordance indicates that assemblies with very few supporting reads are sus
	  #filtered = .addFilter(filtered, "minRead", .genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$SR, c(normalOrdinal, tumourOrdinal)) < gridss.min_direct_read_support)
	  # very high coverage hits assembly threshold; we also need to keep transitive calls so we reallocate them to get the correct VF
	  #filtered = .addFilter(filtered, "NO_ASSEMBLY", str_detect(gr$FILTER, "NO_ASSEMBLY"))
	  # homology FPs handled by normal and/or PON
	  #filtered = .addFilter(filtered, "homologyLength", ihomlen > gridss.max_homology_length)
		# Added ASRP into filter otherwise breakpoint chains don't get called
	  filtered = .addFilter(filtered, "BPI.Filter.PRSupportZero", !isShort & .genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$ASRP,c(normalOrdinal, tumourOrdinal)) == 0)
	  filtered = .addFilter(filtered, "BPI.Filter.SRSupportZero", isShort & .genosum(g$SR,c(normalOrdinal, tumourOrdinal)) == 0)
	}
	if (somatic_filters) {
		#normalaf <- gridss_af(gr, vcf, normalOrdinal)
	  filtered = .addFilter(filtered, "normalSupport", .genosum(g$VF,normalOrdinal) > gridss.allowable_normal_contamination * .genosum(g$VF,tumourOrdinal))
	  filtered = .addFilter(filtered, "SRNormalSupport", isShort & .genosum(g$SR, normalOrdinal) != 0)
	  filtered = .addFilter(filtered, "normalCoverage", .genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$VF, normalOrdinal) < gridss.min_normal_depth)
	}
	return(filtered)
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakend_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, normalOrdinal=1, tumourOrdinal=2,
                                  pon_dir=NULL,
                                  pongr=import(paste(pon_dir, "gridss_pon_single_breakend.bed", sep="/"))) {
  vcf = vcf[names(gr)]
  i = info(vcf)
  g = geno(vcf)
  filtered = rep("", length(gr))
  if (!is.null(pongr)) {
    filtered = .addFilter(filtered, "pon", queryHits(findOverlaps(gr, pongr)))
  }
  if (min_support_filters) {
    filtered = .addFilter(filtered, "af", gridss_somatic_be_af(gr, vcf) < gridss.min_af)
    filtered = .addFilter(filtered, "imprecise", i$IMPRECISE)
    filtered = .addFilter(filtered, "NO_ASSEMBLY", str_detect(gr$FILTER, "NO_ASSEMBLY"))
    # require direct read support as well
    #filtered = .addFilter(filtered, "minRead", .genosum(g$BSC,c(normalOrdinal, tumourOrdinal)) + .genosum(g$BUM, c(normalOrdinal, tumourOrdinal)) < gridss.min_direct_read_support * gridss.single_breakend_multiplier)
    # require at least one read pair included in the assembly
    # this is a relatively strict filter but does filter out most of the
    # noise from microsatellite sequences
    filtered = .addFilter(filtered, "NO_ASRP", i$BASRP == 0)
  }
  if (somatic_filters) {
    filtered = .addFilter(filtered, "normalSupport", .genosum(g$BVF,normalOrdinal) > gridss.allowable_normal_contamination * .genosum(g$BVF,tumourOrdinal))
    filtered = .addFilter(filtered, "normalCoverage", .genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$BVF, normalOrdinal) < gridss.min_normal_depth)
  }
  return(filtered)
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
#' @description filter out 'shadow' calls of strong multi-mapping calls
#' bwa overestimates the MAPQ of some multimapping reads
#' and makes FP calls to alternate locations.
is_shadow_breakpoint = function(bpgr, begr, vcf, breakendQualMultiple=3) {
  combinedgr <- bpgr
  combinedgr$partner <- NULL
  combinedgr <- c(combinedgr, begr)
  # require exact overlaps
  bestOverlap = findOverlaps(bpgr, combinedgr) %>%
    as.data.frame() %>%
    # bpgr offsets are the same as combinedgr ofsets
    filter(queryHits != subjectHits) %>%
    mutate(beQUAL = combinedgr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    summarise(beQUAL = max(beQUAL))
  bpgr$overlapQUAL = 0
  bpgr$overlapQUAL[bestOverlap$queryHits] = bestOverlap$beQUAL
  i <- info(vcf[bpgr$vcfId])
  better_call_filter = !is_short_event(bpgr) &
    bpgr$overlapQUAL > gridss.shadow_breakend_multiple * bpgr$QUAL &
    i$BANSRQ + i$BANRPQ > i$ASQ
  return(as.logical(better_call_filter))
}

is_short_event = function(gr) {
  seqnames(gr) == seqnames(partner(gr)) &
    abs(start(gr)-start(partner(gr))) < gridss.short_event_size_threshold
}
gridss_bp_af = function(gr, vcf, i=c(1)) {
  return(.gridss_af(gr, vcf, i, !is_short_deldup(gr), includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE))
}
gridss_somatic_bp_af = function(gr, vcf) {
	return(gridss_bp_af(gr, vcf, 2))
}
gridss_be_af = function(gr, vcf, i=c(1)) {
  return(.gridss_af(gr, vcf, i, TRUE, includeBreakpointSupport=FALSE, includeBreakendSupport=TRUE))
}
gridss_somatic_be_af = function(gr, vcf) {
  return(gridss_be_af(gr, vcf, 2))
}
.gridss_af = function(gr, vcf, i, includeRefPair, no_coverage_af=0, includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE) {
  genotype = geno(vcf[names(gr)])
  g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], i) } else { genotype[[field]] } })
  names(g) <- names(genotype)
  support = rep(0, length(gr))
  if (includeBreakpointSupport) {
    support = support + g$VF
  }
  if (includeBreakendSupport) {
    support = support + g$BVF
  }
  vf_af = support / (support + g$REF + ifelse(!includeRefPair, 0, g$REFPAIR))
  vf_af[is.nan(vf_af)] = no_coverage_af
  return(vf_af)
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
simpleEventType <- function(gr) {
  same_chr = as.logical(seqnames(gr) != seqnames(partner(gr)))
  insSeq = rep("", length(gr))
  if (!is.null(gr$insSeq)) {
    insSeq = gr$insSeq
  }
  if (!is.null(gr$insertSequence)) {
    insSeq = gr$insertSequence
  }
  more_ins_than_length = str_length(insSeq) >= abs(start(gr)-start(partner(gr))) * 0.7
  same_strand = strand(gr) == strand(partner(gr))
  return(ifelse(!same_chr, "ITX", # inter-chromosomosal
          ifelse(more_ins_than_length, "INS", # TODO: improve classification of complex events
            ifelse(same_strand, "INV",
              ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
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

	resultdf = data.frame(ordinal=integer(), name=character(), bp_path=character(), path_length=integer(), distance=integer(), stringsAsFactors=FALSE)
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

  result_df = data.frame(transitive=character(), bp_path=character(), min_length=integer(), max_length=integer())
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
#' Assumes the input is sorted
#' @param is_higher_breakend record is a breakpoint record and is considered the higher of the two breakends.
#' Default check uses GRIDSS notation. TODO: use breakpointRanges() to make a generic default
align_breakpoints <- function(vcf, align=c("centre"), is_higher_breakend=str_detect(names(vcf), "h$")) {
  align = match.arg(align)
  nominal_start = start(rowRanges(vcf))
  if (!all(elementNROWS(info(vcf)$CIPOS) == 2)) {
    stop("CIPOS not specified for all variants.")
  }
  cipos = t(matrix(unlist(info(vcf)$CIPOS), nrow=2))
  ciwdith = cipos[,2] - cipos[,1]
  if (align == "centre") {
    adjust_by = cipos[,1] + ciwdith / 2.0
    adjust_by = ifelse(round(adjust_by) != adjust_by, adjust_by + ifelse(is_higher_breakend, -0.5, 0.5), adjust_by)
  } else {
    stop("Only centre alignment is currently implemented.")
  }
  rowRanges(vcf) = shift(rowRanges(vcf), ifelse(is.na(adjust_by), 0, adjust_by))
  info(vcf)$CIPOS = info(vcf)$CIPOS - adjust_by
  if (!is.null(info(vcf)$CIEND)) {
    info(vcf)$CIEND = info(vcf)$CIEND - adjust_by
  }
  if (!is.null(info(vcf)$IHOMPOS)) {
    info(vcf)$IHOMPOS = info(vcf)$IHOMPOS - adjust_by
  }
  info(vcf)$CIRPOS = NULL # TODO: remove CIRPOS from GRIDSS entirely
  # left align lower breakend (forces right alignment of the partner if they're in the same orientation)
  #info(vcf[names(gr)])$CIPOS = relist(c(rbind(ifelse(left_align, 0, -bp_width), ifelse(left_align, bp_width, 0))), PartitioningByEnd(seq(2, 2*length(gr), 2)))
  return(vcf)
}

readVcf = function(file, ...) {
  raw_vcf = VariantAnnotation::readVcf(file=file, ...)
  # work-around for https://github.com/Bioconductor/VariantAnnotation/issues/8
  alt = read_tsv(file, comment="#", col_names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", seq_len(ncol(geno(raw_vcf)[[1]]))), cols_only(ALT=col_character()))$ALT
  VariantAnnotation::fixed(raw_vcf)$ALT = CharacterList(lapply(as.character(alt), function(x) x))
  return(raw_vcf)
}

read_gridss_breakpoint_pon = function(file) {
  df = read_tsv(file,
                col_names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2", "IMPRECISE", "samples"),
                col_types="ciiciiccccli")
  gro = GRanges(
    seqnames=df$chr1,
    ranges=IRanges(
      start=df$start1 + 1,
      end=df$end1),
    strand=df$strand1,
    partner=paste0(seq_len(nrow(df)), "h"),
    IMPRECISE = df$IMPRECISE,
    samples=df$samples)
  names(gro) = paste0(seq_len(nrow(df)), "o")
  grh = GRanges(
    seqnames=df$chr2,
    ranges=IRanges(
      start=df$start2 + 1,
      end=df$end2),
    strand=df$strand2,
    partner=paste0(seq_len(nrow(df)), "o"),
    IMPRECISE = df$IMPRECISE,
    samples=df$samples)
  names(grh) = paste0(seq_len(nrow(df)), "h")
  return(c(gro, grh))
}
linked_assemblies = function(vcf) {
  # Assembly-based event linking
  asm_linked_df = data.frame(
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
transitive_calls = function(vcf, bpgr) {
  # transitive calling reduction
  transitive_df = transitive_breakpoints(bpgr, min_segment_length=20, report="all")
  transitive_df = transitive_df %>%
    filter(info(vcf[transitive_df$transitive])$IMPRECISE) %>%
    mutate(full_path=bp_path)
  if (nrow(transitive_df) != 0) {
    transitive_df = transitive_df %>%
      separate_rows(bp_path, sep=";") %>%
      mutate(imprecise=info(vcf[bp_path])$IMPRECISE) %>%
      group_by(transitive, full_path) %>%
      summarise(
        imprecise=sum(imprecise),
        hops=n(),
        min_length=min(min_length),
        max_length=max(max_length)) %>%
      # for now we just link imprecise transitive events
      # via precise calls
      filter(imprecise == 0) %>%
      group_by(transitive) %>%
      top_n(1, min_length) %>%
      ungroup() %>%
      dplyr::select(linked_by = transitive, vcfid = full_path) %>%
      separate_rows(vcfid, sep=";")
  } else {
    # work-around for https://github.com/tidyverse/tidyr/issues/470
    transitive_df = data.frame(linked_by=character(), vcfid=character(), stringsAsFactors=FALSE)
  }
}
#' @description Calculates a somatic score for the given variant.
#' The score is defined as the log likelihood ratio of the observed
#' variant read counts being caused by tumour contamination of the normal,
#' and the variant being heterozygous in the normal.
#' @param vcf GRIDSS VCF object to calculate somatic scores for
#' @param normalOrdinal Name or index of the normal genotype in the VCF
#' @param tumourOrdinal Name or index of the tumour genotype in the VCF
#' @param contamination_rate estimated rate of tumour contamination. Typically,
#' this should be set to the tumour variant allele frequency multiplied by the
#' copy number ratio of the tumour and normal.
#' @param normal_ploidy Expected normal ploidy.
gridss_breakpoint_somatic_llr = function(vcf, normalOrdinal=1, tumourOrdinal=2, contamination_rate, normal_ploidy = 2, bpgr=breakpointRanges(vcf)) {
  g = geno(vcf)
  df = data.frame(
    is_bp = names(vcf) %in% names(bpgr),
    is_short_bp = names(vcf) %in% names(bpgr)[is_short_deldup(bpgr)])
  row.names(df) = names(vcf)
  partner_lookup = match(as.character(bpgr$partner), row.names(df))
  names(partner_lookup) = names(bpgr)
  df$partner = NA_integer_
  df[df$is_bp, "partner"] = partner_lookup[row.names(df)[df$is_bp]]
  df = df %>% mutate(
    n_ref = .genosum(g$REF,c(normalOrdinal)) + ifelse(is_short_bp, 0, .genosum(g$REFPAIR,c(normalOrdinal))),
    n_var = ifelse(is_bp, .genosum(g$VF,c(normalOrdinal)), .genosum(g$BVF,c(normalOrdinal))))
  df = df %>% mutate(
    n_ref = ifelse(!is_bp, n_ref, n_ref + df$n_ref[partner]),
    # assume variant ploidy is 1, then:
    germline_het_af = ifelse(is_bp, 1 / (1 + 2*(normal_ploidy - 1)), 1 / normal_ploidy),
    contamination_af = contamination_rate
  )
  df = df %>% mutate(
    germline_het_log_p = dbinom(n_var, n_ref + n_var, germline_het_af, log=TRUE),
    contamination_log_p = dbinom(n_var, n_ref + n_var, contamination_af, log=TRUE)
  )
  return (df$contamination_log_p - df$germline_het_log_p)
}
