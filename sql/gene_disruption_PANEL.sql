USE hmfpatients;
SELECT
	m.gene,
    m.exon_rank,
    m.total_exons,
    m.is_upstream,
    m.is_exonic,
    m.svType,
	count(DISTINCT m.sampleID) as num_samples,
	group_concat(DISTINCT m.sampleID SEPARATOR ',') as samples

FROM (

	SELECT
		bp1_gene as gene,
		bp1_exon_rank as exon_rank,
        (SELECT MAX(rank) FROM homo_sapiens_core_89_37.exon_transcript as et WHERE et.transcript_id = bp1_transcript_id GROUP BY et.transcript_id) as total_exons,
		COALESCE(bp1_is_upstream, 1) as is_upstream,
        bp1_is_exonic as is_exonic,
		sv.type as svType,
		sv.sampleId as sampleId
	FROM hmfpatients.sv_annotation2 as ann1
		LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = sv_id
		LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
		LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
	WHERE
		(bp1_gene_id is not null) and
		bp1_is_canonical_transcript and
		not (bp1_gene_id = bp2_gene_id and abs(bp1_exon_rank-bp2_exon_rank)=is_strand_compatible and not bp1_is_exonic and not bp2_is_exonic) and   # within 1 intron
		EXISTS (SELECT 1 FROM census.ensembl_panel WHERE ensembl_gene_id = bp1_gene_id) and
        NOT EXISTS (SELECT 1 FROM hmfpatients.sv_annotation2 as ann2 WHERE ann2.bp2_gene_id = ann1.bp1_gene_id)
		
UNION

	SELECT
		bp2_gene as gene,
		bp2_exon_rank as exon_rank,
        (SELECT MAX(rank) FROM homo_sapiens_core_89_37.exon_transcript as et WHERE et.transcript_id = bp2_transcript_id GROUP BY et.transcript_id) as total_exons,
		COALESCE(bp2_is_upstream, 1) as is_upstream,
        bp2_is_exonic as is_exonic,
		sv.type as svType,
		sv.sampleId as sampleId
	FROM hmfpatients.sv_annotation2 as ann1
		LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = sv_id
		LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
		LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
	WHERE
		(bp2_gene_id is not null) and
		bp2_is_canonical_transcript and
		not (bp1_gene_id = bp2_gene_id and abs(bp1_exon_rank-bp2_exon_rank)=is_strand_compatible and not bp1_is_exonic and not bp2_is_exonic) and   # within 1 intron
		EXISTS (SELECT 1 FROM census.ensembl_panel WHERE ensembl_gene_id = bp2_gene_id) and
        NOT EXISTS (SELECT 1 FROM hmfpatients.sv_annotation2 as ann2 WHERE ann2.bp1_gene_id = ann1.bp2_gene_id)

) as m

GROUP BY 1,2,3,4,5,6
ORDER BY 7 DESC
