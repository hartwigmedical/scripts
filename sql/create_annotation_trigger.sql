DROP TRIGGER trigger_sv_annotation;
CREATE TRIGGER trigger_sv_annotation AFTER INSERT ON structuralVariant
FOR EACH ROW
INSERT INTO sv_annotation2
SELECT
	m.*,
    COALESCE(et1.rank, 0) as bp1_exon_rank,
    COALESCE(et2.rank, 0) as bp2_exon_rank,
    (transcript1.transcript_id = gene1.canonical_transcript_id) as bp1_is_canonical_transcript,
    (transcript2.transcript_id = gene2.canonical_transcript_id) as bp2_is_canonical_transcript,
    (ex1.seq_region_strand * startOrientation = ex2.seq_region_strand * endOrientation) as is_strand_compatible,
    (ex1.seq_region_strand * startOrientation > 0) as bp1_is_upstream,
    (ex2.seq_region_strand * endOrientation < 0) as bp2_is_upstream
FROM (
SELECT 
    sv.id AS sv_id,
    gene1.gene_id AS bp1_gene_id,
    gene2.gene_id AS bp2_gene_id,
    xref1.display_label as bp1_gene,
    xref2.display_label as bp2_gene,
    (CASE WHEN sv.startOrientation > 0 THEN ex1_left.exon_id ELSE ex1_right.exon_id END) AS bp1_exon_id,
    (CASE WHEN sv.endOrientation > 0 THEN ex2_right.exon_id ELSE ex2_left.exon_id END) AS bp2_exon_id,
    COALESCE(ex1_left.exon_id = ex1_right.exon_id, FALSE) AS bp1_is_exonic,
    COALESCE(ex2_left.exon_id = ex2_right.exon_id, FALSE) AS bp2_is_exonic,
    transcript1.transcript_id AS bp1_transcript_id,
    transcript2.transcript_id AS bp2_transcript_id
FROM
    hmfpatients.structuralVariant AS sv
    # get the chromosome name
	INNER JOIN homo_sapiens_core_89_37.seq_region AS seq1 ON seq1.name = sv.startChromosome AND seq1.coord_system_id = 2
	INNER JOIN homo_sapiens_core_89_37.seq_region AS seq2 ON seq2.name = sv.endChromosome AND seq2.coord_system_id = 2
    # get the genes
	LEFT JOIN homo_sapiens_core_89_37.gene AS gene1
		ON gene1.seq_region_id = seq1.seq_region_id
			AND gene1.status = 'KNOWN'
            AND CASE
				WHEN gene1.seq_region_strand > 0 THEN gene1.seq_region_start + 10000 <= sv.startPosition
                ELSE gene1.seq_region_start <= sv.startPosition END
			AND CASE
				WHEN gene1.seq_region_strand < 0 THEN gene1.seq_region_end + 10000 >= sv.startPosition
                ELSE gene1.seq_region_end >= sv.startPosition END
	LEFT JOIN homo_sapiens_core_89_37.gene AS gene2 
		ON gene2.seq_region_id = seq2.seq_region_id
			AND gene2.status = 'KNOWN'
            AND CASE
				WHEN gene2.seq_region_strand > 0 THEN gene2.seq_region_start + 10000 <= sv.endPosition
                ELSE gene2.seq_region_start <= sv.endPosition END
			AND CASE
				WHEN gene2.seq_region_strand < 0 THEN gene2.seq_region_end + 10000 >= sv.endPosition
                ELSE gene2.seq_region_end >= sv.endPosition END
	# get the gene name
	LEFT JOIN homo_sapiens_core_89_37.xref AS xref1
		ON gene1.display_xref_id = xref1.xref_id
	LEFT JOIN homo_sapiens_core_89_37.xref AS xref2
		ON gene2.display_xref_id = xref2.xref_id
	LEFT JOIN homo_sapiens_core_89_37.transcript as transcript1
		ON transcript1.gene_id = gene1.gene_id and transcript1.status = 'KNOWN'
	LEFT JOIN homo_sapiens_core_89_37.transcript as transcript2
		ON transcript2.gene_id = gene2.gene_id and transcript2.status = 'KNOWN'
	# get the exons
	LEFT JOIN homo_sapiens_core_89_37.exon AS ex1_left
		ON ex1_left.exon_id = (
			SELECT exon.exon_id FROM homo_sapiens_core_89_37.exon_transcript AS et 
				INNER JOIN homo_sapiens_core_89_37.exon AS exon ON et.exon_id = exon.exon_id
				WHERE et.transcript_id = transcript1.transcript_id
					AND sv.startPosition >= exon.seq_region_start
				ORDER BY exon.seq_region_start DESC
                LIMIT 1
		)
	LEFT JOIN homo_sapiens_core_89_37.exon AS ex1_right
		ON ex1_right.exon_id = (
			SELECT exon.exon_id FROM homo_sapiens_core_89_37.exon_transcript AS et 
				INNER JOIN homo_sapiens_core_89_37.exon AS exon ON et.exon_id = exon.exon_id
				WHERE et.transcript_id = transcript1.transcript_id
					AND sv.startPosition <= exon.seq_region_end
				ORDER BY exon.seq_region_end ASC
                LIMIT 1
		)
	LEFT JOIN homo_sapiens_core_89_37.exon AS ex2_left
		ON ex2_left.exon_id = (
			SELECT exon.exon_id FROM homo_sapiens_core_89_37.exon_transcript AS et 
				INNER JOIN homo_sapiens_core_89_37.exon AS exon ON et.exon_id = exon.exon_id
				WHERE et.transcript_id = transcript2.transcript_id
					AND sv.endPosition >= exon.seq_region_start
				ORDER BY exon.seq_region_start DESC
                LIMIT 1
		)
	LEFT JOIN homo_sapiens_core_89_37.exon AS ex2_right
		ON ex2_right.exon_id = (
			SELECT exon.exon_id FROM homo_sapiens_core_89_37.exon_transcript AS et 
				INNER JOIN homo_sapiens_core_89_37.exon AS exon ON et.exon_id = exon.exon_id
				WHERE et.transcript_id = transcript2.transcript_id
					AND sv.endPosition <= exon.seq_region_end
				ORDER BY exon.seq_region_end ASC
                LIMIT 1
		)
	WHERE sv.id = NEW.id
) as m

LEFT JOIN structuralVariant AS sv ON sv.id = sv_id
LEFT JOIN homo_sapiens_core_89_37.gene AS gene1 ON gene1.gene_id = m.bp1_gene_id
LEFT JOIN homo_sapiens_core_89_37.gene AS gene2 ON gene2.gene_id = m.bp2_gene_id
LEFT JOIN homo_sapiens_core_89_37.transcript AS transcript1 ON transcript1.transcript_id = m.bp1_transcript_id
LEFT JOIN homo_sapiens_core_89_37.transcript AS transcript2 ON transcript2.transcript_id = m.bp2_transcript_id
LEFT JOIN homo_sapiens_core_89_37.exon AS ex1 ON ex1.exon_id = m.bp1_exon_id
LEFT JOIN homo_sapiens_core_89_37.exon AS ex2 ON ex2.exon_id = m.bp2_exon_id
LEFT JOIN homo_sapiens_core_89_37.exon_transcript AS et1 ON et1.transcript_id = transcript1.transcript_id AND et1.exon_id = ex1.exon_id
LEFT JOIN homo_sapiens_core_89_37.exon_transcript AS et2 ON et2.transcript_id = transcript2.transcript_id AND et2.exon_id = ex2.exon_id
;