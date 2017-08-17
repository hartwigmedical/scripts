USE hmfpatients;
SELECT
	bp1_gene,
    bp1_exon_rank,
    bp2_exon_rank,
    bp1_is_exonic,
    bp2_is_exonic,
    sv.type,
    count(DISTINCT sampleId) as sampleCount,
    GROUP_CONCAT(DISTINCT sampleId SEPARATOR ',') as samples,
    round(avg(sv.endPosition-sv.startPosition),0) as avgSVlength
FROM sv_annotation2 as ann
LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = ann.sv_id
LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
WHERE
	bp1_gene_id IS NOT NULL and
	bp1_transcript_id = bp2_transcript_id and # same gene and same transcript
    bp1_is_canonical_transcript and
    bp1_exon_rank > 0 and bp2_exon_rank > 0 and # NOTE: this will skip exon 1 being deleted -- what is the impact w.r.t start codons etc
    # in phase
    (CASE WHEN bp1_is_upstream THEN ex1.end_phase = ex2.phase ELSE ex1.phase = ex2.end_phase END) and
    # gene is cosmic
	(EXISTS (SELECT 1 FROM census.ensembl_census WHERE ensembl_gene_id = bp1_gene_id))
GROUP BY 1,2,3,4,5,6
ORDER BY 7 DESC;
