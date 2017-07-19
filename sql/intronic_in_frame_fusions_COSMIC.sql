USE hmfpatients;
SELECT
	bp1_gene,
	bp2_gene,
	count(DISTINCT sv.sampleID) as num_samples,
	group_concat(DISTINCT sv.sampleID SEPARATOR ',') as samples
FROM hmfpatients.sv_annotation2
	LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = sv_id
	LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
	LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
WHERE
	bp1_gene_id is not null and
	bp2_gene_id is not null and
	bp1_gene_id <> bp2_gene_id and
	not bp1_is_exonic and
	not bp2_is_exonic and
	is_strand_compatible and
    # phase compatible
	(CASE WHEN bp1_is_upstream THEN ex1.end_phase = ex2.phase ELSE ex1.phase = ex2.end_phase END) and
    # one end in COSMIC
	(EXISTS (SELECT 1 FROM census.ensembl_census WHERE ensembl_gene_id = bp1_gene_id) OR EXISTS (SELECT 1 FROM census.ensembl_census WHERE ensembl_gene_id = bp2_gene_id))
GROUP BY bp1_gene, bp2_gene
ORDER BY 3 DESC;
