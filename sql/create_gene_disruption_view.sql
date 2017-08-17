USE hmfpatients;
DROP VIEW IF EXISTS gene_disruption_view;
CREATE VIEW gene_disruption_view AS
SELECT
    m.sampleId,
    m.gene,
    m.exon_rank,
    m.total_exons,
    m.is_upstream,
    m.is_exonic,
    m.svType,
    m.bp_chromosome,
    m.bp_position,
    m.opposing_bp_chromosome,
    m.opposing_bp_position,
    if (m.bp_chromosome = m.opposing_bp_chromosome,abs(m.opposing_bp_position-m.bp_position),0) as sv_length
FROM (
 
    SELECT
        bp1_gene as gene,
        bp1_exon_rank as exon_rank,
        (SELECT MAX(rank) FROM homo_sapiens_core_89_37.exon_transcript as et WHERE et.transcript_id = bp1_transcript_id GROUP BY et.transcript_id) as total_exons,
        COALESCE(bp1_is_upstream, 1) as is_upstream,
        bp1_is_exonic as is_exonic,
        sv.type as svType,
        sv.sampleId as sampleId,
        startChromosome as bp_chromosome,
        startPosition as bp_position,
        endChromosome as opposing_bp_chromosome,
        endPosition as opposing_bp_position
 
    FROM hmfpatients.sv_annotation2 as ann1
        LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = sv_id
        LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
        LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
    WHERE
        (bp1_gene_id is not null) and
        bp1_is_canonical_transcript and
        EXISTS (SELECT 1 FROM census.ensembl_panel WHERE ensembl_gene_id = bp1_gene_id) and
        NOT EXISTS (SELECT 1 FROM hmfpatients.sv_annotation2 as ann2 WHERE ann2.bp1_gene_id = ann2.bp2_gene_id and ann2.sv_id = ann1.sv_id
        and abs(ann2.bp1_exon_rank-ann2.bp2_exon_rank)=ann2.is_strand_compatible and not ann2.bp1_is_exonic and not ann2.bp2_is_exonic)
UNION
 
        SELECT
            bp2_gene as gene,
            bp2_exon_rank as exon_rank,
            (SELECT MAX(rank) FROM homo_sapiens_core_89_37.exon_transcript as et WHERE et.transcript_id = bp2_transcript_id GROUP BY et.transcript_id) as total_exons,
            COALESCE(bp2_is_upstream, 1) as is_upstream,
            bp2_is_exonic as is_exonic,
            sv.type as svType,
            sv.sampleId as sampleId,
            endChromosome as bp_chromosome,
            endPosition as bp_position,
            startChromosome as opposing_bp_chromosome,
            startPosition as opposing_bp_position
       
    FROM hmfpatients.sv_annotation2 as ann1
        LEFT JOIN hmfpatients.structuralVariant as sv ON sv.id = sv_id
        LEFT JOIN homo_sapiens_core_89_37.exon as ex1 ON ex1.exon_id = bp1_exon_id
        LEFT JOIN homo_sapiens_core_89_37.exon as ex2 on ex2.exon_id = bp2_exon_id
    WHERE
        (bp2_gene_id is not null) and
        bp2_is_canonical_transcript and
        EXISTS (SELECT 1 FROM census.ensembl_panel WHERE ensembl_gene_id = bp2_gene_id) and
        NOT EXISTS (SELECT 1 FROM hmfpatients.sv_annotation2 as ann2 WHERE ann2.bp1_gene_id = ann2.bp2_gene_id and ann2.sv_id = ann1.sv_id
        and abs(ann2.bp1_exon_rank-ann2.bp2_exon_rank)=ann2.is_strand_compatible and not ann2.bp1_is_exonic and not ann2.bp2_is_exonic)
) as m;
