import argparse
import logging
import sys
from pathlib import Path
from typing import List

from ref_genome_feature_analysis import ReferenceGenomeFeatureAnalysisConfig, ReferenceGenomeFeatureAnalyzer
from ref_util import set_up_logging


# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.


def main(config: ReferenceGenomeFeatureAnalysisConfig) -> None:
    set_up_logging()
    config.validate()

    analysis = ReferenceGenomeFeatureAnalyzer.do_analysis(config)

    logging.info(f"FEATURES GENOME:")

    logging.info(f"Unplaced contigs: {analysis.has_unplaced_contigs}")
    logging.info(f"Unlocalized contigs: {analysis.has_unlocalized_contigs}")
    logging.info(f"Alts: {analysis.has_alts}")
    logging.info(f"Decoys (hs38d1): {analysis.has_decoys}")
    logging.info(f"Patches: {analysis.has_patches}")
    logging.info(f"EBV: {analysis.has_ebv}")
    logging.info(f"rCRS mitochondrial sequence: {analysis.has_rcrs}")
    logging.info(f"Uses canonical contig names, so 'chr1' etc.: {analysis.uses_canonical_chrom_names}")
    logging.info(f"PAR hardmask (not fully accurate): {analysis.has_only_hardmasked_nucleotides_at_y_par1}")
    logging.info(f"Semi ambiguous IUB codes: {analysis.has_semi_ambiguous_iub_codes}")
    logging.info(f"Has softmasked nucleotides: {analysis.has_softmasked_nucleotides}")
    logging.info(f"Alts are padded with N: {analysis.alts_are_padded}")
    logging.info(f"PhiX: False?")
    logging.info(f"")
    logging.info(f"For easy copy-paste:")
    answers: List[bool] = [
        analysis.has_unplaced_contigs,
        analysis.has_unlocalized_contigs,
        analysis.has_alts,
        analysis.has_decoys,
        analysis.has_patches,
        analysis.has_ebv,
        analysis.has_rcrs,
        analysis.uses_canonical_chrom_names,
        analysis.has_only_hardmasked_nucleotides_at_y_par1,
        analysis.has_semi_ambiguous_iub_codes,
        analysis.has_softmasked_nucleotides,
        analysis.alts_are_padded,
    ]
    value_to_answer = {True: "Yes", False: "No", None: "?"}
    print("\n".join([value_to_answer[answer] for answer in answers]))


def parse_args(sys_args: List[str]) -> ReferenceGenomeFeatureAnalysisConfig:
    parser = argparse.ArgumentParser(
        prog="check_ref_genome_features",
        description=(
            "Check important features for ref genome FASTA file. Correctness is not guaranteed, especially for hg19. "
            "See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files."
        ),
    )
    parser.add_argument("--fasta", "-i", type=Path, required=True, help="Fasta file for ref genome.")
    parser.add_argument(
        "--contig_alias",
        "-c",
        type=str,
        required=True,
        help="Bucket path to TSV file with contig name translations. Source: create_contig_translation_file.py.",
    )

    parser.add_argument(
        "--rcrs",
        "-r",
        type=Path,
        required=False,
        help="(Optional) Fasta file for rCRS mitochondrial genome.",
    )

    args = parser.parse_args(sys_args)

    config = ReferenceGenomeFeatureAnalysisConfig(args.fasta, args.rcrs, args.contig_alias)
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
