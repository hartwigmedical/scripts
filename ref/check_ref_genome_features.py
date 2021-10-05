import argparse
import concurrent.futures
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple, Tuple, Optional, Set

import pysam

DESIRED_AUTOSOME_CONTIG_NAMES = {f"chr{n}" for n in range(1, 23)}
AUTOSOME_CONTIG_NAMES = DESIRED_AUTOSOME_CONTIG_NAMES.union({str(n) for n in range(1, 23)})
X_CHROMOSOME_CONTIG_NAMES = {"chrX", "X"}
Y_CHROMOSOME_CONTIG_NAMES = {"chrY", "Y"}
MITOCHONDRIAL_CONTIG_NAMES = {"chrM", "M", "MT"}
EBV_CONTIG_NAMES = {"chrEBV"}

DECOY_SUFFIX = "_decoy"
UNLOCALIZED_SUFFIX = "_random"
UNPLACED_PREFIX = "chrUn_"

STANDARD_NUCLEOTIDES = {"A", "C", "G", "T", "N"}
SOFTMASKED_NUCLEOTIDES = {"a", "c", "g", "t", "n"}
UNKNOWN_NUCLEOTIDES = {"N", "n"}
Y_PAR1_TEST_REGION = (20000, 2640000)  # this region lies within the Y PAR1 for both GRCh37 and GRCh38


class Config(NamedTuple):
    ref_genome_path: Path
    rcrs_path: Path

    def validate(self) -> None:
        assert_file_exists(self.ref_genome_path)
        assert_file_exists(self.rcrs_path)


class CategorizedContigNames(NamedTuple):
    autosomes: Tuple[str]
    x: str
    y: str
    mitochondrial: str
    ebv: Optional[str]
    decoys: Tuple[str]
    unlocalized_contigs: Tuple[str]
    unplaced_contigs: Tuple[str]
    uncategorized_contigs: Tuple[str]


def main(config: Config) -> None:
    set_up_logging()
    config.validate()

    categorized_contig_names = get_categorized_contig_names(config)
    logging.info(categorized_contig_names)

    has_rcrs = mitochondrial_sequence_is_rcrs(config, categorized_contig_names.mitochondrial)
    y_test_nucleotides = get_nucleotides_from_string(get_y_test_sequence(config, categorized_contig_names.y))
    logging.info(f"nucleotides at y par1 test region: {y_test_nucleotides}")

    nucleotides_file = Path(f"{config.ref_genome_path}.nucleotides")
    if not nucleotides_file.exists():
        nucleotides = get_nucleotides(config)
        nucleotides_string = "".join(sorted(list(nucleotides)))
        with open(nucleotides_file, "w+") as f:
            f.write(nucleotides_string)

    with open(nucleotides_file, "r+") as f:
        nucleotides = set(f.readline().split(""))

    logging.info(f"nucleotides: {nucleotides}")

    uses_desired_chrom_names = bool(set(categorized_contig_names.autosomes).intersection(DESIRED_AUTOSOME_CONTIG_NAMES))
    has_only_hardmasked_nucleotides_at_y_par1 = not bool(y_test_nucleotides.difference(UNKNOWN_NUCLEOTIDES))
    has_semi_ambiguous_iub_codes = bool(nucleotides.difference(STANDARD_NUCLEOTIDES).difference(SOFTMASKED_NUCLEOTIDES))
    has_softmasked_nucleotides = bool(nucleotides.intersection(SOFTMASKED_NUCLEOTIDES))

    logging.info(f"FEATURES GENOME:")
    logging.info(f"Unplaced contigs: {bool(categorized_contig_names.unplaced_contigs)}")
    logging.info(f"Unlocalized contigs: {bool(categorized_contig_names.unlocalized_contigs)}")
    logging.info(f"ALTS: ?")
    logging.info(f"rCRS mitochondrial sequence: {has_rcrs}")
    logging.info(f"Accession numbers: False?")
    logging.info(f"Uses 'chr1' chrom names: {uses_desired_chrom_names}")
    logging.info(f"PAR hardmask (not fully accurate): {has_only_hardmasked_nucleotides_at_y_par1}")
    logging.info(f"Decoys (hs38d1): {bool(categorized_contig_names.decoys)}")
    logging.info(f"EBV: {categorized_contig_names.ebv is not None}")
    logging.info(f"Patches: ?")
    logging.info(f"PhiX: ?")
    logging.info(f"Semi ambiguous IUB codes: {has_semi_ambiguous_iub_codes}")
    logging.info(f"Has softmasked nucleotides: {has_softmasked_nucleotides}")


def get_categorized_contig_names(config: Config) -> CategorizedContigNames:
    with pysam.Fastafile(config.ref_genome_path) as genome_f:
        contig_names = list(genome_f.references)

    autosome_contigs: List[str] = []
    x_contigs: List[str] = []
    y_contigs: List[str] = []
    mitochondrial_contigs: List[str] = []
    ebv_contigs: List[str] = []
    decoy_contigs: List[str] = []
    unlocalized_contigs: List[str] = []
    unplaced_contigs: List[str] = []
    uncategorized_contigs: List[str] = []

    for contig_name in contig_names:
        if contig_name in AUTOSOME_CONTIG_NAMES:
            autosome_contigs.append(contig_name)
        elif contig_name in X_CHROMOSOME_CONTIG_NAMES:
            x_contigs.append(contig_name)
        elif contig_name in Y_CHROMOSOME_CONTIG_NAMES:
            y_contigs.append(contig_name)
        elif contig_name in MITOCHONDRIAL_CONTIG_NAMES:
            mitochondrial_contigs.append(contig_name)
        elif contig_name in EBV_CONTIG_NAMES:
            ebv_contigs.append(contig_name)
        elif contig_name[-6:] == DECOY_SUFFIX and contig_name[:6] == UNPLACED_PREFIX:
            decoy_contigs.append(contig_name)
        elif contig_name[-7:] == UNLOCALIZED_SUFFIX:
            unlocalized_contigs.append(contig_name)
        elif contig_name[:6] == UNPLACED_PREFIX and contig_name[-6:] != DECOY_SUFFIX:
            unplaced_contigs.append(contig_name)
        else:
            uncategorized_contigs.append(contig_name)

    if uncategorized_contigs:
        raise ValueError(f"Uncategorized contigs: {uncategorized_contigs}")
    if len(mitochondrial_contigs) != 1:
        error_msg = (
            f"Did not find exactly one mitochondrial contig: "
            f"mitochondrial={mitochondrial_contigs}, uncategorized={uncategorized_contigs}"
        )
        raise ValueError(error_msg)
    if len(x_contigs) != 1:
        error_msg = (
            f"Did not find exactly one X contig: "
            f"x={x_contigs}, uncategorized={uncategorized_contigs}"
        )
        raise ValueError(error_msg)
    if len(y_contigs) != 1:
        error_msg = (
            f"Did not find exactly one Y contig: "
            f"y={y_contigs}, uncategorized={uncategorized_contigs}"
        )
        raise ValueError(error_msg)
    if len(ebv_contigs) > 1:
        raise ValueError(f"Found more than one EBV contig: {ebv_contigs}")

    categorized_contigs = CategorizedContigNames(
        tuple(autosome_contigs),
        x_contigs[0],
        y_contigs[0],
        mitochondrial_contigs[0],
        ebv_contigs[0] if len(ebv_contigs) == 1 else None,
        tuple(decoy_contigs),
        tuple(unlocalized_contigs),
        tuple(unplaced_contigs),
        tuple(uncategorized_contigs),
    )
    return categorized_contigs


def mitochondrial_sequence_is_rcrs(config: Config, mitochondrial_contig_name: str) -> bool:
    with pysam.Fastafile(config.rcrs_path) as rcrs_f:
        if rcrs_f.nreferences != 1:
            raise ValueError(f"rCRS FASTA file has more than one contig")
        rcrs_genome = rcrs_f.fetch(rcrs_f.references[0])
    with pysam.Fastafile(config.ref_genome_path) as genome_f:
        mitochondrial_from_ref = genome_f.fetch(mitochondrial_contig_name)
    return rcrs_genome == mitochondrial_from_ref


def get_y_test_sequence(config: Config, y_contig_name: str) -> str:
    with pysam.Fastafile(config.ref_genome_path) as genome_f:
        return genome_f.fetch(y_contig_name, Y_PAR1_TEST_REGION[0], Y_PAR1_TEST_REGION[1])


def get_nucleotides(config: Config) -> Set[str]:
    futures = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        with pysam.Fastafile(config.ref_genome_path) as genome_f:
            for contig_name in genome_f.references:
                contig = genome_f.fetch(contig_name)
                futures.append(executor.submit(get_nucleotides_from_string, contig))

    nucleotides = set()
    for future in futures:
        try:
            nucleotides = nucleotides.union(future.result())
        except Exception as exc:
            raise ValueError(exc)
    return nucleotides


def get_nucleotides_from_string(sequence: str) -> Set[str]:
    nucleotides = set()
    for nucleotide in sequence:
        nucleotides.add(nucleotide)
    return nucleotides


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def delete_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()


def assert_file_exists(path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"File does not exist: {path}")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="check_hg38_ref_genome_features",
        description=(
            "Check important features for ref genome FASTA file. Correctness is not guaranteed, especially for hg19."
        ),
    )
    parser.add_argument("--fasta", "-i", type=str, required=True, help="Fasta file for ref genome.")
    parser.add_argument("--rcrs", "-r", type=str, required=True, help="Fasta file for rCRS mitochondrial genome.")

    args = parser.parse_args(sys_args)

    config = Config(Path(args.fasta), Path(args.rcrs))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
