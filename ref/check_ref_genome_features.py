import argparse
import concurrent.futures
import logging
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import List, NamedTuple, Tuple, Optional, Set, Dict

from google.cloud import storage
import pysam

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

DESIRED_AUTOSOME_CONTIG_NAMES = {f"chr{n}" for n in range(1, 23)}
X_CHROMOSOME_CONTIG_NAMES = {"chrX"}
Y_CHROMOSOME_CONTIG_NAMES = {"chrY"}
MITOCHONDRIAL_CONTIG_NAMES = {"chrM"}
EBV_CONTIG_NAMES = {"chrEBV"}

UNPLACED_PREFIX = "chrUn_"
DECOY_SUFFIX = "_decoy"
UNLOCALIZED_SUFFIX = "_random"
FIX_PATCH_SUFFIX = "_fix"
NOVEL_PATCH_SUFFIX = "_novel"
ALT_SUFFIX = "_alt"

STANDARD_NUCLEOTIDES = {"A", "C", "G", "T", "N"}
SOFTMASKED_NUCLEOTIDES = {"a", "c", "g", "t", "n"}
UNKNOWN_NUCLEOTIDES = {"N", "n"}
Y_PAR1_TEST_REGION = (20000, 2640000)  # this region lies within the Y PAR1 for both GRCh37 and GRCh38


class Config(NamedTuple):
    ref_genome_path: Path
    rcrs_path: Path
    contig_alias_bucket_path: str

    def validate(self) -> None:
        assert_file_exists(self.ref_genome_path)
        assert_file_exists(self.rcrs_path)
        if not re.fullmatch(r"gs://.+", self.contig_alias_bucket_path):
            raise ValueError(f"Contig alias bucket path is not of the form 'gs://some/file/path'")
        assert_file_exists_in_bucket(self.contig_alias_bucket_path)


class CategorizedContigNames(NamedTuple):
    autosomes: Tuple[str]
    x: str
    y: str
    mitochondrial: str
    ebv: Optional[str]
    decoys: Tuple[str]
    unlocalized_contigs: Tuple[str]
    unplaced_contigs: Tuple[str]
    alt_contigs: Tuple[str]
    fix_patch_contigs: Tuple[str]
    novel_patch_contigs: Tuple[str]
    uncategorized_contigs: Tuple[str]


class ContigNameTranslator(object):
    """Standardizes names if it can. Returns argument as is if it cannot."""
    def __init__(self, contig_name_to_canonical_name: Dict[str, str]) -> None:
        self._contig_name_to_canonical_name = deepcopy(contig_name_to_canonical_name)

    def standardize(self, contig_name: str) -> str:
        if contig_name in self._contig_name_to_canonical_name:
            return self._contig_name_to_canonical_name[contig_name]
        else:
            return contig_name


def main(config: Config) -> None:
    set_up_logging()
    config.validate()

    contig_name_translator = get_contig_name_translator(config.contig_alias_bucket_path)

    categorized_contig_names = get_categorized_contig_names(config.ref_genome_path, contig_name_translator)
    logging.debug(categorized_contig_names)

    has_rcrs = mitochondrial_sequence_is_rcrs(config, categorized_contig_names.mitochondrial)
    y_test_nucleotides = get_nucleotides_from_string(
        get_y_test_sequence(config.ref_genome_path, categorized_contig_names.y)
    )
    logging.info(f"nucleotides at y par1 test region: {y_test_nucleotides}")

    nucleotides_file = Path(f"{config.ref_genome_path}.nucleotides")
    if not nucleotides_file.exists():
        nucleotides = get_nucleotides(config.ref_genome_path)
        nucleotides_string = "".join(sorted(list(nucleotides)))
        with open(nucleotides_file, "w+") as f:
            f.write(nucleotides_string)

    with open(nucleotides_file, "r+") as f:
        line = f.readline()
        nucleotides = {char for char in line}

    logging.info(f"nucleotides: {nucleotides}")

    uses_desired_chrom_names = bool(set(categorized_contig_names.autosomes).intersection(DESIRED_AUTOSOME_CONTIG_NAMES))
    has_only_hardmasked_nucleotides_at_y_par1 = not bool(y_test_nucleotides.difference(UNKNOWN_NUCLEOTIDES))
    has_semi_ambiguous_iub_codes = bool(nucleotides.difference(STANDARD_NUCLEOTIDES).difference(SOFTMASKED_NUCLEOTIDES))
    has_softmasked_nucleotides = bool(nucleotides.intersection(SOFTMASKED_NUCLEOTIDES))

    logging.info(f"FEATURES GENOME:")
    logging.info(f"Unplaced contigs: {bool(categorized_contig_names.unplaced_contigs)}")
    logging.info(f"Unlocalized contigs: {bool(categorized_contig_names.unlocalized_contigs)}")
    logging.info(f"ALTS: {bool(categorized_contig_names.alt_contigs)}")
    logging.info(f"rCRS mitochondrial sequence: {has_rcrs}")
    logging.info(f"Accession numbers: ?")
    logging.info(f"Uses 'chr1' chrom names: {uses_desired_chrom_names}")  # TODO: do this properly !!
    logging.info(f"PAR hardmask (not fully accurate): {has_only_hardmasked_nucleotides_at_y_par1}")
    logging.info(f"Decoys (hs38d1): {bool(categorized_contig_names.decoys)}")
    logging.info(f"EBV: {categorized_contig_names.ebv is not None}")
    logging.info(f"Patches: {bool(categorized_contig_names.fix_patch_contigs) or bool(categorized_contig_names.novel_patch_contigs)}")
    logging.info(f"PhiX: ?")
    logging.info(f"Semi ambiguous IUB codes: {has_semi_ambiguous_iub_codes}")
    logging.info(f"Has softmasked nucleotides: {has_softmasked_nucleotides}")


def get_contig_name_translator(contig_alias_bucket_path: str) -> ContigNameTranslator:
    contig_alias_text = get_blob(contig_alias_bucket_path).download_as_text()

    contig_name_to_canonical_name = {}
    for line in contig_alias_text.split("\n"):
        split_line = line.split("\t")
        if len(split_line) != 2:
            raise ValueError(f"Incorrect length line: {line}")
        contig_name, canonical_name = split_line
        if contig_name in contig_name_to_canonical_name:
            raise ValueError(f"Encountered contig name multiple times: {contig_name}")
        contig_name_to_canonical_name[contig_name] = canonical_name
    return ContigNameTranslator(contig_name_to_canonical_name)


def get_categorized_contig_names(
    ref_genome_path: Path,
    contig_name_translator: ContigNameTranslator,
) -> CategorizedContigNames:
    with pysam.Fastafile(ref_genome_path) as genome_f:
        contig_names = list(genome_f.references)

    autosome_contigs: List[str] = []
    x_contigs: List[str] = []
    y_contigs: List[str] = []
    mitochondrial_contigs: List[str] = []
    ebv_contigs: List[str] = []
    decoy_contigs: List[str] = []
    unlocalized_contigs: List[str] = []
    unplaced_contigs: List[str] = []
    alt_contigs: List[str] = []
    fix_patch_contigs: List[str] = []
    novel_patch_contigs: List[str] = []
    uncategorized_contigs: List[str] = []

    for contig_name in contig_names:
        standardized_contig_name = contig_name_translator.standardize(contig_name)
        if standardized_contig_name in DESIRED_AUTOSOME_CONTIG_NAMES:
            autosome_contigs.append(contig_name)
        elif standardized_contig_name in X_CHROMOSOME_CONTIG_NAMES:
            x_contigs.append(contig_name)
        elif standardized_contig_name in Y_CHROMOSOME_CONTIG_NAMES:
            y_contigs.append(contig_name)
        elif standardized_contig_name in MITOCHONDRIAL_CONTIG_NAMES:
            mitochondrial_contigs.append(contig_name)
        elif standardized_contig_name in EBV_CONTIG_NAMES:
            ebv_contigs.append(contig_name)
        elif is_decoy_contig_name(standardized_contig_name):
            decoy_contigs.append(contig_name)
        elif is_unlocalized_contig_name(standardized_contig_name):
            unlocalized_contigs.append(contig_name)
        elif is_unplaced_contig_name(standardized_contig_name):
            unplaced_contigs.append(contig_name)
        elif is_alt_contig_name(standardized_contig_name):
            alt_contigs.append(contig_name)
        elif is_fix_patch_contig_name(standardized_contig_name):
            fix_patch_contigs.append(contig_name)
        elif is_novel_patch_contig_name(standardized_contig_name):
            novel_patch_contigs.append(contig_name)
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
        tuple(alt_contigs),
        tuple(fix_patch_contigs),
        tuple(novel_patch_contigs),
        tuple(uncategorized_contigs),
    )
    return categorized_contigs


def is_decoy_contig_name(contig_name: str) -> bool:
    return contig_name.startswith(UNPLACED_PREFIX) and contig_name.endswith(DECOY_SUFFIX)


def is_unlocalized_contig_name(contig_name: str) -> bool:
    return not contig_name.startswith(UNPLACED_PREFIX) and contig_name.endswith(UNLOCALIZED_SUFFIX)


def is_unplaced_contig_name(contig_name: str) -> bool:
    return contig_name.startswith(UNPLACED_PREFIX) and not contig_name.endswith(DECOY_SUFFIX)


def is_alt_contig_name(contig_name: str) -> bool:
    return not contig_name.startswith(UNPLACED_PREFIX) and contig_name.endswith(ALT_SUFFIX)


def is_fix_patch_contig_name(contig_name: str) -> bool:
    return not contig_name.startswith(UNPLACED_PREFIX) and contig_name.endswith(FIX_PATCH_SUFFIX)


def is_novel_patch_contig_name(contig_name: str) -> bool:
    return not contig_name.startswith(UNPLACED_PREFIX) and contig_name.endswith(NOVEL_PATCH_SUFFIX)


def mitochondrial_sequence_is_rcrs(config: Config, ref_mitochondrial_contig_name: str) -> bool:
    with pysam.Fastafile(config.rcrs_path) as rcrs_f:
        if rcrs_f.nreferences != 1:
            raise ValueError(f"rCRS FASTA file has more than one contig")
        rcrs_genome = rcrs_f.fetch(rcrs_f.references[0])
    with pysam.Fastafile(config.ref_genome_path) as genome_f:
        mitochondrial_from_ref = genome_f.fetch(ref_mitochondrial_contig_name)
    return rcrs_genome == mitochondrial_from_ref


def get_y_test_sequence(ref_genome_path: Path, y_contig_name: str) -> str:
    with pysam.Fastafile(ref_genome_path) as genome_f:
        return genome_f.fetch(y_contig_name, Y_PAR1_TEST_REGION[0], Y_PAR1_TEST_REGION[1])


def get_nucleotides(ref_genome_path: Path) -> Set[str]:
    futures = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        with pysam.Fastafile(ref_genome_path) as genome_f:
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


def assert_file_exists_in_bucket(path: str) -> None:
    if not get_blob(path).exists():
        raise ValueError(f"File in bucket does not exist: {path}")


def get_blob(path: str) -> storage.Blob:
    bucket_name = path.split("/")[2]
    relative_path = "/".join(path.split("/")[3:])
    return storage.Client().get_bucket(bucket_name).get_blob(relative_path)


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="check_hg38_ref_genome_features",
        description=(
            "Check important features for ref genome FASTA file. Correctness is not guaranteed, especially for hg19. "
            "See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files."
        ),
    )
    parser.add_argument("--fasta", "-i", type=Path, required=True, help="Fasta file for ref genome.")
    parser.add_argument("--rcrs", "-r", type=Path, required=True, help="Fasta file for rCRS mitochondrial genome.")
    parser.add_argument(
        "--contig_alias",
        "-c",
        type=str,
        required=True,
        help=(
            "Bucket path to TSV file with contig name translations. Source: create_contig_translation_file.py"
        )
    )

    args = parser.parse_args(sys_args)

    config = Config(args.fasta, args.rcrs, args.contig_alias)
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
