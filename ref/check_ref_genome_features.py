import argparse
import concurrent.futures
import logging
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import List, NamedTuple, Tuple, Set, Dict

from google.cloud import storage
import pysam

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

AUTOSOME_CONTIG_NAMES = {f"chr{n}" for n in range(1, 23)}
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
    x_contigs: Tuple[str]
    y_contigs: Tuple[str]
    mitochondrial_contigs: Tuple[str]
    ebv_contigs: Tuple[str]
    decoys: Tuple[str]
    unlocalized_contigs: Tuple[str]
    unplaced_contigs: Tuple[str]
    alt_contigs: Tuple[str]
    fix_patch_contigs: Tuple[str]
    novel_patch_contigs: Tuple[str]
    uncategorized_contigs: Tuple[str]

    def get_contig_names(self) -> Tuple[str]:
        contig_names = []
        contig_names.extend(self.autosomes)
        contig_names.extend(self.x_contigs)
        contig_names.extend(self.y_contigs)
        contig_names.extend(self.mitochondrial_contigs)
        contig_names.extend(self.ebv_contigs)
        contig_names.extend(self.decoys)
        contig_names.extend(self.unlocalized_contigs)
        contig_names.extend(self.unplaced_contigs)
        contig_names.extend(self.alt_contigs)
        contig_names.extend(self.fix_patch_contigs)
        contig_names.extend(self.novel_patch_contigs)
        contig_names.extend(self.uncategorized_contigs)
        return tuple(contig_names)


class ContigNameTranslator(object):
    """Standardizes names if it can. Returns argument as is if it cannot."""
    def __init__(self, contig_name_to_canonical_name: Dict[str, str]) -> None:
        self._contig_name_to_canonical_name = deepcopy(contig_name_to_canonical_name)

    def standardize(self, contig_name: str) -> str:
        if contig_name in self._contig_name_to_canonical_name:
            return self._contig_name_to_canonical_name[contig_name]
        else:
            raise ValueError(f"Could not standardize '{contig_name}'")

    def is_canonical(self, contig_name: str) -> bool:
        return contig_name in self._contig_name_to_canonical_name.values()


def main(config: Config) -> None:
    set_up_logging()
    config.validate()

    contig_name_translator = get_contig_name_translator(config.contig_alias_bucket_path)

    categorized_contig_names = get_categorized_contig_names(config.ref_genome_path, contig_name_translator)
    logging.debug(categorized_contig_names)

    nucleotides_file = Path(f"{config.ref_genome_path}.nucleotides")
    if not nucleotides_file.exists():
        nucleotides = get_nucleotides(config.ref_genome_path)
        nucleotides_string = "".join(sorted(list(nucleotides)))
        with open(nucleotides_file, "w+") as f:
            f.write(nucleotides_string)

    with open(nucleotides_file, "r+") as f:
        line = f.readline()
        nucleotides = {char for char in line}

    logging.info(f"nucleotides: {sorted(nucleotides)}")

    if len(categorized_contig_names.autosomes) != 22:
        warn_msg = (
            f"Did not find exactly 22 autosome contigs: "
            f"autosomes={categorized_contig_names.autosomes}"
        )
        logging.warning(warn_msg)
    if len(categorized_contig_names.x_contigs) != 1:
        warn_msg = (
            f"Did not find exactly one X contig: "
            f"x={categorized_contig_names.x_contigs}"
        )
        logging.warning(warn_msg)

    has_unplaced_contigs = bool(categorized_contig_names.unplaced_contigs)
    has_unlocalized_contigs = bool(categorized_contig_names.unlocalized_contigs)
    has_alts = bool(categorized_contig_names.alt_contigs)
    has_decoys = bool(categorized_contig_names.decoys)
    has_patches = bool(categorized_contig_names.fix_patch_contigs) or bool(categorized_contig_names.novel_patch_contigs)
    if len(categorized_contig_names.ebv_contigs) == 1:
        has_ebv = True
    elif len(categorized_contig_names.ebv_contigs) == 0:
        has_ebv = False
    else:
        logging.warning(f"Found more than one EBV contig: {categorized_contig_names.ebv_contigs}")
        has_ebv = None
    if len(categorized_contig_names.mitochondrial_contigs) == 1:
        has_rcrs = mitochondrial_sequence_is_rcrs(config, categorized_contig_names.mitochondrial_contigs[0])
    else:
        warn_msg = (
            f"Did not find exactly one mitochondrial contig: "
            f"mitochondrial={categorized_contig_names.mitochondrial_contigs}"
        )
        logging.warning(warn_msg)
        has_rcrs = None
    uses_canonical_chrom_names = all(
        contig_name_translator.is_canonical(contig_name)
        for contig_name in categorized_contig_names.get_contig_names()
    )
    if len(categorized_contig_names.y_contigs) == 1:
        y_test_nucleotides = get_nucleotides_from_string(
            get_y_test_sequence(categorized_contig_names.y_contigs[0], config.ref_genome_path)
        )
        logging.info(f"nucleotides at y par1 test region: {sorted(y_test_nucleotides)}")
        has_only_hardmasked_nucleotides_at_y_par1 = not bool(y_test_nucleotides.difference(UNKNOWN_NUCLEOTIDES))
    else:
        warn_msg = (
            f"Did not find exactly one Y contig: "
            f"y={categorized_contig_names.y_contigs}"
        )
        logging.warning(warn_msg)
        has_only_hardmasked_nucleotides_at_y_par1 = None
    has_semi_ambiguous_iub_codes = bool(nucleotides.difference(STANDARD_NUCLEOTIDES).difference(SOFTMASKED_NUCLEOTIDES))
    has_softmasked_nucleotides = bool(nucleotides.intersection(SOFTMASKED_NUCLEOTIDES))
    if all(is_definitely_padded_with_n(contig, config.ref_genome_path) for contig in categorized_contig_names.alt_contigs):
        alts_are_padded = True
    elif all(not is_definitely_padded_with_n(contig, config.ref_genome_path) for contig in categorized_contig_names.alt_contigs):
        alts_are_padded = False
    else:
        logging.warning(f"Could not determine whether alts are padded with N's (or n's)")
        alts_are_padded = None

    logging.info(f"FEATURES GENOME:")

    logging.info(f"Unplaced contigs: {has_unplaced_contigs}")
    logging.info(f"Unlocalized contigs: {has_unlocalized_contigs}")
    logging.info(f"Alts: {has_alts}")
    logging.info(f"Decoys (hs38d1): {has_decoys}")
    logging.info(f"Patches: {has_patches}")
    logging.info(f"EBV: {has_ebv}")
    logging.info(f"rCRS mitochondrial sequence: {has_rcrs}")
    logging.info(f"Uses canonical contig names, so 'chr1' etc.: {uses_canonical_chrom_names}")
    logging.info(f"PAR hardmask (not fully accurate): {has_only_hardmasked_nucleotides_at_y_par1}")
    logging.info(f"Semi ambiguous IUB codes: {has_semi_ambiguous_iub_codes}")
    logging.info(f"Has softmasked nucleotides: {has_softmasked_nucleotides}")
    logging.info(f"Alts are padded with N: {alts_are_padded}")
    logging.info(f"PhiX: False?")
    logging.info(f"")
    logging.info(f"For easy copy-paste:")
    answers: List[bool] = [
        has_unplaced_contigs,
        has_unlocalized_contigs,
        has_alts,
        has_decoys,
        has_patches,
        has_ebv,
        has_rcrs,
        uses_canonical_chrom_names,
        has_only_hardmasked_nucleotides_at_y_par1,
        has_semi_ambiguous_iub_codes,
        has_softmasked_nucleotides,
        alts_are_padded,
    ]
    value_to_answer = {True: "Yes", False: "No", None: "?"}
    print("\n".join([value_to_answer[answer] for answer in answers]))


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
        if standardized_contig_name in AUTOSOME_CONTIG_NAMES:
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

    categorized_contigs = CategorizedContigNames(
        tuple(autosome_contigs),
        tuple(x_contigs),
        tuple(y_contigs),
        tuple(mitochondrial_contigs),
        tuple(ebv_contigs),
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


def get_y_test_sequence(y_contig_name: str, ref_genome_path: Path) -> str:
    with pysam.Fastafile(ref_genome_path) as genome_f:
        return genome_f.fetch(y_contig_name, Y_PAR1_TEST_REGION[0], Y_PAR1_TEST_REGION[1])


def is_definitely_padded_with_n(contig_name: str, ref_genome_path: Path) -> bool:
    with pysam.Fastafile(ref_genome_path) as genome_f:
        contig_sequence = genome_f.fetch(contig_name)

    first_1000_nucleotides = get_nucleotides(contig_sequence[:1000])
    last_1000_nucleotides = get_nucleotides(contig_sequence[-1000:])

    return first_1000_nucleotides.issubset(UNKNOWN_NUCLEOTIDES) and last_1000_nucleotides.issubset(UNKNOWN_NUCLEOTIDES)


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
