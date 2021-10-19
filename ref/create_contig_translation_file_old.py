import argparse
import logging
import re
import sys
from collections import defaultdict
from concurrent import futures
from pathlib import Path
from typing import List, NamedTuple, Tuple, Dict, DefaultDict, Set

from ref_util import assert_file_exists, set_up_logging

SCRIPT_NAME = "create_chrom_translation_file_old"


class Config(NamedTuple):
    fasta_paths: Tuple[Path]

    def validate(self) -> None:
        for fasta_path in self.fasta_paths:
            assert_file_exists(fasta_path)


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    logging.info("Start reading FASTA files")
    contig_name_to_info_list: List[Dict[str, str]] = []
    with futures.ProcessPoolExecutor() as executor:
        future_to_fasta_path = {
            executor.submit(get_contig_name_to_info, fasta_path): fasta_path for fasta_path in config.fasta_paths
        }
        for future in futures.as_completed(future_to_fasta_path.keys()):
            fasta_path = future_to_fasta_path[future]
            try:
                contig_name_to_info = future.result()
                contig_name_to_info_list.append(contig_name_to_info)
                logging.info(f"Finished reading {fasta_path}")
            except Exception as exc:
                logging.warning(f"Could not process {fasta_path}: {exc}")

    logging.info("Finished reading FASTA files")

    seen_contig_names: Set[str] = set()
    canonical_contig_name_to_aliases: DefaultDict[str, Set[str]] = defaultdict(set)

    add_hardcoded_aliases(canonical_contig_name_to_aliases)

    for contig_name_to_info in contig_name_to_info_list:
        for contig_name, info in contig_name_to_info.items():
            seen_contig_names.add(contig_name)
            if "AC:" in info:
                # From GCA/GCF analysis set
                alias_matches = re.findall(r"AC:[\S]+", info)
                if len(alias_matches) != 1:
                    raise SyntaxError(f"Incomplete match with 'AC:...' substring: {info}")
                canonical_contig_name_to_aliases[contig_name].add(contig_name)
                canonical_contig_name_to_aliases[contig_name].add(alias_matches[0][3:])
            elif re.fullmatch(r"[A-Z]{2}[0-9]{6}\.[0-9]", contig_name) or re.fullmatch(r"NC_[0-9]{6}\.[0-9]", contig_name):
                contig_name_without_dot = contig_name.replace(".", "v")
                if "GRCh38 reference primary assembly" in info:
                    # From GCA/GCF non-analysis set
                    if "unplaced genomic contig" in info:
                        canonical_contig_name_to_aliases[f"chrUn_{contig_name_without_dot}"].add(contig_name)
                    elif "unlocalized genomic contig" in info:
                        regex_matches = re.findall(r"Homo sapiens chromosome [0-9A-Z]{1,2} unlocalized genomic contig", info)
                        if len(regex_matches) != 1:
                            raise ValueError(f"Could not find chromosome for unlocalized genomic contig: ({contig_name}, {info})")
                        chromosome = regex_matches[0].split(" ")[3]
                        canonical_contig_name_to_aliases[f"chr{chromosome}_{contig_name_without_dot}_random"].add(contig_name)
                    else:
                        regex_matches = re.findall(r"Homo sapiens chromosome [0-9A-Z]{1,2},", info)
                        if len(regex_matches) != 1:
                            raise ValueError(f"Could not find chromosome for normal genomic contig: ({contig_name}, {info})")
                        chromosome = regex_matches[0].strip(",").split(" ")[3]
                        canonical_contig_name_to_aliases[f"chr{chromosome}"].add(contig_name)
                elif "PATCH for GRCh38" in info:
                    # From GCA/GCF non-analysis set
                    regex_matches = re.findall(r"Homo sapiens chromosome [0-9A-Z]{1,2} genomic contig", info)
                    if len(regex_matches) != 1:
                        raise ValueError(f"Could not find chromosome for patch genomic contig: ({contig_name}, {info})")
                    chromosome = regex_matches[0].split(" ")[3]
                    if "FIX PATCH for GRCh38" in info:
                        canonical_contig_name_to_aliases[f"chr{chromosome}_{contig_name_without_dot}_fix"].add(contig_name)
                    elif "NOVEL PATCH for GRCh38" in info:
                        canonical_contig_name_to_aliases[f"chr{chromosome}_{contig_name_without_dot}_novel"].add(contig_name)
                    else:
                        raise ValueError(f"Found patch that is not fix or novel: ({contig_name}, {info})")
                elif "GRCh38 reference assembly alternate locus group ALT_REF_LOCI" in info:
                    # From GCA/GCF non-analysis set
                    regex_matches = re.findall(r"Homo sapiens chromosome [0-9A-Z]{1,2} genomic contig", info)
                    if len(regex_matches) != 1:
                        raise ValueError(f"Could not find chromosome for patch genomic contig: ({contig_name}, {info})")
                    chromosome = regex_matches[0].split(" ")[3]
                    canonical_contig_name_to_aliases[f"chr{chromosome}_{contig_name_without_dot}_alt"].add(contig_name)
                elif "mitochondrion" in info:
                    # From GCA/GCF non-analysis set
                    canonical_contig_name_to_aliases[f"chrM"].add(contig_name)
                elif "Homo sapiens unplaced genomic scaffold decoy" in info:
                    # From decoy file
                    canonical_contig_name_to_aliases[f"chrUn_{contig_name_without_dot}_decoy"].add(contig_name)

    logging.info(canonical_contig_name_to_aliases)

    handled_contig_names = set.union(set(), *canonical_contig_name_to_aliases.values())
    unsolved_contig_names = seen_contig_names.difference(handled_contig_names)
    logging.info(f"Unsolved contig names: {sorted(unsolved_contig_names)}")

    total_categorized_contigs_count_separate = sum(len(aliases) for aliases in canonical_contig_name_to_aliases.values())
    total_categorized_contigs_count_combined = len(set.union(*canonical_contig_name_to_aliases.values()))
    if total_categorized_contigs_count_separate != total_categorized_contigs_count_combined:
        raise ValueError(
            f"Counts don't match. Some aliases are present for multiple canonical contig names: "
            f"{total_categorized_contigs_count_separate}, {total_categorized_contigs_count_combined}"
        )


def get_contig_name_to_info(fasta_path: Path) -> Dict[str, str]:
    contig_name_to_info: Dict[str, str] = {}
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip("\n")
            line_is_contig_header = line.startswith(">")
            if line_is_contig_header:
                contig_name = line.split(" ")[0][1:]
                info = line[len(contig_name)+1:]  # strip out contig name
                if contig_name in contig_name_to_info:
                    raise SyntaxError(f"Fasta file {fasta_path} contains contig {contig_name} more than once.")
                contig_name_to_info[contig_name] = info
    return contig_name_to_info


def add_hardcoded_aliases(canonical_contig_name_to_aliases: DefaultDict[str, Set[str]]) -> None:
    for i in range(1, 23):
        canonical_contig_name_to_aliases[f"chr{i}"].add(str(i))
    canonical_contig_name_to_aliases["chrX"].add("X")
    canonical_contig_name_to_aliases["chrY"].add("Y")
    canonical_contig_name_to_aliases["chrM"].add("MT")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description=(
            "Create translation file for ref genome fasta file comparison."
        ),
    )
    parser.add_argument(
        "--fasta",
        "-i",
        type=Path,
        required=True,
        action="append",
        help="Fasta file for ref genome. Option can be used multiple times.",
    )

    args = parser.parse_args(sys_args)

    config = Config(tuple(args.fasta))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
