#!/usr/bin/env python3
import logging
from pathlib import Path

GENE_INPUT_FILE = Path("/data/experiments/230620_erik_INN-104_AltTranscript/ensembl_data_cache/37/ensembl_gene_data.csv")
TRANSCRIPT_INPUT_FILE = Path("/data/experiments/230620_erik_INN-104_AltTranscript/ensembl_data_cache/37/ensembl_trans_exon_data.csv")

def main() -> None:
    logging.info("Start comparing gene and transcript/coding start and end")

    gene_id_to_gene_start_end = {}
    with open(GENE_INPUT_FILE, "r") as in_f:
        # skip header
        next(in_f)

        for line in in_f:
            if line:
                split_line = line.strip().split(",")
                gene_id = split_line[0]
                gene_start = int(split_line[4])
                gene_end = int(split_line[5])
                if gene_id in gene_id_to_gene_start_end.keys():
                    raise ValueError(f"Already encountered gene id: {gene_id}")
                gene_id_to_gene_start_end[gene_id] = (gene_start, gene_end)

    warnings = []
    warning_set = set()

    with open(TRANSCRIPT_INPUT_FILE, "r") as in_f:
        # skip header
        next(in_f)
        for line in in_f:
            if line:
                split_line = line.strip().split(",")
                gene_id = split_line[0]
                transcript_id = split_line[4]
                transcript_start = int(split_line[6]) if split_line[6] != "NULL" else None
                transcript_end = int(split_line[7]) if split_line[7] != "NULL" else None
                exon_start = int(split_line[9]) if split_line[9] != "NULL" else None
                exon_end = int(split_line[10]) if split_line[10] != "NULL" else None
                coding_start = int(split_line[13]) if split_line[13] != "NULL" else None
                coding_end = int(split_line[14]) if split_line[14] != "NULL" else None

                gene_start, gene_end = gene_id_to_gene_start_end[gene_id]
                # if transcript_start is not None and transcript_start < gene_start:
                #     warning = f"Early transcript start!: {gene_id}:{transcript_id}, {transcript_start} < {gene_start}"
                #     if warning not in warning_set:
                #         warnings.append(warning)
                #         warning_set.add(warning)
                # if transcript_end is not None and transcript_end > gene_end:
                #     warning = f"Late transcript end!: {gene_id}:{transcript_id}, {transcript_end} > {gene_end}"
                #     if warning not in warning_set:
                #         warnings.append(warning)
                #         warning_set.add(warning)
                # if exon_start is not None and exon_start < gene_start:
                #     warning = f"Early exon start!: {gene_id}:{transcript_id}, {exon_start} < {gene_start}"
                #     if warning not in warning_set:
                #         warnings.append(warning)
                #         warning_set.add(warning)
                # if exon_end is not None and exon_end > gene_end:
                #     warning = f"Late exon end!: {gene_id}:{transcript_id}, {exon_end} > {gene_end}"
                #     if warning not in warning_set:
                #         warnings.append(warning)
                #         warning_set.add(warning)
                if coding_start is not None and coding_start < gene_start:
                    warning = f"Early coding start!: {gene_id}:{transcript_id}, {coding_start} < {gene_start}"
                    if warning not in warning_set:
                        warnings.append(warning)
                        warning_set.add(warning)
                if coding_end is not None and coding_end > gene_end:
                    warning = f"Late coding end!: {gene_id}:{transcript_id}, {coding_end} > {gene_end}"
                    if warning not in warning_set:
                        warnings.append(warning)
                        warning_set.add(warning)

    for warning in warnings:
        logging.warning(warning)

    logging.info("Finished comparing gene and transcript/coding start and end")


if __name__ == '__main__':
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    main()

