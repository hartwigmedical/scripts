from copy import deepcopy
from pathlib import Path
from typing import NamedTuple, Tuple

import pandas as pd


class AllGenesData(object):
    CHROM = "chrom"
    GENE_START = "gene_start"
    GENE_END = "gene_end"
    GENE_ENSEMBL_ID = "gene_id"
    GENE_NAME = "gene"
    TRANSCRIPT_ENSEMBL_ID = "transcript_id"
    GENE_START2 = "gene_start2"
    GENE_END2 = "gene_end2"
    EXON_ENSEMBL_ID = "exon_id"
    EXON_START_WITH_UTR = "exon_start_with_utr"
    EXON_END_WITH_UTR = "exon_end_with_utr"
    TRANSCRIPT_START = "transcript_start"
    TRANSCRIPT_END = "transcript_end"

    COLUMN_NAMES = [
        CHROM, GENE_START, GENE_END, GENE_ENSEMBL_ID, GENE_NAME, TRANSCRIPT_ENSEMBL_ID, GENE_START2, GENE_END2,
        EXON_ENSEMBL_ID, EXON_START_WITH_UTR, EXON_END_WITH_UTR, TRANSCRIPT_START, TRANSCRIPT_END,
    ]
    USED_COLS = [0, 1, 2, 3, 4, 7, 9, 10, 11, 12, 13, 15, 16]
    COLUMN_TO_TYPE = {
        CHROM: str,
        GENE_START: int,
        GENE_END: int,
        GENE_ENSEMBL_ID: str,
        GENE_NAME: str,
        TRANSCRIPT_ENSEMBL_ID: str,
        GENE_START2: int,
        GENE_END2: int,
        EXON_ENSEMBL_ID: str,
        EXON_START_WITH_UTR: int,
        EXON_END_WITH_UTR: int,
        TRANSCRIPT_START: pd.Int64Dtype(),
        TRANSCRIPT_END: pd.Int64Dtype(),
    }

    def __init__(self, df: pd.DataFrame) -> None:
        self.df = deepcopy(df)

    @classmethod
    def from_file(cls, path: Path) -> "AllGenesData":
        df = pd.read_csv(
            path, sep="\t", names=cls.COLUMN_NAMES, index_col=False, usecols=cls.USED_COLS, dtype=cls.COLUMN_TO_TYPE)
        return AllGenesData(df)


class Interval(NamedTuple):
    chromosome: str
    start_position: int
    end_position: int


class Exon(NamedTuple):
    gene: str
    gene_ensembl_id: str
    exon_ensembl_id: str
    interval: Interval


class Exome(object):
    """Only includes transcripts"""

    CHROM = "chrom"
    GENE_NAME = "gene"
    GENE_ENSEMBL_ID = "gene_id"
    EXON_ENSEMBL_ID = "exon_id"
    EXON_START = "exon_start"
    EXON_END = "exon_end"

    COLUMN_NAMES = [
        CHROM, GENE_NAME, GENE_ENSEMBL_ID, EXON_ENSEMBL_ID, EXON_START, EXON_END,
    ]

    def __init__(self, df: pd.DataFrame) -> None:
        self.df = deepcopy(df)

    @classmethod
    def from_all_genes_data(cls, all_genes: AllGenesData) -> "Exome":
        all_genes_with_transcripts_df = all_genes.df.dropna(axis=0)
        exon_start = all_genes_with_transcripts_df[[
            AllGenesData.EXON_START_WITH_UTR,
            AllGenesData.GENE_START,
            AllGenesData.GENE_START2,
            AllGenesData.TRANSCRIPT_START,
        ]].astype(int).max(axis=1)
        exon_end = all_genes_with_transcripts_df[[
            AllGenesData.EXON_END_WITH_UTR,
            AllGenesData.GENE_END,
            AllGenesData.GENE_END2,
            AllGenesData.TRANSCRIPT_END,
        ]].astype(int).min(axis=1)
        exome_df = pd.concat([
            all_genes_with_transcripts_df[AllGenesData.CHROM],
            all_genes_with_transcripts_df[AllGenesData.GENE_NAME],
            all_genes_with_transcripts_df[AllGenesData.GENE_ENSEMBL_ID],
            all_genes_with_transcripts_df[AllGenesData.EXON_ENSEMBL_ID],
            exon_start,
            exon_end,
        ], axis=1)
        exome_df = exome_df[exon_start <= exon_end]
        exome_df.columns = cls.COLUMN_NAMES
        exome_df = exome_df.reset_index(drop=True)
        return Exome(exome_df)

    def get_exons(self, gene_list: Tuple[str, ...]) -> Tuple[Exon, ...]:
        relevant_df = self.df[self.df[self.GENE_NAME].isin(gene_list)]
        result = {
            Exon(gene, gene_id, exon_id, Interval(chrom, int(exon_start), int(exon_end)))
            for index, chrom, gene, gene_id, exon_id, exon_start, exon_end in
            relevant_df.itertuples()
        }
        return tuple(result)


class Position(NamedTuple):
    chromosome: str
    position: int

    def get_interval(self) -> Interval:
        return Interval(self.chromosome, self.position, self.position)


class MsiSite(NamedTuple):
    chromosome: str
    position: int
    repeat_count: int
    probe3_start: int
    probe3_end: int
    probe3_id: str
    probe5_start: int
    probe5_end: int
    probe5_id: str

    def get_site_interval(self) -> Interval:
        return Interval(self.chromosome, self.position, self.position + self.repeat_count)


class BafSite(NamedTuple):
    chromosome: str
    position: int
    label: str
    probe_start: int
    probe_end: int

    def get_site_interval(self) -> Interval:
        return Interval(self.chromosome, self.position, self.position)


class PgxSite(NamedTuple):
    interval: Interval
    gene: str
    label: str


class FusionSite(NamedTuple):
    interval: Interval
    gene: str
    intron_start: int
    intron_end: int
