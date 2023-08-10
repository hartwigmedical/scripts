from pathlib import Path
from typing import Tuple, NamedTuple

from gcp.base import GCPPath
from genome import BafSite, FusionSite, Position, MsiSite, PgxSite, Interval, Exon
from util import assert_file_exists


class ProgramConfig(NamedTuple):
    panel_config_dir: Path
    output_dir: GCPPath
    samtools: Path
    working_dir: Path
    min_coverages: Tuple[int, ...]
    bams: Tuple[GCPPath, ...]


class PanelFileConfig(NamedTuple):
    all_genes_tsv: Path
    baf_sites_list: Path
    fusion_sites_list: Path
    gene_list: Path
    hotspot_list: Path
    msi_sites_list: Path
    pgx_sites_list: Path
    tert_site: Path

    def validate(self) -> None:
        assert_file_exists(self.all_genes_tsv)
        assert_file_exists(self.baf_sites_list)
        assert_file_exists(self.fusion_sites_list)
        assert_file_exists(self.gene_list)
        assert_file_exists(self.hotspot_list)
        assert_file_exists(self.msi_sites_list)
        assert_file_exists(self.pgx_sites_list)
        assert_file_exists(self.tert_site)


class AnalysisTypeConfig(NamedTuple):
    baf: bool
    exome: bool
    fusion: bool
    hotspot: bool
    msi: bool
    pgx: bool
    tert: bool


class Panel(NamedTuple):
    baf_sites: Tuple[BafSite, ...]
    exons: Tuple[Exon, ...]
    fusion_sites: Tuple[FusionSite, ...]
    hotspots: Tuple[Position, ...]
    msi_sites: Tuple[MsiSite, ...]
    pgx_sites: Tuple[PgxSite, ...]
    tert_site: Interval
