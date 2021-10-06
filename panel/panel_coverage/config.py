from pathlib import Path
from typing import Tuple, NamedTuple

from genome import BafSite, FusionSite, Position, MsiSite, PgxSite, Interval, Exon


class PanelFileConfig(NamedTuple):
    all_genes_tsv: Path
    baf_sites_list: Path
    fusion_sites_list: Path
    gene_list: Path
    hotspot_list: Path
    msi_sites_list: Path
    pgx_sites_list: Path
    tert_site: Path


class AnalysisTodoConfig(NamedTuple):
    baf: bool
    exome: bool
    fusion: bool
    hotspot: bool
    msi: bool
    pgx: bool
    tert: bool


class Config(NamedTuple):
    sample_with_depth_file_pairs: Tuple[Tuple[str, Path]]
    panel_file_config: PanelFileConfig
    analysis_todo_config: AnalysisTodoConfig
    min_coverage: int
    output_dir: Path


class Panel(NamedTuple):
    baf_sites: Tuple[BafSite]
    exons: Tuple[Exon]
    fusion_sites: Tuple[FusionSite]
    hotspots: Tuple[Position]
    msi_sites: Tuple[MsiSite]
    pgx_sites: Tuple[PgxSite]
    tert_site: Interval
