from typing import Tuple, NamedTuple

from genome import BafSite, FusionSite, Position, MsiSite, PgxSite, Interval, Exon


class PanelFileConfig(NamedTuple):
    all_genes_tsv: str
    baf_sites_list: str
    fusion_sites_list: str
    gene_list: str
    hotspot_list: str
    msi_sites_list: str
    pgx_sites_list: str
    tert_site: str


class AnalysisTodoConfig(NamedTuple):
    baf: bool
    exome: bool
    fusion: bool
    hotspot: bool
    msi: bool
    pgx: bool
    tert: bool


class Config(NamedTuple):
    sample_with_depth_file_pairs: Tuple[Tuple[str, str]]
    panel_file_config: PanelFileConfig
    analysis_todo_config: AnalysisTodoConfig
    min_coverage: int
    output_dir: str


class Panel(NamedTuple):
    baf_sites: Tuple[BafSite]
    exons: Tuple[Exon]
    fusion_sites: Tuple[FusionSite]
    hotspots: Tuple[Position]
    msi_sites: Tuple[MsiSite]
    pgx_sites: Tuple[PgxSite]
    tert_site: Interval
