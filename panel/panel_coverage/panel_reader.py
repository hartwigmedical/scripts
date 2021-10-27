from typing import Tuple

from config import PanelFileConfig, Panel
from genome import AllGenesData, Exome, BafSite, FusionSite, Interval, Position, MsiSite, PgxSite, Exon


class PanelReader(object):
    @classmethod
    def get_panel(cls, config: PanelFileConfig) -> Panel:
        baf_sites = cls.__get_baf_sites(config)
        fusion_sites = cls.__get_fusion_sites(config)
        hotspot = cls.__get_hotspot(config)
        msi_sites = cls.__get_msi_sites(config)
        panel_exons = cls.__get_panel_exons(config)
        pgx_sites = cls.__get_pgx_sites(config)
        tert_site = cls.__get_tert_site(config)
        panel = Panel(baf_sites, panel_exons, fusion_sites, hotspot, msi_sites, pgx_sites, tert_site)
        return panel

    @classmethod
    def __get_baf_sites(cls, config: PanelFileConfig) -> Tuple[BafSite, ...]:
        with open(config.baf_sites_list) as baf_list_f:
            split_lines = [line.split("\t") for line in baf_list_f.read().split("\n") if line != ""]
        baf_sites_list = [
            BafSite(
                f"chr{chromosome}",
                int(position),
                label,
                int(probe3_start),
                int(probe3_end),
            )
            for chromosome, position, label, probe3_start, probe3_end in split_lines
        ]
        return tuple(baf_sites_list)

    @classmethod
    def __get_fusion_sites(cls, config: PanelFileConfig) -> Tuple[FusionSite, ...]:
        with open(config.fusion_sites_list) as fusion_list_f:
            split_lines = [line.split("\t") for line in fusion_list_f.read().split("\n") if line != ""]
        fusion_sites_list = [
            FusionSite(
                Interval(f"chr{chromosome}", int(start_position), int(end_position)),
                gene,
                int(intron_start),
                int(intron_end),
            )
            for
            chromosome, start_position, end_position, gene, intron_start, intron_end in split_lines
        ]
        return tuple(fusion_sites_list)

    @classmethod
    def __get_genes(cls, config: PanelFileConfig) -> Tuple[str, ...]:
        with open(config.gene_list) as gene_list_f:
            gene_list = [gene for gene in gene_list_f.read().split("\n") if gene != ""]
        return tuple(gene_list)

    @classmethod
    def __get_hotspot(cls, config: PanelFileConfig) -> Tuple[Position, ...]:
        with open(config.hotspot_list) as hotspot_list_f:
            split_lines = [line.split("\t") for line in hotspot_list_f.read().split("\n") if line != ""]
        hotspot_list = [Position(f"chr{chromosome}", int(position_str)) for chromosome, position_str in split_lines]
        return tuple(hotspot_list)

    @classmethod
    def __get_msi_sites(cls, config: PanelFileConfig) -> Tuple[MsiSite, ...]:
        with open(config.msi_sites_list) as msi_list_f:
            split_lines = [line.split("\t") for line in msi_list_f.read().split("\n") if line != ""]
        msi_sites_list = [
            MsiSite(
                f"chr{chromosome}",
                int(position),
                int(repeat_length),
                int(probe3_start),
                int(probe3_end),
                probe3_id,
                int(probe5_start),
                int(probe5_end),
                probe5_id
            )
            for
            chromosome, position, repeat_length, probe3_start, probe3_end, probe3_id, probe5_start, probe5_end, probe5_id
            in split_lines
        ]
        return tuple(msi_sites_list)

    @classmethod
    def __get_panel_exons(cls, config: PanelFileConfig) -> Tuple[Exon, ...]:
        all_genes_data = AllGenesData.from_file(config.all_genes_tsv)
        exome = Exome.from_all_genes_data(all_genes_data)
        gene = cls.__get_genes(config)
        return exome.get_exons(gene)

    @classmethod
    def __get_pgx_sites(cls, config: PanelFileConfig) -> Tuple[PgxSite, ...]:
        with open(config.pgx_sites_list) as pgx_list_f:
            split_lines = [line.split("\t") for line in pgx_list_f.read().split("\n") if line != ""]
        pgx_sites_list = [
            PgxSite(
                Interval(f"chr{chromosome}", int(start), int(end)),
                gene,
                label,
            )
            for chromosome, start, end, gene, label in split_lines
        ]
        return tuple(pgx_sites_list)

    @classmethod
    def __get_tert_site(cls, config: PanelFileConfig) -> Interval:
        with open(config.tert_site) as tert_site_f:
            split_lines = [line.split("\t") for line in tert_site_f.read().split("\n") if line != ""]
        tert_sites = [
            Interval(f"chr{chromosome}", int(start), int(end))
            for chromosome, start, end in split_lines
        ]
        if len(tert_sites) != 1:
            error_msg = f"TERT site list file should only contain one interval."
            raise SyntaxError(error_msg)
        return tert_sites[0]
