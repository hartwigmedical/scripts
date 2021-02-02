from typing import NamedTuple, Dict, Any, List


class Panel(NamedTuple):
    haplotypes_info: Dict[str, Dict[str, Any]]  # Nested dict: gene, haplotype, rsinfo
    rs_id_to_position: Dict[str, str]  # Position is expressed as {chromosome number}:{base count}

    @classmethod
    def from_json(cls, data: Dict[str, Any]):
        haplotypes_info = {}
        rs_id_to_position = {}
        for gene_info in data['genes']:
            if gene_info['genomeBuild'] != "GRCh37":
                raise ValueError("Exiting, we only support GRCh37, not " + str(data['orientation']))
            if gene_info['gene'] not in haplotypes_info:
                haplotypes_info[gene_info['gene']] = gene_info
            for variant in gene_info['variants']:
                if variant['rsid'] not in rs_id_to_position:
                    rs_id_to_position[variant['rsid']] = str(variant['chromosome']) + ":" + str(variant['position'])
        return Panel(haplotypes_info, rs_id_to_position)

    def is_empty(self) -> bool:
        return not self.haplotypes_info and not self.rs_id_to_position

    def get_genes(self) -> List[str]:
        return list(self.haplotypes_info.keys())
