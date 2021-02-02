from typing import NamedTuple, Dict, Any, List, Optional


class SNP(NamedTuple):
    rs_id: str
    chromosome: int
    position: int

    def get_position_string(self):
        return f"{self.chromosome}:{self.position}"

    def matches_position_string(self, position_string: str) -> bool:
        return self.get_position_string() == position_string


class Panel(NamedTuple):
    haplotypes_info: Dict[str, Dict[str, Any]]  # Nested dict: gene, haplotype, rsinfo
    rs_id_to_snp: Dict[str, SNP]  # Position is expressed as {chromosome number}:{base count}

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Panel":
        haplotypes_info = {}
        rs_id_to_snp = {}
        for gene_info in data['genes']:
            if gene_info['genomeBuild'] != "GRCh37":
                raise ValueError("Exiting, we only support GRCh37, not " + str(data['orientation']))
            if gene_info['gene'] not in haplotypes_info:
                haplotypes_info[gene_info['gene']] = gene_info
            for variant in gene_info['variants']:
                rs_id = variant['rsid']
                if rs_id not in rs_id_to_snp:
                    rs_id_to_snp[rs_id] = SNP(rs_id, variant['chromosome'], variant['position'])
        return Panel(haplotypes_info, rs_id_to_snp)

    def contains_snp_with_position(self, position_string: str) -> bool:
        for snp in self.rs_id_to_snp.values():
            if snp.matches_position_string(position_string):
                return True
        return False

    def get_snp_with_position(self, position_string: str) -> Optional[SNP]:
        matching_snps = []
        for snp in self.rs_id_to_snp.values():
            if snp.matches_position_string(position_string):
                matching_snps.append(snp)

        if matching_snps and len(matching_snps) == 1:
            return matching_snps.pop()
        elif not matching_snps:
            raise ValueError("No snps match position")
        else:
            raise ValueError("Multiple snps match position")

    def is_empty(self) -> bool:
        return not self.haplotypes_info and not self.rs_id_to_snp

    def get_genes(self) -> List[str]:
        return list(self.haplotypes_info.keys())



