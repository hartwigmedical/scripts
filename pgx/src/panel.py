from collections import defaultdict
from typing import NamedTuple, Dict, Any, List, Optional, Set


class SNP(NamedTuple):
    rs_id: str
    chromosome: int
    position: int

    @classmethod
    def from_json(cls, data: Dict[str, str]) -> "SNP":
        rs_id = data['rsid']
        chromosome = int(data['chromosome'])
        position = int(data['position'])
        return SNP(rs_id, chromosome, position)

    def get_position_string(self):
        return f"{self.chromosome}:{self.position}"

    def matches_position_string(self, position_string: str) -> bool:
        return self.get_position_string() == position_string


class GeneInfo(NamedTuple):
    dict: Dict[str, Any]


class Panel(NamedTuple):
    haplotypes_info: Dict[str, GeneInfo]  # Nested dict: gene, haplotype, rsinfo
    snps: Set[SNP]

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Panel":
        haplotypes_info = {}
        snps = set()
        for gene_info_json in data['genes']:
            gene_info = GeneInfo(gene_info_json)
            if gene_info.dict['genomeBuild'] != "GRCh37":
                raise ValueError("Exiting, we only support GRCh37, not " + str(gene_info.dict['genomeBuild']))
            if gene_info.dict['gene'] not in haplotypes_info:
                haplotypes_info[gene_info.dict['gene']] = gene_info
            for variant in gene_info_json['variants']:
                snps.add(SNP.from_json(variant))

        if len({snp.rs_id for snp in snps}) != len(snps):
            # Error: some snps share rs_id but are different
            rs_id_to_snps = defaultdict(list)
            for snp in snps:
                rs_id_to_snps[snp.rs_id].append(snp)

            rs_id_to_duplicate_snps = {rs_id: snps for rs_id, snps in rs_id_to_snps.items() if len(snps) > 1}

            raise ValueError(
                ("Panel json contains snps with the same rs id but different positions. "
                 "Duplicates: {rs_id_to_duplicate_snps}").format(
                    rs_id_to_duplicate_snps=rs_id_to_duplicate_snps
                )
            )

        return Panel(haplotypes_info, snps)

    def contains_snp_with_position(self, position_string: str) -> bool:
        for snp in self.snps:
            if snp.matches_position_string(position_string):
                return True
        return False

    def get_snp_with_position(self, position_string: str) -> SNP:
        matching_snps = []
        for snp in self.snps:
            if snp.matches_position_string(position_string):
                matching_snps.append(snp)

        if matching_snps and len(matching_snps) == 1:
            return matching_snps.pop()
        elif not matching_snps:
            raise ValueError("No snps match position")
        else:
            raise ValueError("Multiple snps match position")

    def contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.get_rs_ids()

    def get_rs_ids(self) -> Set[str]:
        return {snp.rs_id for snp in self.snps}

    def is_empty(self) -> bool:
        return not self.haplotypes_info and not self.snps

    def get_genes(self) -> List[str]:
        return list(self.haplotypes_info.keys())



