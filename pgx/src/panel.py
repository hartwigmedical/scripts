from typing import NamedTuple, Dict, Any


class Panel(NamedTuple):
    haplotypes_info: Dict[str, Dict[str, Any]]
    rs_id_to_position: Dict[str, str]

    def is_empty(self) -> bool:
        return not self.haplotypes_info and not self.rs_id_to_position
