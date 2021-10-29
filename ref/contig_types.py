from enum import unique, Enum, auto
from typing import NamedTuple, FrozenSet, Set


@unique
class Assembly(Enum):
    GRCH38 = auto()
    HS38D1 = auto()
    OTHER = auto()


@unique
class ContigType(Enum):
    AUTOSOME = auto()
    X = auto()
    Y = auto()
    MITOCHONDRIAL = auto()
    EBV = auto()
    UNLOCALIZED = auto()
    UNPLACED = auto()
    ALT = auto()
    FIX_PATCH = auto()
    NOVEL_PATCH = auto()
    DECOY = auto()
    UNCATEGORIZABLE = auto()

    def get_assembly(self) -> Assembly:
        grch38_contigs = {
            self.AUTOSOME,
            self.X,
            self.Y,
            self.MITOCHONDRIAL,
            self.UNLOCALIZED,
            self.UNPLACED,
            self.ALT,
            self.FIX_PATCH,
            self.NOVEL_PATCH,
        }
        if self in grch38_contigs:
            return Assembly.GRCH38
        elif self == self.DECOY:
            return Assembly.HS38D1
        elif self in {self.EBV, self.UNCATEGORIZABLE}:
            return Assembly.OTHER
        else:
            raise ValueError(f"Unrecognized contig type: {self}")


class ContigTypeDesirabilities(NamedTuple):
    desired_contig_types: FrozenSet[ContigType]
    undesired_contig_types: FrozenSet[ContigType]
    unexpected_contig_types: FrozenSet[ContigType]

    @classmethod
    def create(cls) -> "ContigTypeDesirabilities":
        desired_contig_types = frozenset({
            ContigType.AUTOSOME,
            ContigType.X,
            ContigType.Y,
            ContigType.MITOCHONDRIAL,
            ContigType.EBV,
            ContigType.DECOY,
            ContigType.UNLOCALIZED,
            ContigType.UNPLACED,
        })
        undesired_contig_types = frozenset({
            ContigType.ALT,
            ContigType.FIX_PATCH,
            ContigType.NOVEL_PATCH,
        })
        unexpected_contig_types = frozenset({
            ContigType.UNCATEGORIZABLE,
        })
        result = ContigTypeDesirabilities(desired_contig_types, undesired_contig_types, unexpected_contig_types)
        result.validate()
        return result

    def validate(self) -> None:
        handleable_contig_types = set.union(
            set(self.desired_contig_types),
            set(self.undesired_contig_types),
            set(self.unexpected_contig_types),
        )
        if set(contig_type for contig_type in ContigType) != handleable_contig_types:
            raise ValueError(f"Not all contig types can be handled: handleable_contig_types={handleable_contig_types}")

        contig_types_categories_overlap = (
                self.desired_contig_types.intersection(self.undesired_contig_types)
                or self.desired_contig_types.intersection(self.unexpected_contig_types)
                or self.undesired_contig_types.intersection(self.unexpected_contig_types)
        )
        if contig_types_categories_overlap:
            raise ValueError(f"Contig type categories overlap (so desired, undesired and unexpected are not disjoint)")

    def get_expected_contig_types(self) -> Set[ContigType]:
        return set(self.desired_contig_types.union(self.undesired_contig_types))
