from enum import unique, Enum, auto


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
    DECOY = auto()
    UNLOCALIZED = auto()
    UNPLACED = auto()
    ALT = auto()
    FIX_PATCH = auto()
    NOVEL_PATCH = auto()
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
