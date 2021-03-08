from enum import Enum, auto


class Filter(Enum):
    PASS = auto()
    NO_CALL = auto()
    INFERRED_GRCH37_REF_CALL = auto()
    PASS_BUT_REF_GRCH38 = auto()
