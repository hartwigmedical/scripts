from enum import Enum, auto


class Filter(Enum):
    PASS = auto()
    FILTERED = auto()
    NO_CALL = auto()
    INFERRED_REF_CALL = auto()
