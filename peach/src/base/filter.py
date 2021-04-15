from enum import Enum, auto


class Filter(Enum):
    PASS = auto()
    NO_CALL = auto()
    UNKNOWN = auto()
    INFERRED_PASS = auto()
