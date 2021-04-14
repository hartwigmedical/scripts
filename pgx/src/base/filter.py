from enum import Enum, auto


class Filter(Enum):
    PASS = auto()
    NO_CALL = auto()
    INFERRED_V37_REF_CALL = auto()
    PASS_BUT_REF_V38 = auto()
    V37_PASS_BUT_UNKNOWN = auto()
