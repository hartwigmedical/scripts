from typing import NamedTuple, Tuple, List

from contig_name_translation import ContigNameTranslator
from contig_types import ContigType


class CategorizedContigNames(NamedTuple):
    autosomes: Tuple[str, ...]
    x_contigs: Tuple[str, ...]
    y_contigs: Tuple[str, ...]
    mitochondrial_contigs: Tuple[str, ...]
    ebv_contigs: Tuple[str, ...]
    decoys: Tuple[str, ...]
    unlocalized_contigs: Tuple[str, ...]
    unplaced_contigs: Tuple[str, ...]
    alt_contigs: Tuple[str, ...]
    fix_patch_contigs: Tuple[str, ...]
    novel_patch_contigs: Tuple[str, ...]

    def get_contig_names(self) -> Tuple[str, ...]:
        contig_names: List[str] = []
        contig_names.extend(self.autosomes)
        contig_names.extend(self.x_contigs)
        contig_names.extend(self.y_contigs)
        contig_names.extend(self.mitochondrial_contigs)
        contig_names.extend(self.ebv_contigs)
        contig_names.extend(self.decoys)
        contig_names.extend(self.unlocalized_contigs)
        contig_names.extend(self.unplaced_contigs)
        contig_names.extend(self.alt_contigs)
        contig_names.extend(self.fix_patch_contigs)
        contig_names.extend(self.novel_patch_contigs)
        return tuple(contig_names)


class ContigCategorizer(object):
    AUTOSOME_CONTIG_NAMES = {f"chr{n}" for n in range(1, 23)}
    X_CHROMOSOME_CONTIG_NAME = "chrX"
    Y_CHROMOSOME_CONTIG_NAME = "chrY"
    MITOCHONDRIAL_CONTIG_NAME = "chrM"
    EBV_CONTIG_NAME = "chrEBV"
    UNPLACED_PREFIX = "chrUn_"
    DECOY_SUFFIX = "_decoy"
    UNLOCALIZED_SUFFIX = "_random"
    FIX_PATCH_SUFFIX = "_fix"
    NOVEL_PATCH_SUFFIX = "_novel"
    ALT_SUFFIX = "_alt"

    def __init__(self, contig_name_translator: ContigNameTranslator) -> None:
        self._contig_name_translator = contig_name_translator

    def get_categorized_contig_names(
            self,
            contig_names: List[str],
    ) -> CategorizedContigNames:
        autosome_contigs: List[str] = []
        x_contigs: List[str] = []
        y_contigs: List[str] = []
        mitochondrial_contigs: List[str] = []
        ebv_contigs: List[str] = []
        decoy_contigs: List[str] = []
        unlocalized_contigs: List[str] = []
        unplaced_contigs: List[str] = []
        alt_contigs: List[str] = []
        fix_patch_contigs: List[str] = []
        novel_patch_contigs: List[str] = []
        uncategorized_contigs: List[str] = []

        for contig_name in contig_names:
            contig_type = self.get_contig_type(contig_name)

            if contig_type == ContigType.AUTOSOME:
                autosome_contigs.append(contig_name)
            elif contig_type == ContigType.X:
                x_contigs.append(contig_name)
            elif contig_type == ContigType.Y:
                y_contigs.append(contig_name)
            elif contig_type == ContigType.MITOCHONDRIAL:
                mitochondrial_contigs.append(contig_name)
            elif contig_type == ContigType.EBV:
                ebv_contigs.append(contig_name)
            elif contig_type == ContigType.DECOY:
                decoy_contigs.append(contig_name)
            elif contig_type == ContigType.UNLOCALIZED:
                unlocalized_contigs.append(contig_name)
            elif contig_type == ContigType.UNPLACED:
                unplaced_contigs.append(contig_name)
            elif contig_type == ContigType.ALT:
                alt_contigs.append(contig_name)
            elif contig_type == ContigType.FIX_PATCH:
                fix_patch_contigs.append(contig_name)
            elif contig_type == ContigType.NOVEL_PATCH:
                novel_patch_contigs.append(contig_name)
            elif contig_type == ContigType.UNCATEGORIZABLE:
                uncategorized_contigs.append(contig_name)
            else:
                raise ValueError(f"Unrecognized contig type: {contig_type}")

        if uncategorized_contigs:
            raise ValueError(f"Uncategorized contigs: {uncategorized_contigs}")

        categorized_contigs = CategorizedContigNames(
            tuple(autosome_contigs),
            tuple(x_contigs),
            tuple(y_contigs),
            tuple(mitochondrial_contigs),
            tuple(ebv_contigs),
            tuple(decoy_contigs),
            tuple(unlocalized_contigs),
            tuple(unplaced_contigs),
            tuple(alt_contigs),
            tuple(fix_patch_contigs),
            tuple(novel_patch_contigs),
        )

        if set(categorized_contigs.get_contig_names()) != set(contig_names):
            error_msg = (
                f"Categorizing contigs failed. Mismatch between input and output: "
                f"input={contig_names}, output={categorized_contigs.get_contig_names()}"
            )
            raise ValueError(error_msg)

        return categorized_contigs

    def get_contig_type(self, contig_name: str) -> ContigType:
        matching_contig_types = []
        if self.is_autosome_contig_name(contig_name):
            matching_contig_types.append(ContigType.AUTOSOME)
        if self.is_x_contig_name(contig_name):
            matching_contig_types.append(ContigType.X)
        if self.is_y_contig_name(contig_name):
            matching_contig_types.append(ContigType.Y)
        if self.is_mitochondrial_contig_name(contig_name):
            matching_contig_types.append(ContigType.MITOCHONDRIAL)
        if self.is_ebv_contig_name(contig_name):
            matching_contig_types.append(ContigType.EBV)
        if self.is_decoy_contig_name(contig_name):
            matching_contig_types.append(ContigType.DECOY)
        if self.is_unlocalized_contig_name(contig_name):
            matching_contig_types.append(ContigType.UNLOCALIZED)
        if self.is_unplaced_contig_name(contig_name):
            matching_contig_types.append(ContigType.UNPLACED)
        if self.is_alt_contig_name(contig_name):
            matching_contig_types.append(ContigType.ALT)
        if self.is_fix_patch_contig_name(contig_name):
            matching_contig_types.append(ContigType.FIX_PATCH)
        if self.is_novel_patch_contig_name(contig_name):
            matching_contig_types.append(ContigType.NOVEL_PATCH)

        if len(matching_contig_types) == 1:
            return matching_contig_types[0]
        elif not matching_contig_types:
            return ContigType.UNCATEGORIZABLE
        else:
            raise ValueError(f"Contig matches multiple contig types: {matching_contig_types}")

    def is_autosome_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        return standardized_contig_name in self.AUTOSOME_CONTIG_NAMES

    def is_x_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        return standardized_contig_name == self.X_CHROMOSOME_CONTIG_NAME

    def is_y_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        return standardized_contig_name == self.Y_CHROMOSOME_CONTIG_NAME

    def is_mitochondrial_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        return standardized_contig_name == self.MITOCHONDRIAL_CONTIG_NAME

    def is_ebv_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        return standardized_contig_name == self.EBV_CONTIG_NAME

    def is_decoy_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and standardized_contig_name.endswith(self.DECOY_SUFFIX)
        )
        return result

    def is_unlocalized_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            not standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and standardized_contig_name.endswith(self.UNLOCALIZED_SUFFIX)
        )
        return result

    def is_unplaced_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and not standardized_contig_name.endswith(self.DECOY_SUFFIX)
        )
        return result

    def is_alt_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            not standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and standardized_contig_name.endswith(self.ALT_SUFFIX)
        )
        return result

    def is_fix_patch_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            not standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and standardized_contig_name.endswith(self.FIX_PATCH_SUFFIX)
        )
        return result

    def is_novel_patch_contig_name(self, contig_name: str) -> bool:
        standardized_contig_name = self._contig_name_translator.standardize(contig_name)
        result = (
            not standardized_contig_name.startswith(self.UNPLACED_PREFIX)
            and standardized_contig_name.endswith(self.NOVEL_PATCH_SUFFIX)
        )
        return result
