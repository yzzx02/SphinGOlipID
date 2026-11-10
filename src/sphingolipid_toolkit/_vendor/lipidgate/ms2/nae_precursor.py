"""NAE evidence policy: one retained precursor peak, with binary scoring."""
from dataclasses import replace

from .models import FragmentRecord, LibraryRecord


def nae_precursor_fragment(record: LibraryRecord) -> FragmentRecord:
    return FragmentRecord(record.precursor_mz, record.adduct, "Precursor Ion")


def normalize_nae_record(record: LibraryRecord) -> LibraryRecord:
    if record.compound_class.strip().upper() != "NAE":
        return record
    fragment = nae_precursor_fragment(record)
    if record.fragments == [fragment]:
        return record
    return replace(record, fragments=[fragment])
