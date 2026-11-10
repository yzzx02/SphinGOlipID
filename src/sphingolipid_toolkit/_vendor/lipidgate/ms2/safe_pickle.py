"""Read existing library encodings without accepting executable pickle globals.

Indexes contain primitive values, array columns and occasionally NumPy scalars.
Legacy caches additionally contain the two library record classes. Keep their
numeric types and occurrence order intact while rejecting other constructors.
"""
from __future__ import annotations

import array
import io
import pickle

from .models import FragmentRecord, LibraryRecord


def _numpy_dtype(*args):
    import numpy as np

    dtype = np.dtype(*args)
    if dtype.hasobject:
        raise pickle.UnpicklingError("Object dtypes are not supported in library caches")
    return dtype


def _numpy_scalar(dtype, value):
    import numpy as np

    if not isinstance(dtype, np.dtype) or dtype.hasobject:
        raise pickle.UnpicklingError("Object scalars are not supported in library caches")
    try:
        from numpy._core.multiarray import scalar
    except ImportError:  # NumPy 1.x caches remain readable.
        from numpy.core.multiarray import scalar
    return scalar(dtype, value)


def _legacy_bytes(value, encoding, errors="strict"):
    if not isinstance(value, str) or encoding not in {"latin1", "latin-1"} or errors != "strict":
        raise pickle.UnpicklingError("Unsupported legacy byte encoding in library cache")
    return value.encode("latin1")


_GLOBALS = {
    ("sphingolipid_toolkit._vendor.lipidgate.ms2.models", "FragmentRecord"): FragmentRecord,
    ("sphingolipid_toolkit._vendor.lipidgate.ms2.models", "LibraryRecord"): LibraryRecord,
    ("array", "array"): array.array,
    ("array", "_array_reconstructor"): array._array_reconstructor,
    ("builtins", "set"): set,
    ("builtins", "frozenset"): frozenset,
    ("builtins", "complex"): complex,
    ("__builtin__", "set"): set,
    ("__builtin__", "frozenset"): frozenset,
    ("__builtin__", "complex"): complex,
    ("_codecs", "encode"): _legacy_bytes,
    ("numpy", "dtype"): _numpy_dtype,
    ("numpy.core.multiarray", "scalar"): _numpy_scalar,
    ("numpy._core.multiarray", "scalar"): _numpy_scalar,
}


class _LibraryUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            return _GLOBALS[module, name]
        except KeyError:
            raise pickle.UnpicklingError(
                f"Unsupported library pickle global: {module}.{name}"
            ) from None


def load(handle):
    return _LibraryUnpickler(handle).load()


def loads(payload):
    return load(io.BytesIO(payload))
