"""Python-converted replacement for ``version3.2plus.ipynb``.

That notebook is a path-specific MS2 workflow variant. Its core matching logic is
covered by ``ms2_legacy_core.py``; the hard-coded batch cell is replaced by the
configurable pipeline imported below.
"""

from sphingolipid_toolkit.ms2_legacy_core import *  # noqa: F401,F403
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig, parse_file_indices, run_batch, run_one_file
