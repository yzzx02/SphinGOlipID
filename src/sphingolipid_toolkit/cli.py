"""Command-line interface for SphinGOlipID project modules."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import LEGACY_FRAGMENT_PPM, SphinGOlipIDConfig
from .logging_utils import configure_logging
from .ms2_pipeline import parse_file_indices, run_batch


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sphingolipid-ms2", description="Run SphinGOlipID MS2 matching workflow.")
    parser.add_argument(
        "--input-dir",
        "--raw-ms2-dir",
        dest="input_dir",
        required=True,
        type=Path,
        help="Folder containing raw MS/MS txt files. Historically this also held precursor Excel files.",
    )
    parser.add_argument(
        "--precursor-dir",
        type=Path,
        default=None,
        help="Folder containing precursor target Excel files. Defaults to --input-dir.",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Folder where result_msms-{i}.xlsx files will be written.")
    parser.add_argument("--files", default="1-6", help="File indices to process, e.g. 1-6 or 1,3,5.")
    parser.add_argument("--encoding", default="GBK", help="Encoding of exported MS/MS txt files. Default: GBK.")
    parser.add_argument("--fragment-ppm", type=float, default=LEGACY_FRAGMENT_PPM, help=f"Fragment m/z tolerance in ppm. Default: {LEGACY_FRAGMENT_PPM:g}.")
    parser.add_argument("--top-n", type=int, default=3, help="Keep top N candidates per precursor. Default: 3.")
    parser.add_argument("--ms2-library", type=Path, help="Optional validated sphingolipid SQLite index.")
    parser.add_argument("--precursor-ppm", type=float, default=10.)
    parser.add_argument("--fragment-da", type=float, help="Use an absolute MS2 tolerance instead of ppm.")
    parser.add_argument("--min-total-score", type=float, default=50., help="Reference engine total score threshold (0-100).")
    parser.add_argument(
        "--keep-intermediate",
        "--save-intermediate",
        dest="save_intermediate",
        action="store_true",
        help="Keep intermediate files under output/intermediate.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level. Default: INFO.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logger = configure_logging(args.log_level.upper(), logger_name="sphingolipid_toolkit.ms2")
    config = SphinGOlipIDConfig(
        raw_ms2_dir=args.input_dir,
        precursor_dir=args.precursor_dir or args.input_dir,
        output_dir=args.output_dir,
        ms1_library_path=None,
        file_indices=parse_file_indices(args.files),
        encoding=args.encoding,
        fragment_ppm=args.fragment_ppm,
        top_n=args.top_n,
        save_intermediate=args.save_intermediate,
        ms2_library_path=args.ms2_library,
        precursor_ppm=args.precursor_ppm,
        fragment_da=args.fragment_da,
        min_total_score=args.min_total_score,
    )
    results = run_batch(config, logger=logger)
    print("Generated result files:")
    for path in results:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
