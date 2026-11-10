"""Identify ordinary mzML files without a targetlist or the desktop interface."""
import argparse
from pathlib import Path

from .engine import IdentificationConfig,SphingolipidIdentifier
from .result_writer import write_identifications


def main(argv=None):
    parser=argparse.ArgumentParser(description='SphinGOlipID sphingolipid identification')
    parser.add_argument('mzml',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--precursor-ppm',type=float,default=10.)
    tolerance=parser.add_mutually_exclusive_group()
    tolerance.add_argument('--fragment-ppm',type=float)
    tolerance.add_argument('--fragment-da',type=float)
    parser.add_argument('--min-score',type=float,default=50.)
    parser.add_argument('--top-n',type=int,default=3)
    parser.add_argument('--library',type=Path)
    args=parser.parse_args(argv)
    if not args.mzml.is_file():parser.error('Input mzML file does not exist')
    if args.output.suffix.lower()!='.xlsx':parser.error('Output must be an .xlsx file')
    if args.output.resolve()==args.mzml.resolve():parser.error('Output would overwrite the input')
    config=IdentificationConfig(precursor_ppm=args.precursor_ppm,
        fragment_ppm=None if args.fragment_da is not None else args.fragment_ppm if args.fragment_ppm is not None else 20.,
        fragment_da=args.fragment_da,min_total_score=args.min_score,top_n=args.top_n,library_path=args.library)
    with SphingolipidIdentifier(config) as identifier:
        results=identifier.identify_mzml(args.mzml)
    write_identifications(results,args.output)
    print(f'{len(results)} candidates written to {args.output}')
    return 0


if __name__=='__main__':raise SystemExit(main())
