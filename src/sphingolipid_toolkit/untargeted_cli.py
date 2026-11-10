"""Non-targeted mzML -> XCMS local features -> current MS2 engine -> XLSX."""
import argparse
import logging
from pathlib import Path

from .untargeted import UntargetedConfig, run_untargeted


def main(argv=None):
    parser=argparse.ArgumentParser(description='SphinGOlipID non-targeted XCMS mzML workflow')
    parser.add_argument('mzml',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--polarity',choices=('+','-'),default='+')
    parser.add_argument('--precursor-ppm',type=float,default=10.)
    tolerance=parser.add_mutually_exclusive_group()
    tolerance.add_argument('--fragment-ppm',type=float)
    tolerance.add_argument('--fragment-da',type=float)
    parser.add_argument('--min-score',type=float,default=50.)
    parser.add_argument('--top-n',type=int,default=1)
    parser.add_argument('--peakwidth-min',type=float,default=5.,help='seconds')
    parser.add_argument('--peakwidth-max',type=float,default=60.,help='seconds')
    parser.add_argument('--noise',type=float,default=1000.)
    parser.add_argument('--snthresh',type=float,default=3.)
    parser.add_argument('--rscript',type=Path)
    parser.add_argument('--save-intermediate',action='store_true')
    args=parser.parse_args(argv)
    try:
        config=UntargetedConfig(tuple(args.mzml),args.output,polarity=args.polarity,
            precursor_ppm=args.precursor_ppm,
            fragment_ppm=None if args.fragment_da is not None else args.fragment_ppm if args.fragment_ppm is not None else 20.,
            fragment_da=args.fragment_da,min_total_score=args.min_score,top_n=args.top_n,
            peakwidth_min=args.peakwidth_min,peakwidth_max=args.peakwidth_max,
            noise=args.noise,snthresh=args.snthresh,rscript_path=args.rscript,save_intermediate=args.save_intermediate)
    except ValueError as exc:
        parser.error(str(exc))
    logger=logging.getLogger('sphingolipid_toolkit.untargeted')
    logger.setLevel(logging.INFO)
    handler=logging.StreamHandler()
    logger.addHandler(handler)
    try:
        run_untargeted(config,logger)
    finally:
        logger.removeHandler(handler)
        handler.close()
    return 0


if __name__=='__main__':raise SystemExit(main())
