"""Command-line interface for variability-search."""

import os
import argparse
import logging

from src.variability_search import __version__
from src.variability_search.prepare_catalogues import PrepareCatalogues
from src.variability_search.variability_search import VariabilitySearch

from src.variability_search.variability_analyser import VariabilityAnalyser


def logger(logfile=None, loglevel=logging.INFO):
    logger = logging.getLogger(__name__)

    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] @%(module)s.%(funcName)s() %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.setLevel(loglevel)

    if logfile is not None:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def parse_args():
    """Main entry point for the variability-search CLI."""
    parser = argparse.ArgumentParser(
        description="Search for variability in astronomical observations"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "inputdir",
        nargs="?",
        help="Path to input list of catalogues file (CSV, ASCII, or FITS)",
    )
    parser.add_argument(
        "--workdir",
        default=None,
        help="Path to working directory. Default is inputdir if not specified.",
    )
    parser.add_argument(
        "--outputdir",
        default=None,
        help="Path to output directory. Default is 'output'.",
    )
    parser.add_argument(
        "--input_catalogue",
        default="variability_results.csv",
        help="Path to input variability results file (CSV, ASCII, or FITS)",
    )
    parser.add_argument(
        "--raunit",
        default="deg",
        help="Unit for Right Ascension (default: deg)",
    )
    parser.add_argument(
        "--np",
        type=int,
        default=1,
        help="Number of processes to use for parallel processing (default: 1)",
    )
    parser.add_argument(
        "--min_epocs",
        type=int,
        default=20,
        help="Minimum number of epochs required for variability analysis (default: 20)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output file",
    )
    parser.add_argument(
        "--save_output",
        action="store_true",
        help="Save output files",
    )
    parser.add_argument(
        "--save_plots",
        action="store_true",
        help="Save diagnostic plots",
    )
    parser.add_argument(
        "--compute_fourier",
        action="store_true",
        help="Compute Fourier analysis",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        help="Path to log file",
    )
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode",
    )

    args = parser.parse_args()

    return args


def main():
    """Main entry point for the variability-search CLI."""
    args = parse_args()
    log = logger(logfile=args.logfile,
                 loglevel=getattr(logging, args.loglevel))
    workdir = args.workdir if args.workdir else args.inputdir

    if args.compute_fourier:
        if not os.path.join(workdir, "variability_results.csv"):
            # Step 1: Run prepare_catalogues.py
            log.info("Step 1: Preparing catalogues...")
            preparer = PrepareCatalogues(args, logger=log)
            unified_catalogue = preparer.run()

            # Step 2: Run variability_search.py
            log.info("Step 2: Searching for variability...")
            variability_searcher = VariabilitySearch(
                args, unified_catalogue=unified_catalogue, logger=log)
            variability_searcher.run()

            # Step 3: Run variability_analyser.py
            log.info("Step 3: Analysing variability...")
            variability_analyser = VariabilityAnalyser(
                args, logger=log)
            variability_analyser.run()
        else:
            log.info("Starting with step 3: Analysing variability...")
            variability_analyser = VariabilityAnalyser(
                args, logger=log)
            variability_analyser.run()


if __name__ == "__main__":
    main()
