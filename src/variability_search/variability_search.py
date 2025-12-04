#!/bin/python3

import os
import sys
import argparse
import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Search for variability in astronomical observations"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0",
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
        default="output",
        help="Path to output directory. Default is 'output'.",
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
        "-o",
        "--output",
        help="Path to output file",
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


class VariabilitySearch:
    def __init__(self, args, logger=logging.getLogger(__name__)):
        self.logger = logger
        self.inputdir = args.inputdir
        self.workdir = args.workdir if args.workdir else args.inputdir
        self.outputdir = args.outputdir
        self.raunit = args.raunit
        self.np = args.np

    def define_reference_star(self):
        """
        Define a reference star from the catalogues.
        The method uses the saturation level as a criterion for the brightest
        stars and cap at the 50th percentile to avoid low signal-to-noise
        stars.
        """
        # TODO: Implement the logic to define the reference star.
        # NOTE: Could use the average flux of the more stable and well constrained stars.
        # load unified dataset
        self.logger.info("Defining reference star from unified catalogue.")
        unified_catalogue = pd.read_csv(os.path.join(
            self.workdir, "unified_catalogue.csv"))
        # Identify candidate reference stars based on brightness criteria
        mask = (unified_catalogue['FLUX'] < 50000) & (
            unified_catalogue['FLUX'] > np.percentile(unified_catalogue['FLUX'], 50))
        import pdb
        pdb.set_trace()

    def select_comparison_stars(self, threshold=0.01):
        """
        Select comparison stars based on variability threshold.
        The method is iterative to remove higher-than-threshold variability
        stars. The comparison stars correspond to the median magnitude or
        flux within bins of similar brightness.
        """

    def search_variability(self):
        """
        Search for variability in the unified dataset.
        The method computes the variability metrics for each star and
        identifies variable stars based on predefined criteria.
        The stars meeting the variability criteria are saved to an output file
        per star both as a CSV with FLUX, FLUX_ERR and/or MAG and MAG_ERR
        columns and MAG - MAG_COMP and a plot showing the light curve.
        """

    def run(self):
        self.logger.info("Starting variability search process.")
        self.define_reference_star()


if __name__ == "__main__":
    args = parse_args()
    variability_search = VariabilitySearch(args)
    variability_search.run()
