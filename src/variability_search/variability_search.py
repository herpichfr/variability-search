#!/bin/python3

import argparse
import logging
import os

import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Search for variability in astronomical observations"
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
    def __init__(self,
                 args,
                 unified_catalogue: pd.DataFrame | None = None,
                 logger=logging.getLogger(__name__)
                 ):
        self.logger = logger
        self.inputdir = args.inputdir
        self.workdir = args.workdir if args.workdir else args.inputdir
        self.outputdir = args.outputdir
        self.raunit = args.raunit
        self.np = args.np
        self.output = args.output if args.output else os.path.join(
            self.outputdir, "variability_results.csv")
        self.debug = args.debug
        self.unified_catalogue = unified_catalogue

    def define_reference_star(self):
        """
        Define a reference star from the catalogues.
        The method uses the saturation level as a criterion for the brightest
        stars and cap at the 50th percentile to avoid low signal-to-noise
        stars.
        """
        # load unified dataset
        self.logger.info("Defining reference star from unified catalogue.")
        if self.unified_catalogue is None:
            unified_catalogue = pd.read_csv(os.path.join(
                self.workdir, "unified_catalogue.csv"))
        # Identify candidate reference stars based on brightness criteria
        _ref_column = 'MAG'
        _min_mag = 4.0
        _mag_columns = [
            col for col in self.unified_catalogue.columns if col.startswith(f'{_ref_column}_') and 'ERR' not in col]
        mask = (self.unified_catalogue[f'{_ref_column}_0'] >= _min_mag) & (
            self.unified_catalogue[f'{_ref_column}_0'] < np.percentile(
                self.unified_catalogue[f'{_ref_column}_0'][~np.isnan(self.unified_catalogue[f'{_ref_column}_0'])], 95))
        print('min mag for column', _mag_columns[0], min(self.unified_catalogue[_mag_columns[0]]), 'and max mag',
              np.percentile(self.unified_catalogue[_mag_columns[0]][~np.isnan(self.unified_catalogue[_mag_columns[0]])], 95))
        for col in _mag_columns[1:]:
            mask &= (self.unified_catalogue[col] >= 8.0) & (
                self.unified_catalogue[col] < np.percentile(
                    self.unified_catalogue[col][~np.isnan(self.unified_catalogue[col])], 95))
        candidate_stars = self.unified_catalogue[mask]
        self.logger.info(f"Found {len(candidate_stars)
                                  } candidate reference stars.")

        # Iterate over the selected candidate stars to remove variable stars or
        # outliers based on their light curves

        # Create a matrix with the light curves of the candidate stars
        light_curves = candidate_stars[_mag_columns].to_numpy()
        # Compute the standard deviation of each light curve
        std_devs = np.nanstd(light_curves, axis=1)
        # Remove the 5% of stars with the highest standard deviation
        threshold = np.percentile(std_devs, 95)
        stable_stars = candidate_stars[std_devs <= threshold]
        self.logger.info(f"Selected {len(stable_stars)
                                     } stable reference stars after variability filtering.")

        # Calculate the median light curve of the stable stars
        median_light_curve = np.nanmedian(
            stable_stars[_mag_columns].to_numpy(), axis=0)

        if self.debug:
            self.logger.debug(
                f"Median light curve of reference stars: {median_light_curve}")
            import pdb
            pdb.set_trace()

        return stable_stars, median_light_curve

    def select_comparison_stars(self,
                                stable_stars: pd.DataFrame,
                                reference_light_curve: np.ndarray,
                                variability_threshold: float = 0.01
                                ):
        """
        Select comparison stars based on variability threshold.
        The method is iterative to remove higher-than-threshold variability
        stars. The comparison stars correspond to the median magnitude or
        flux within bins of similar brightness.
        """
        # Find the 5% more stable stars after subtracting the reference light curve
        self.logger.info("Selecting comparison stars based on variability.")
        _ref_column = 'MAG'
        _mag_columns = [
            col for col in self.unified_catalogue.columns if col.startswith(f'{_ref_column}_') and 'ERR' not in col]
        light_curves = stable_stars[_mag_columns].to_numpy()
        # Subtract the reference light curve
        adjusted_light_curves = light_curves - reference_light_curve
        # Compute the standard deviation of each adjusted light curve
        std_devs = np.nanstd(adjusted_light_curves, axis=1)
        # Select stars below the variability threshold
        comparison_stars = stable_stars[std_devs <= variability_threshold]
        self.logger.info(f"Selected {len(comparison_stars)
                                     } comparison stars after applying variability threshold of {variability_threshold}.")
        if self.debug:
            import pdb
            pdb.set_trace()

        return comparison_stars

    def search_variability(self,
                           reference_light_curve: np.ndarray,
                           comparison_stars: pd.DataFrame
                           ):
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
        stable_stars, median_light_curve = self.define_reference_star()
        self.select_comparison_stars(stable_stars, median_light_curve)


if __name__ == "__main__":
    args = parse_arguments()
    variability_search = VariabilitySearch(args)
    variability_search.run()
