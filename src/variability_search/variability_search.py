#!/bin/python3

import argparse
import logging
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
        self.outputdir = args.outputdir if args.outputdir else self.workdir
        self.raunit = args.raunit
        self.np = args.np
        self.output = args.output if args.output else os.path.join(
            self.outputdir, "variability_results.csv")
        self.save_output = True if args.output else False
        self.debug = args.debug
        self.unified_catalogue = unified_catalogue
        self.ref_column_name = 'FLUX'
        self.flux_column_names = list()

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
            self.unified_catalogue = pd.read_csv(os.path.join(
                self.workdir, "unified_catalogue.csv"))
        # Identify candidate reference stars based on brightness criteria
        # Use mags to determine the limits in fluxes
        _min_mag = 4.0
        _max_mag = np.percentile(
            self.unified_catalogue['MAG_0'][~np.isnan(self.unified_catalogue['MAG_0'])], 95)
        mask = (self.unified_catalogue['MAG_0'] >= _min_mag) & (
            self.unified_catalogue['MAG_0'] < _max_mag)
        _min_flux, _max_flux = self.unified_catalogue['FLUX_0'][mask].min(
        ), self.unified_catalogue['FLUX_0'][mask].max()

        self.flux_columnames = [
            col for col in self.unified_catalogue.columns if col.startswith(f'{self.ref_column_name}_') and 'ERR' not in col]

        mask = (self.unified_catalogue[f'{self.ref_column_name}_0'] > _min_flux) & (
            self.unified_catalogue[f'{self.ref_column_name}_0'] <= _max_flux)

        print('min flux for column',
              self.flux_columnames[0],
              min(self.unified_catalogue[self.flux_columnames[0]]),
              'and max flux',
              np.percentile(self.unified_catalogue[self.flux_columnames[0]][~np.isnan(self.unified_catalogue[self.flux_columnames[0]])], 95))

        for col in self.flux_columnames[1:]:
            mask &= (self.unified_catalogue[col] > _min_flux) & (
                self.unified_catalogue[col] <= _max_flux)
        candidate_stars = self.unified_catalogue[mask]
        self.logger.info(f"Found {len(candidate_stars)
                                  } candidate reference stars.")

        # Iterate over the selected candidate stars to remove variable stars or
        # outliers based on their light curves

        # Create a matrix with the light curves of the candidate stars
        light_curves = candidate_stars[self.flux_columnames].to_numpy()
        # Compute the standard deviation of each light curve
        std_devs = np.nanstd(light_curves, axis=1)
        # Remove the 5% of stars with the highest standard deviation
        threshold = np.percentile(std_devs, 95)
        stable_stars = candidate_stars[std_devs <= threshold]
        self.logger.info(f"Selected {len(stable_stars)
                                     } stable reference stars after variability filtering.")

        # Calculate the median light curve of the stable stars
        median_light_curve = np.nanmedian(
            stable_stars[self.flux_columnames].to_numpy(), axis=0)

        if self.debug:
            self.logger.debug(
                f"Median light curve of reference stars: {median_light_curve}")
            import pdb
            pdb.set_trace()

        return stable_stars, median_light_curve

    def select_comparison_stars(self,
                                stable_stars: pd.DataFrame,
                                reference_light_curve: np.ndarray,
                                variability_threshold: float = 0.02  # TODO: Verify appropriate threshold
                                ):
        """
        Select comparison stars based on variability threshold.
        The method is iterative to remove higher-than-threshold variability
        stars. The comparison stars correspond to the median magnitude or
        flux within bins of similar brightness.
        """
        # Find the 5% more stable stars after subtracting the reference light curve
        self.logger.info("Selecting comparison stars based on variability.")
        light_curves = stable_stars[self.flux_columnames].to_numpy()
        # Subtract the reference light curve
        adjusted_light_curves = light_curves / reference_light_curve
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
        # Clean unified_catalogue to exclude stars that have less than two observations that are not NaN
        self.logger.info("Searching for variability in the unified catalogue.")
        unified_catalogue = self.unified_catalogue.copy()
        valid_obs_mask = unified_catalogue[self.flux_columnames].notna().sum(
            axis=1) >= 2
        unified_catalogue = unified_catalogue[valid_obs_mask]
        # Compute the average light curve of the comparison stars
        comparison_light_curve = np.nanmedian(
            comparison_stars[self.flux_columnames].to_numpy(), axis=0)
        # Iterate over all stars in the unified catalogue to compute variability
        variability_results = pd.DataFrame(
            columns=[f'REL_FLUX_{i}' for i in range(len(self.flux_columnames))] + ['STD_DEV'])
        for index, star in unified_catalogue[:2].iterrows():
            star_light_curve = star[self.flux_columnames].to_numpy()
            # Subtract the reference light curve from the star's light curve
            star_light_curve = star_light_curve / reference_light_curve
            # Compute the differential light curve
            differential_light_curve = star_light_curve - comparison_light_curve
            # Compute variability metrics (e.g., standard deviation)
            std_dev = np.nanstd(differential_light_curve)
            variability_results.loc[index] = list(
                differential_light_curve) + [std_dev]
            fig = plt.figure(figsize=(10, 6))
            plt.errorbar(range(len(star_light_curve)), star_light_curve,
                         yerr=[star[f'FLUX_ERR_{i}'] if f'FLUX_ERR_{
                             i}' in star else 0 for i in range(len(self.flux_columnames))],
                         fmt='o', label='Star Light Curve')
            plt.errorbar(range(len(reference_light_curve)), reference_light_curve,
                         fmt='o', label='Reference Light Curve')
            plt.errorbar(range(len(comparison_light_curve)), comparison_light_curve,
                         fmt='o', label='Comparison Light Curve')
            plt.xlabel('Observation Index')
            plt.ylabel('Relative Flux')
            plt.title(f'Variability Analysis for Star index {index}')
            plt.legend()
            plt.grid()
            plt.show()
            if self.debug:
                ask = input("Continue showing plots? (y/n): ")
                if ask.lower() != 'y':
                    plt.close('all')
                    break
            plt.close(fig)

        # Merge the two dataframes
        final_results = pd.concat(
            [unified_catalogue.reset_index(drop=True), variability_results.reset_index(drop=True)], axis=1)
        import pdb
        pdb.set_trace()
        if self.save_output:
            # Save the results to the output file
            final_results.to_csv(self.output, index=False)
            self.logger.info(
                f"Variability results saved to {self.output}.")

    def run(self):
        self.logger.info("Starting variability search process.")
        stable_stars, median_light_curve = self.define_reference_star()
        comparison_stars = self.select_comparison_stars(
            stable_stars, median_light_curve)
        self.search_variability(median_light_curve, comparison_stars)
        self.logger.info("Variability search process completed.")


if __name__ == "__main__":
    args = parse_arguments()
    variability_search = VariabilitySearch(args)
    variability_search.run()
