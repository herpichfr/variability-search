#!/bin/python3

import argparse
import logging
import os

from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle
from scipy.signal import find_peaks


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
        self.save_output = args.save_output
        self.save_plots = args.save_plots
        self.debug = args.debug
        self.unified_catalogue = unified_catalogue
        self.ref_column_name = 'FLUX'
        self.flux_columnames = []
        self.istest = args.istest

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
        _ref_flux = self.unified_catalogue[f'{self.ref_column_name}_0']
        _min_flux = np.percentile(_ref_flux[~np.isnan(_ref_flux)], 5)
        _max_flux = np.percentile(_ref_flux[~np.isnan(_ref_flux)], 99)

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
                                variability_threshold: float = 0.05  # TODO: Verify appropriate threshold
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

    def calculate_lombscargle(self, star_light_curve, star_ra, star_dec, obs_dates):
        frequency = np.linspace(
            1.16e-6, 0.0056, len(self.flux_columnames)) * u.Hz
        times = (np.float64(obs_dates['OBS_DATE']) -
                 np.min(np.float64(obs_dates['OBS_DATE']))) * u.day
        # star_light_curve = star[self.flux_columnames].to_numpy()
        nan_mask = np.isnan(star_light_curve)
        star_light_curve = np.ma.masked_array(
            star_light_curve, mask=nan_mask)
        # star_std = np.nanstd(star_light_curve, axis=1)
        ls_model = LombScargle(
            times[~nan_mask],
            star_light_curve.data[~nan_mask],
            star_light_curve.data[~nan_mask].size * [0.07],
        )
        power_spectrum = ls_model.power(frequency)
        peaks = find_peaks(power_spectrum, height=0.4)
        if peaks[0].size == 0:
            self.logger.warning(
                f"No peaks found for star RA: {star_ra:.4f}, DEC: {star_dec:.4f}.")
            t_fit, y_fit, best_period, best_peak_height = None, None, None, None
        else:
            self.logger.info(
                f"Found {peaks[0].size} peaks for star RA: {star_ra:.4f}, DEC: {star_dec:.4f}.")
            best_peak_height = np.max(peaks[1]['peak_heights'])
            best_peak = peaks[0][np.argmax(peaks[1]['peak_heights'])]
            best_frequency = frequency[best_peak]
            best_period = 1.0 / best_frequency.value * u.day
            t_fit = np.linspace(
                np.min(times), np.max(times), 1000)
            # phase_fit = t_fit / best_period
            try:
                y_fit = ls_model.model(t_fit, frequency[best_peak])
            except ValueError:
                self.logger.warning(
                    f"Failed to fit LombScargle model for star RA: {star_ra:.4f}, DEC: {star_dec:.4f}.")
                t_fit, y_fit, best_period, best_peak_height = None, None, None, None
                import pdb
                pdb.set_trace()

        return t_fit, y_fit, best_period, best_peak_height

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

        # Load obs_info catalogue to get observation dates
        obs_dates = pd.read_csv(os.path.join(
            self.workdir, "observations_dates.csv"))
        if len(obs_dates) == len(self.flux_columnames):
            self.logger.info(
                "Observation dates loaded successfully for variability analysis.")
            x_labels = np.float64(obs_dates['OBS_DATE'])
        else:
            x_labels = range(len(self.flux_columnames))
            self.logger.warning(
                "Observation dates not found or do not match number of observations. Using index as x-axis labels.")

        # If stars have less than one third of valid observations, remove them
        min_valid_obs = len(self.flux_columnames) // 3
        valid_obs_mask = unified_catalogue[self.flux_columnames].notna().sum(
            axis=1) >= min_valid_obs
        unified_catalogue = unified_catalogue[valid_obs_mask]
        self.logger.info(f"Filtered unified catalogue to {len(unified_catalogue)
                                                          } stars with at least {min_valid_obs} valid observations.")
        # Dive comparison stars light curve by reference light curve
        comparison_stars[self.flux_columnames] = comparison_stars[self.flux_columnames].div(
            reference_light_curve, axis=1)

        # Divide all stars by the reference light curve
        stars = unified_catalogue[self.flux_columnames].div(
            reference_light_curve, axis=1)

        stars_dfs = []
        for comp_star in comparison_stars[self.flux_columnames].itertuples(index=False):
            comp_star_series = pd.Series(comp_star, index=self.flux_columnames)
            stars_div_comp = stars.div(comp_star_series, axis=1)
            stars_dfs.append(stars_div_comp)

        median_comparison_light_curve = np.nanmedian(
            comparison_stars[self.flux_columnames].to_numpy(), axis=0)

        # Compute the median of each star in stars_dfs
        median_stars_lcs = np.nanmedian(
            np.array([df[self.flux_columnames].to_numpy() for df in stars_dfs]), axis=0)
        median_stars_df = pd.DataFrame(
            unified_catalogue.copy())
        for i, col in enumerate(self.flux_columnames):
            median_stars_df[col] = median_stars_lcs[:, i]

        best_periods = np.zeros(len(median_stars_df))
        amplitudes = np.zeros(len(median_stars_df))
        best_peaks = np.zeros(len(median_stars_df))
        # Iterate over all stars in the unified catalogue to compute variability
        j = 0
        for index, star in median_stars_df.iterrows():
            star_ra = star['RA']
            star_dec = star['DEC']
            output_plotsdir = os.path.join(self.outputdir, 'plots')
            os.makedirs(output_plotsdir, exist_ok=True)
            plot_path = os.path.join(
                output_plotsdir, f'star_{star_ra:.5f}_{star_dec:.5f}_{index:05d}.png')
            if os.path.exists(plot_path) and self.save_plots:
                self.logger.info(
                    f"Plot for star index {index} already exists. Skipping."
                )
                continue
            # NOTE: delete the following lines to consider all stars
            # target_var_star_coords = SkyCoord(
            #     ra='23 34 15.0857248317', dec='-42 03 41.047972591', unit=('hour', 'deg')) # WASP 4b
            if self.istest:
                target_var_star_coords = SkyCoord(
                    # WASP 145A
                    ra='21 29 00.6895374360', dec='-58 50 10.299523848', unit=('hour', 'deg'))
                star_coords = SkyCoord(
                    ra=star_ra, dec=star_dec, unit=(self.raunit, 'deg'))
                separation = star_coords.separation(
                    target_var_star_coords).arcsec
                if separation > 5.0:
                    continue

            star_median_mag = np.nanmedian(
                star[[f'MAG_{i}' for i in range(len(self.flux_columnames))]].to_numpy())
            median_star_light_curve = star[self.flux_columnames].to_numpy()
            scaled_median_star_light_curve = median_star_light_curve / \
                np.nanmedian(median_star_light_curve)
            # Compute a floating median of the median star light curve with a window of 5
            window_size = 50 if len(self.flux_columnames) >= 50 else max(
                5, len(self.flux_columnames) // 10)
            floating_median = pd.Series(median_star_light_curve).rolling(
                window=window_size, center=True, min_periods=1).median().to_numpy()

            # Create subplot of a thirds of the figure height for the light curve
            _, (a0, a1) = plt.subplots(2, 1, sharex=True, gridspec_kw={
                'height_ratios': [3, 1]}, figsize=(10, 8))
            a0.errorbar(x_labels,
                        reference_light_curve /
                        np.nanmedian(reference_light_curve) + 1,
                        c='r', fmt='s', label='Reference Light Curve')
            a0.errorbar(x_labels,
                        median_comparison_light_curve /
                        np.nanmedian(median_comparison_light_curve) + 0.5,
                        c='g', fmt='*', label=f'N comparisons: {len(comparison_stars)}')
            _star_std = np.std(scaled_median_star_light_curve[~np.isnan(
                scaled_median_star_light_curve)])
            a0.errorbar(x_labels,
                        scaled_median_star_light_curve,
                        c='k', fmt='o', label=f'Star LC. STD: {_star_std:.3f}')
            a1.errorbar(x_labels,
                        scaled_median_star_light_curve,
                        c='k', fmt='o')
            for i, star_df in enumerate(stars_dfs):
                comp_light_curve = comparison_stars[self.flux_columnames].iloc[i].to_numpy(
                )
                a0.errorbar(x_labels,
                            comp_light_curve /
                            np.nanmedian(comp_light_curve) + 0.5,
                            c='g', fmt='*', alpha=0.3, zorder=-1)
                star_light_curve = star_df.loc[index].to_numpy()
                a0.errorbar(x_labels,
                            star_light_curve / np.nanmedian(star_light_curve),
                            c='gray', fmt='o', alpha=0.3, zorder=-2)
                a1.errorbar(x_labels,
                            star_light_curve / np.nanmedian(star_light_curve),
                            c='gray', fmt='o', alpha=0.3, zorder=-2)
            # Plot the floating median over the star light curve
            _floating_median_scaked = floating_median[~np.isnan(
                scaled_median_star_light_curve)] / np.nanmedian(floating_median)
            a1.plot(x_labels[~np.isnan(scaled_median_star_light_curve)],
                    _floating_median_scaked,
                    c='orange', ls='-', lw=2,
                    label='Floating Median', zorder=3)
            t_fit, y_fit, best_period, best_peak = self.calculate_lombscargle(
                _floating_median_scaked, star_ra, star_dec, obs_dates[~np.isnan(scaled_median_star_light_curve)])
            if t_fit is not None and y_fit is not None:
                t_fit = t_fit / u.day + \
                    np.min(obs_dates['OBS_DATE'][~np.isnan(
                        scaled_median_star_light_curve)])
                amplitude = (np.nanmax(y_fit) - np.nanmin(y_fit)) / 2
                a1.plot(t_fit, y_fit, lw=3, c='m', ls='--',
                        label=f'Period: {best_period / u.day /
                                         3600:.2f} h Amp: {amplitude:.3f} Peak: {best_peak:.2f}',
                        zorder=4)
                best_periods[j] = best_period / u.day
                amplitudes[j] = amplitude
                best_peaks[j] = best_peak
            j += 1

            a0.set_ylim(0.75, 2.25)
            a0.set_ylabel('Relative Flux')
            # Remove x labels from a0
            a0.set_xticklabels([])
            a0.set_title(f'Variability Analysis for Star index {index:05d} (RA: {
                star_ra:.5f}, DEC: {star_dec:.5f}) Mag: {star_median_mag:.1f}')
            a0.legend(loc='upper right')
            a0.grid()

            # Set y limits for a1 as the fifth and 95th percentiles of the star light curve
            _min_ylim = np.percentile(
                scaled_median_star_light_curve[~np.isnan(scaled_median_star_light_curve)], 3)
            _max_ylim = np.percentile(
                scaled_median_star_light_curve[~np.isnan(scaled_median_star_light_curve)], 97)
            a1.set_ylim(_min_ylim, _max_ylim)
            a1.set_ylabel('Relative Flux')
            a1.set_xlabel('Observation Date' if len(obs_dates) ==
                          len(self.flux_columnames) else 'Observation Index')
            a1.set_xticklabels(a1.get_xticks())
            a1.legend(loc='upper right')
            a1.grid()

            # Remove spaces between subplots
            plt.subplots_adjust(hspace=0.01)
            if self.save_plots:
                print(f'Saving plot {index} of {len(median_stars_df)}')
                plt.savefig(plot_path)
                self.logger.info(f"Saved variability plot for star index {
                                 index} to {plot_path}.")
                plt.close()
            else:
                plt.show()
                ask = input("Press Enter to continue, or 'q' to quit: ")
                if ask.lower() == 'q':
                    plt.close()
                    break

        if self.debug:
            import pdb
            pdb.set_trace()

        # Save median_stars_df to output file
        if self.save_output:
            os.makedirs(self.outputdir, exist_ok=True)
            median_stars_df['BEST_PERIOD_DAYS'] = best_periods
            median_stars_df['AMPLITUDE'] = amplitudes
            median_stars_df['BEST_PEAK'] = best_peaks
            median_stars_df.to_csv(self.output, index=False)
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
