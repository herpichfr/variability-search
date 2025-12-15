#!/bin/python3

import os
import logging
import argparse

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle
from astropy import units as u


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Analyse variability in astronomical observations"
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

    args = parser.parse_args()

    return args


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


class VariabilityAnalyser:
    def __init__(self, args, logger=logging.getLogger(__name__)):
        self.logger = logger
        self.inputdir = args.inputdir
        self.workdir = args.workdir if args.workdir else args.inputdir
        self.outputdir = args.outputdir if args.outputdir else self.workdir
        self.input_catalogue = os.path.join(self.workdir, args.input_catalogue)
        self.raunit = args.raunit
        self.np = args.np
        self.save_output = args.save_output
        self.save_plots = args.save_plots
        self.debug = args.debug

    def run_fourier_analysis(self, variability_results):
        """
            Run Fourier analysis on variability results

            Parameters
            ----------
            variability_results : pd.DataFrame
                DataFrame containing variability results

            Returns
            -------
            pd.DataFrame
                DataFrame containing Fourier analysis results
        """
        _flux_columnames = [col for col in variability_results.columns if col.startswith(
            'FLUX') and 'ERR' not in col]
        # Load observations dates
        obs_dates = pd.read_csv(os.path.join(
            self.workdir, "observations_dates.csv"))
        if len(obs_dates) == len(_flux_columnames):
            self.logger.info(
                "Observation dates loaded successfully for variability analysis.")
            x_labels = np.float64(obs_dates['OBS_DATE'])
        else:
            x_labels = range(len(_flux_columnames))
            self.logger.warning(
                "Observation dates not found or do not match number of observations. Using index as x-axis labels.")
        stars_light_curve = variability_results[_flux_columnames].to_numpy()
        frequency = np.linspace(3.47e-5, 0.006, len(_flux_columnames)) * u.Hz
        times = (x_labels - np.min(x_labels))
        # Compute the power spectrum of the light curves considering that many values may be NaN
        nan_mask = np.isnan(stars_light_curve)
        stars_light_curve = np.ma.masked_array(
            stars_light_curve, mask=nan_mask)
        stars_std = np.nanstd(stars_light_curve, axis=1)
        for i in range(len(stars_light_curve)):
            ls_model = LombScargle(
                times[~nan_mask[i]] * u.day,
                stars_light_curve[i].data[~nan_mask[i]],
                stars_light_curve[i].data[~nan_mask[i]] * [0.07],
            )
            power_spectrum = ls_model.power(frequency)
            peaks = find_peaks(power_spectrum, height=0.4)
            if peaks[0].size == 0:
                self.logger.warning(
                    f"No peaks found for star {i}. Skipping plot.")
            else:
                self.logger.info(
                    f"Found {peaks[0].size} peaks for star {i}.")
                best_peak = peaks[0][np.argmax(peaks[1]['peak_heights'])]
                best_frequency = frequency[best_peak]
                best_period = 1.0 / best_frequency.value * u.day
                t_fit = np.linspace(
                    0, np.max(times) - np.min(times), 1000) * u.day
                # phase_fit = t_fit / best_period
                y_fit = ls_model.model(t_fit, frequency[best_peak])
                plt.plot(t_fit, y_fit, lw=2, zorder=2,
                         label=f'Period: {best_period / u.day / 3600:.2f} h')
            plt.errorbar(
                times[~nan_mask[i]],
                stars_light_curve[i].data[~nan_mask[i]],
                yerr=stars_std[i],
                fmt='o', c='k', zorder=1)
            plt.legend()
            plt.show()
            import pdb
            pdb.set_trace()

        return

    def run(self):
        self.logger.info("Starting variability analyser process.")
        variability_results = pd.read_csv(self.input_catalogue)
        self.run_fourier_analysis(variability_results)
        self.logger.info("Variability analyser process completed.")


if __name__ == "__main__":
    args = parse_arguments()
    variability_analyser = VariabilityAnalyser(args)
