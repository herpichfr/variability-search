#!/bin/python3

import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from multiprocessing import Pool


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


class VariabilityFinder:
    def __init__(self, input_dir, output_dir, work_dir=None):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.work_dir = work_dir if work_dir else input_dir
        self.logger = logger(os.path.join(
            self.work_dir, 'variability_finder.log'))
        self.unified_data = None
        self.reference_star = None

    def load_catalogs(self):
        self.logger.info("Loading catalogs from input directory.")
        catalog_files = [f for f in os.listdir(
            self.input_dir) if f.endswith('.fits')]
        catalogues = []
        for file in catalog_files:
            path = os.path.join(self.input_dir, file)
            self.logger.info(f"Loading catalog: {file}")
            with fits.open(path) as hdul:
                data = hdul[1].data
                catalogues.append(pd.DataFrame(data))
        return catalogues

    def organize_data(self, catalogues):
        """
        Match stars across catalogues and produce a unified dataset.
        The generated catalogue shall contain the following columns:
        - Star RA
        - Star DEC
        - MAG_i for i in N catalogues (if magnitudes are present)
        - MAG_ERR_i for i in N catalogues (if magnitude errors are present)
        - FLUX_i for i in N catalogues (if fluxes are present)
        - FLUX_ERR_i for i in N catalogues (if flux errors are present)
        If a star is not found in a catalogue, it is marked as NaN.
        """
        self.logger.info("Organizing data from catalogs.")
        # TODO: Implement the logic to match stars across catalogues
        # and create a unified dataset to be returned.

    def define_reference_star(self):
        """
        Define a reference star from the catalogues.
        The method uses the saturation level as a criterion for the brightest
        stars and cap at the 50th percentile to avoid low signal-to-noise
        stars.
        """
        # TODO: Implement the logic to define the reference star.
        # NOTE: Could use the average flux of the more stable and well constrained stars.
        # TODO: Return all catalogues with all stars divided by the reference star flux
        # or subtracted in case of magnitudes.
        # NOTE: The following definition is a placeholder.

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
