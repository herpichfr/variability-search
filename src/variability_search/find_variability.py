#!/bin/python3

import os
from astropy.io import fits, ascii
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
        self.raunit = u.deg
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
            # Check if FITS or ASCII file
            file_has_lodaded = False
            try:
                with fits.open(path) as hdul:
                    data = hdul[1].data
                    catalogue = pd.DataFrame(data.byteswap().newbyteorder())
                    catalogues.append(catalogue)
                file_has_lodaded = True
            except Exception as e:
                self.logger.warning(
                    f"Failed to load {file} as FITS. Trying ASCII. Error: {e}")
                try:
                    catalogue = ascii.read(path).to_pandas()
                    catalogues.append(catalogue)
                    file_has_lodaded = True
                except Exception as e:
                    self.logger.error(
                        f"Failed to load {file} as ASCII. Skipping file. Error: {e}")
            if not file_has_lodaded:
                self.logger.error(f"Could not load file: {file}")
                return
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
        If a star is not found in a catalogue, its corresponding columns
        shall be filled with NaN values.
        """
        self.logger.info("Organizing data from catalogs.")
        column_names_equiv = {
            'RA': ['RA', 'RightAscension', 'RAJ2000', 'ALPHA_J2000'],
            'DEC': ['DEC', 'Declination', 'DEJ2000', 'DELTA_J2000'],
            'MAG': ['MAG', 'Magnitude', 'MAG_AUTO'],
            'MAG_ERR': ['MAG_ERR', 'Magnitude_Error', 'MAGERR_AUTO'],
            'FLUX': ['FLUX', 'Flux', 'FLUX_AUTO'],
            'FLUX_ERR': ['FLUX_ERR', 'Flux_Error', 'FLUXERR_AUTO'],
            'FWHM': ['FWHM', 'FWHM_IMAGE'],
            'FLAGS': ['FLAGS', 'FLAGS_IMAGE'],
            'CLASS_STAR': ['CLASS_STAR', 'STAR_CLASS']
        }
        # Standardize column names across catalogues
        for i, catalogue in enumerate(catalogues):
            for standard_name, variants in column_names_equiv.items():
                for variant in variants:
                    if variant in catalogue.columns:
                        catalogue.rename(
                            columns={variant: standard_name}, inplace=True)
            required_columns = ['RA', 'DEC', 'MAG', 'MAG_ERR',
                                'FLUX', 'FLUX_ERR']
            missing_columns = [
                col for col in required_columns if col not in catalogue.columns]
            if missing_columns:
                self.logger.warning(
                    f"Catalog {i} is missing columns: {missing_columns}. Filling with NaN.")
                for col in missing_columns:
                    catalogue[col] = np.nan

        ra_dec = catalogues[0][['RA', 'DEC']].values
        for catalogue in catalogues[1:]:
            # gather all stars based on RA and DEC
            ra_dec = np.vstack(ra_dec, catalogue[['RA', 'DEC']].values)

        coords = SkyCoord(
            ra=ra_dec[:, 0], dec=ra_dec[:, 1], unit=(self.raunit, u.deg))
        # Internal match with a tolerance of 1 arcsec
        idx, d2d, _ = coords.match_to_catalog_sky(coords)
        matched = d2d.arcsec < 1.0
        unified_coordinates = coords[matched]
        for i, catalogue in enumerate(catalogues):
            # Match each catalogue to the unified coordinates
            cat_coords = SkyCoord(
                ra=catalogue['RA'].values * self.raunit,
                dec=catalogue['DEC'].values * u.deg)
            idx, d2d, _ = cat_coords.match_to_catalog_sky(unified_coordinates)
            matched = d2d.arcsec < 1.0
            # Add MAG, MAG_ERR, FLUX, FLUX_ERR columns to unified dataset
            self.unified_data[f'MAG_{i}'] = catalogue['MAG'].values[idx]
            self.unified_data[f'MAG_ERR_{
                i}'] = catalogue['MAG_ERR'].values[idx]
            self.unified_data[f'FLUX_{i}'] = catalogue['FLUX'].values[idx]
            self.unified_data[f'FLUX_ERR_{
                i}'] = catalogue['FLUX_ERR'].values[idx]

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
