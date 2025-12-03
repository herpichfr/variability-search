#!/bin/python3

import os
import argparse
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from multiprocessing import Pool
from src.variability_search import __version__


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
    def __init__(self, args, logger=None):
        self.input_dir = args.inputdir
        self.output_dir = args.outputdir
        self.work_dir = args.workdir if args.workdir else self.input_dir
        self.raunit = u.Unit(args.raunit)
        self.numproc = args.np
        self.debug = args.debug
        self.logger = logger if logger else logging.getLogger(__name__)
        self.unified_catalogue = dict()
        self.reference_star = None
        self.coords = None

    def load_catalogs(self):
        self.logger.info("Loading catalogs from input directory.")
        extensions = ('.fits', '.fit', '.fts', '.asc',
                      '.txt', '.csv', '.dat', '.cat')
        catalog_files = [f for f in os.listdir(
            self.input_dir) if f.lower().endswith(extensions)]
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
                self.logger.error(
                    f"Failed to load {file} as FITS. Trying ASCII. Error: {e}")
                try:
                    catalogue = ascii.read(path).to_pandas()
                    catalogues.append(catalogue)
                    file_has_lodaded = True
                except Exception as e:
                    self.logger.error(
                        f"Failed to load {file} as ASCII. Skipping file. Error: {e}")
            if not file_has_lodaded:
                self.logger.critical(f"Could not load file: {file}")
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

        # if CLASS_STAR is present, filter out non-stellar objects
        if 'CLASS_STAR' in catalogues[0].columns:
            mask = catalogues[0]['CLASS_STAR'] >= 0.8
            catalogue = catalogues[0][mask]
        else:
            self.logger.info(
                "CLASS_STAR column not found. Skipping star/galaxy separation.")

        ref_coords = SkyCoord(
            ra=catalogues[0]['RA'].values * self.raunit,
            dec=catalogues[0]['DEC'].values * u.deg)

        for catalogue in catalogues[1:]:
            if 'CLASS_STAR' in catalogue.columns:
                mask = catalogue['CLASS_STAR'] >= 0.8
                catalogue = catalogue[mask]
            cat_coords = SkyCoord(
                ra=catalogue['RA'].values * self.raunit,
                dec=catalogue['DEC'].values * u.deg)
            idx, d2d, _ = cat_coords.match_to_catalog_sky(ref_coords)
            matched = d2d.arcsec < 1.0
            # Add unmatched coordinates to ref_coords
            unmatched_coords = cat_coords[~matched]
            ref_coords = SkyCoord(
                ra=np.concatenate(
                    [ref_coords.ra.deg, unmatched_coords.ra.deg]) * self.raunit,
                dec=np.concatenate(
                    [ref_coords.dec.deg, unmatched_coords.dec.deg]) * u.deg
            )

        self.unified_catalogue['RA'] = ref_coords.ra.deg
        self.unified_catalogue['DEC'] = ref_coords.dec.deg
        _array_size = ref_coords.ra.size
        # fill up the columns in the unified catalogue
        for i, catalogue in enumerate(catalogues):
            self.unified_catalogue[f'MAG_{i}'] = np.full(_array_size, np.nan)
            self.unified_catalogue[f'MAG_ERR_{
                i}'] = np.full(_array_size, np.nan)
            self.unified_catalogue[f'FLUX_{i}'] = np.full(_array_size, np.nan)
            self.unified_catalogue[f'FLUX_ERR_{
                i}'] = np.full(_array_size, np.nan)
            cat_coords = SkyCoord(
                ra=catalogue['RA'].values * self.raunit,
                dec=catalogue['DEC'].values * u.deg)
            idx, d2d, _ = cat_coords.match_to_catalog_sky(ref_coords)
            matched = d2d.arcsec < 1.0
            self.unified_catalogue[f'MAG_{i}'][idx[matched]
                                               ] = catalogue['MAG'].values[matched]
            self.unified_catalogue[f'MAG_ERR_{
                i}'][idx[matched]] = catalogue['MAG_ERR'].values[matched]
            self.unified_catalogue[f'FLUX_{i}'][idx[matched]
                                                ] = catalogue['FLUX'].values[matched]
            self.unified_catalogue[f'FLUX_ERR_{
                i}'][idx[matched]] = catalogue['FLUX_ERR'].values[matched]

        self.unified_catalogue = pd.DataFrame(self.unified_catalogue)
        self.logger.info("Data organization complete.")
        import pdb
        pdb.set_trace()

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

    def run(self):
        self.logger.info("Starting variability search process.")
        catalogues = self.load_catalogs()
        if not catalogues:
            self.logger.critical("No catalogs loaded. Exiting.")
            return
        self.organize_data(catalogues)
        import pdb
        pdb.set_trace()
        # self.define_reference_star()
        # self.select_comparison_stars()
        # self.search_variability()
        # self.logger.info("Variability search process completed.")


if __name__ == "__main__":
    args = parse_args()
    variability_search = VariabilityFinder(
        args,
        logger=logger(logfile=args.logfile,
                      loglevel=getattr(logging, args.loglevel))
    )
    variability_search.run()
