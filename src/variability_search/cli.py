"""Command-line interface for variability-search."""

import argparse

from variability_search import __version__


def main():
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
        Default=None,
        help="Path to working directory. Default is inputdir if not specified.",
    )
    parser.add_argument(
        "--outputdir",
        Default="output",
        help="Path to output directory. Default is 'output'.",
    )
    parser.add_argument(
        "--raunit",
        default="deg",
        help="Unit for Right Ascension (default: deg)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output file",
    )

    args = parser.parse_args()

    if args.catalog:
        print(f"Processing catalog: {args.catalog}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
