"""Command-line interface for variability-search."""

import argparse


def main():
    """Main entry point for the variability-search CLI."""
    parser = argparse.ArgumentParser(
        description="Search for variability in astronomical observations"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )
    parser.add_argument(
        "catalog",
        nargs="?",
        help="Path to input catalog file (CSV, ASCII, or FITS)",
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
