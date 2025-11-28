# variability-search

Project dedicated to search for variability in different types of astronomical observations.

## Description

This project analyzes photometry catalogues (CSV, ASCII, or FITS) from astronomical images to identify variable sources.

## Installation

### From source (development)

1. Clone the repository:
```bash
git clone https://github.com/herpichfr/variability-search.git
cd variability-search
```

2. Create and activate a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install in development mode:
```bash
pip install -e ".[dev]"
```

### From PyPI (when available)

```bash
pip install variability-search
```

## Usage

```bash
variability-search --help
variability-search catalog.fits -o output.csv
```

## Project Structure

```
variability-search/
├── src/
│   └── variability_search/    # Main package
├── data/                       # Data files (not tracked)
├── tests/                      # Test suite
├── pyproject.toml              # Project configuration
└── README.md
```

## Dependencies

- astropy
- numpy
- matplotlib

## Development

Run tests:
```bash
pytest
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
