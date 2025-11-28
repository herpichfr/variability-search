[![Release](https://img.shields.io/github/v/release/herpichfr/variability-search)](https://github.com/herpichfr/variability-search/releases)
[![Issues](https://img.shields.io/github/issues/herpichfr/variability-search)](https://github.com/herpichfr/variability-search/issues)
[![License](https://img.shields.io/github/license/herpichfr/variability-search)](https://github.com/herpichfr/variability-search/blob/main/LICENSE)
[![Python](https://img.shields.io/badge/python-3.x-blue)](https://www.python.org/)
[![Stars](https://img.shields.io/github/stars/herpichfr/variability-search)](https://github.com/herpichfr/variability-search/stargazers)

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
