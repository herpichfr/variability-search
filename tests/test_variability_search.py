"""Basic tests for variability_search package."""

import variability_search


def test_version():
    """Test that version is defined."""
    assert variability_search.__version__ == "0.1.0"


def test_variability_finder_imports():
    """Test that variability_finder module can be imported."""
    from variability_search import variability_finder

    # Check that the main function exists
    assert hasattr(variability_finder, 'main')
    assert hasattr(variability_finder, 'parse_args')
    assert hasattr(variability_finder, 'logger')


def test_prepare_catalogues_import():
    """Test that PrepareCatalogues can be imported."""
    from variability_search.prepare_catalogues import PrepareCatalogues

    assert PrepareCatalogues is not None


def test_variability_search_import():
    """Test that VariabilitySearch can be imported."""
    from variability_search.variability_search import VariabilitySearch

    assert VariabilitySearch is not None
