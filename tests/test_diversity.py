"""Tests for diversity metric calculations."""

import numpy as np
import pandas as pd
import pytest

from src.features.diversity import alpha_diversity, beta_diversity_matrix, pcoa


@pytest.fixture
def sample_abundances():
    """Create sample abundance data."""
    return pd.Series({
        "Species_A": 100,
        "Species_B": 50,
        "Species_C": 25,
        "Species_D": 10,
        "Species_E": 5,
    })


@pytest.fixture
def abundance_matrix():
    """Create multi-sample abundance matrix."""
    np.random.seed(42)
    return pd.DataFrame(
        np.random.dirichlet(np.ones(10), size=5) * 1000,
        index=[f"Sample_{i}" for i in range(5)],
        columns=[f"Species_{i}" for i in range(10)],
    )


class TestAlphaDiversity:
    def test_shannon_positive(self, sample_abundances):
        h = alpha_diversity(sample_abundances, metric="shannon")
        assert h > 0

    def test_shannon_max_for_even_community(self):
        even = pd.Series([100, 100, 100, 100])
        uneven = pd.Series([370, 10, 10, 10])
        assert alpha_diversity(even, "shannon") > alpha_diversity(uneven, "shannon")

    def test_simpson_between_0_and_1(self, sample_abundances):
        d = alpha_diversity(sample_abundances, metric="simpson")
        assert 0 < d < 1

    def test_chao1_geq_observed(self, sample_abundances):
        chao = alpha_diversity(sample_abundances, metric="chao1")
        observed = alpha_diversity(sample_abundances, metric="observed")
        assert chao >= observed

    def test_observed_counts_species(self, sample_abundances):
        assert alpha_diversity(sample_abundances, metric="observed") == 5

    def test_empty_returns_zero(self):
        empty = pd.Series(dtype=float)
        assert alpha_diversity(empty, metric="shannon") == 0.0

    def test_unknown_metric_raises(self, sample_abundances):
        with pytest.raises(ValueError, match="Unknown metric"):
            alpha_diversity(sample_abundances, metric="invalid")


class TestBetaDiversity:
    def test_bray_curtis_square_matrix(self, abundance_matrix):
        dm = beta_diversity_matrix(abundance_matrix, metric="bray_curtis")
        assert dm.shape == (5, 5)

    def test_diagonal_is_zero(self, abundance_matrix):
        dm = beta_diversity_matrix(abundance_matrix, metric="bray_curtis")
        np.testing.assert_array_almost_equal(np.diag(dm.values), 0)

    def test_symmetric(self, abundance_matrix):
        dm = beta_diversity_matrix(abundance_matrix, metric="bray_curtis")
        np.testing.assert_array_almost_equal(dm.values, dm.values.T)


class TestPCoA:
    def test_returns_correct_shape(self, abundance_matrix):
        dm = beta_diversity_matrix(abundance_matrix)
        coords = pcoa(dm, n_components=2)
        assert coords.shape == (5, 2)
        assert list(coords.columns) == ["PC1", "PC2"]
