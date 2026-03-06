"""Alpha and beta diversity metrics from taxonomic abundance profiles.

Computes ecological diversity measures that serve as features for
disease prediction models.

Example:
    >>> from src.features.diversity import alpha_diversity, beta_diversity_matrix
    >>> shannon = alpha_diversity(abundance_series, metric="shannon")
    >>> dm = beta_diversity_matrix(abundance_df, metric="bray_curtis")
"""

import logging
from typing import Literal

import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis, pdist, squareform

logger = logging.getLogger(__name__)


def alpha_diversity(
    abundances: pd.Series,
    metric: Literal["shannon", "simpson", "chao1", "observed"] = "shannon",
) -> float:
    """Calculate alpha diversity for a single sample.

    Args:
        abundances: Series of species abundances (counts or relative).
        metric: Diversity metric to compute.
            - shannon: Shannon entropy (H')
            - simpson: Simpson's diversity index (1 - D)
            - chao1: Chao1 richness estimator (requires counts)
            - observed: Number of observed species (richness)

    Returns:
        Diversity value as a float.

    Raises:
        ValueError: If metric is unknown or data is invalid.
    """
    # Remove zeros and NaN
    ab = abundances[abundances > 0].dropna()

    if ab.empty:
        logger.warning("Empty abundance data — returning 0.0")
        return 0.0

    if metric == "shannon":
        # Convert to relative abundances
        props = ab / ab.sum()
        return float(-np.sum(props * np.log(props)))

    elif metric == "simpson":
        props = ab / ab.sum()
        return float(1 - np.sum(props ** 2))

    elif metric == "chao1":
        # Chao1 requires integer counts
        counts = ab.astype(int)
        s_obs = (counts > 0).sum()
        f1 = (counts == 1).sum()  # singletons
        f2 = (counts == 2).sum()  # doubletons

        if f2 > 0:
            return float(s_obs + (f1 ** 2) / (2 * f2))
        elif f1 > 0:
            return float(s_obs + f1 * (f1 - 1) / 2)
        else:
            return float(s_obs)

    elif metric == "observed":
        return float((ab > 0).sum())

    else:
        raise ValueError(f"Unknown metric: {metric}. Use shannon, simpson, chao1, or observed.")


def alpha_diversity_table(
    abundance_matrix: pd.DataFrame,
    metrics: list[str] | None = None,
) -> pd.DataFrame:
    """Compute alpha diversity metrics for all samples.

    Args:
        abundance_matrix: DataFrame with samples as rows, species as columns.
        metrics: List of metrics to compute. Default: all four.

    Returns:
        DataFrame with samples as index and diversity metrics as columns.
    """
    if metrics is None:
        metrics = ["shannon", "simpson", "chao1", "observed"]

    results = {}
    for sample_id, row in abundance_matrix.iterrows():
        results[sample_id] = {
            metric: alpha_diversity(row, metric=metric)
            for metric in metrics
        }

    df = pd.DataFrame(results).T
    df.index.name = "sample_id"

    logger.info(f"Computed {len(metrics)} diversity metrics for {len(df)} samples")
    return df


def beta_diversity_matrix(
    abundance_matrix: pd.DataFrame,
    metric: Literal["bray_curtis", "jaccard"] = "bray_curtis",
) -> pd.DataFrame:
    """Compute pairwise beta diversity distance matrix.

    Args:
        abundance_matrix: DataFrame with samples as rows, species as columns.
        metric: Distance metric to use.
            - bray_curtis: Bray-Curtis dissimilarity
            - jaccard: Jaccard distance (presence/absence)

    Returns:
        Square DataFrame of pairwise distances (samples × samples).
    """
    if metric == "bray_curtis":
        distances = pdist(abundance_matrix.values, metric="braycurtis")
    elif metric == "jaccard":
        binary = (abundance_matrix > 0).astype(float)
        distances = pdist(binary.values, metric="jaccard")
    else:
        raise ValueError(f"Unknown metric: {metric}")

    dm = squareform(distances)
    dm_df = pd.DataFrame(
        dm,
        index=abundance_matrix.index,
        columns=abundance_matrix.index,
    )

    logger.info(f"Computed {metric} distance matrix: {dm_df.shape}")
    return dm_df


def pcoa(
    distance_matrix: pd.DataFrame,
    n_components: int = 3,
) -> pd.DataFrame:
    """Principal Coordinates Analysis on a distance matrix.

    Args:
        distance_matrix: Square pairwise distance matrix.
        n_components: Number of principal coordinates to return.

    Returns:
        DataFrame with PC coordinates for each sample.
    """
    # Classical MDS / PCoA
    n = len(distance_matrix)
    D = distance_matrix.values

    # Double centering
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(B)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Take top components (only positive eigenvalues)
    positive = eigenvalues > 0
    eigenvalues = eigenvalues[positive][:n_components]
    eigenvectors = eigenvectors[:, positive][:, :n_components]

    # Scale by sqrt of eigenvalues
    coords = eigenvectors * np.sqrt(eigenvalues)

    # Variance explained
    total_var = np.sum(eigenvalues[eigenvalues > 0])
    var_explained = eigenvalues / total_var

    result = pd.DataFrame(
        coords,
        index=distance_matrix.index,
        columns=[f"PC{i+1}" for i in range(n_components)],
    )

    logger.info(
        f"PCoA: {n_components} components explain "
        f"{sum(var_explained):.1%} of variance"
    )

    return result
