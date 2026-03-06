"""ML-based disease prediction from metagenomic features.

Trains and evaluates classifiers to predict disease status from
microbiome composition and functional features.

TODO:
    - [ ] Implement full training pipeline with nested CV
    - [ ] Add SHAP-based feature importance
    - [ ] Add model comparison report generation
    - [ ] Integrate with diversity features
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def train_disease_model(
    feature_matrix: pd.DataFrame,
    target_col: str = "disease_status",
    model_type: str = "random_forest",
    n_splits: int = 5,
    random_state: int = 42,
) -> dict:
    """Train a disease prediction model with cross-validation.

    Args:
        feature_matrix: DataFrame with features and target column.
        target_col: Name of binary target column.
        model_type: Model to train. Options: random_forest, gradient_boosting,
            logistic_regression.
        n_splits: Number of cross-validation folds.
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary with trained model, CV scores, and feature importances.
    """
    raise NotImplementedError(
        "ML training is planned for Phase 5. "
        "See README.md roadmap for current progress."
    )


def evaluate_model(
    model,
    X_test: pd.DataFrame,
    y_test: pd.Series,
) -> dict:
    """Evaluate a trained model on held-out test data.

    Args:
        model: Trained sklearn-compatible model.
        X_test: Test feature matrix.
        y_test: True labels.

    Returns:
        Dictionary with AUC, accuracy, precision, recall, F1, and
        confusion matrix.
    """
    raise NotImplementedError("Planned for Phase 5.")
